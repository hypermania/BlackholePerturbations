# High-precision SdS evolution optimization

## Goal

We benchmark the sourced Schwarzschild--de Sitter evolution used by the
areal-source scan. The production grid has 50,001 points, the scalar type is
IEEE binary128, and a fixed Dopri5 step has size `0.01`. The first pilot took
5435.98 seconds, or 90.6 minutes, for 100,000 steps on six CPU cores.

The acceptance condition is that the homogeneous and sourced right-hand sides
both sustain at least 90% of a matched minimum-kernel ceiling. We also require
the optimized full-profile evolution to agree exactly with the previous
implementation on the GCC 15.2 benchmark build. Since fast-math reassociation
is compiler dependent, builds using other GCC versions must agree within
binary128 roundoff.

## Performance model

The CPU implements binary128 arithmetic in software. A hardware peak-FLOP
bound is therefore not useful. We instead compare the production equation with
a compact interior kernel that has the same

- binary128 scalar backend,
- fourth-order stencil,
- potential and source streams,
- state layout, and
- static OpenMP partition.

The ceiling omits only equation-interface checks and the four boundary points.
For the sourced case, it includes the source callback and workspace fill. This
is an attainable upper bound for the selected discretization and generic source
interface.

The benchmark alternates production and ceiling calculations in adjacent
five-call blocks. This pairing is needed on the virtualized VPS: minute-scale
samples can otherwise differ by more than 10% because of host scheduling. The
benchmark uses seven samples and reports the median paired efficiency. Values
slightly above 100% represent sampling noise and are interpreted as saturation
of the ceiling.

## Bottleneck

The spatial operator was not the cause of the long run. It already uses the
compact symmetric stencil, restricted pointers, one OpenMP loop, and a reusable
source workspace. The translated source evaluation takes about 0.14 ms per
call in the accepted benchmark, which is also a small part of one Dopri5 step.

The remaining cost was the stage algebra. Dopri5 forms six full-state linear
combinations per step. Each combination had the form

```text
v0 = 1 * v1 + alpha2 * v2 + ...
```

The generic Eigen/Odeint adapter evaluated the multiplication by one for every
binary128 element. On the 100,002-component state, this performed six
unnecessary full-vector software multiplications per step.

## Changes

1. The binary128 `scale_sum1` through `scale_sum6` operations now test the
   first coefficient once per vector operation. When it is exactly one, they
   copy or add the first input directly. The fallback path is unchanged for
   all other scalar types and coefficients.
2. Under `-ffast-math`, GCC reassociates the old six-term expression, and the
   selected tree differs between GCC 11 and GCC 15. The unit path therefore
   uses explicit libgcc binary128 additions to impose a stable tree while
   omitting only `1 * v1`. It is bitwise identical to the old expression on
   the compiler used for the benchmark; portable cross-version checks require
   binary128 roundoff equivalence because the compiler-generated reference is
   itself version dependent.
3. A direct `-march=native` versus `-march=alderlake` check found the existing
   Alder Lake setting slightly faster for this software-binary128 workload,
   despite the virtual AMD CPU label. The original production setting was
   therefore retained rather than changed on architectural assumptions.
4. The SdS benchmark now uses the pilot parameters
   `q=0.1`, `s=0`, `l=0`, and `beta=0`, compiles with the production
   `-ffast-math` and architecture settings, and measures a paired complete
   Dopri5 step against a frozen pre-optimization operation policy. Wall and
   process CPU time are both reported.
5. The areal-scan runner disables dynamic OpenMP team sizing and records the
   maximum thread count in its metadata.

The PDE layout, source interface, stencil, boundary conditions, and integrator
are unchanged.

## Results

The accepted command is

```bash
make check-sds-precise-performance
```

on the six-core AMD EPYC VPS. The 50,001-point, 600-iteration result is

| Quantity | Result |
|---|---:|
| Homogeneous RHS | 20.59 million points/s |
| Homogeneous paired-ceiling efficiency | **100.5%** |
| Sourced RHS | 18.34 million points/s |
| Sourced paired-ceiling efficiency | **100.3%** |
| Source evaluations | 6543/s |
| Optimized sourced Dopri5 step | 45.40 ms |
| Frozen old sourced Dopri5 step | 53.26 ms |
| Paired wall-time speedup | **1.173 times** |
| Paired process-CPU speedup | **1.168 times** |
| Estimated 100,000-step time | 75.7 min |

The frozen-old and optimized steppers are alternated in the same benchmark
process. This avoids comparing separate VPS runs. The optimized stage algebra
reduces the measured step time by 14.8%. The original 90.6-minute pilot is
consistent with the old-step estimate of 88.8 minutes; the remaining difference
is normal VPS and observer overhead.

The wall-time estimate remains sensitive to VPS load. It is not a replacement
for recording the wall time of every production run. The current Teukolsky
benchmark also varies substantially with host scheduling. In the post-change
reference run, its sourced right-hand side reached 100.8% of its paired ceiling
and its sourced Dopri5 step took 62.74 ms, compared with 45.40 ms for SdS. The
earlier 30-minute Teukolsky run therefore does not indicate an SdS-specific
stencil regression on the current host and build.

## Numerical equivalence

The correctness test retains a direct implementation of the previous sourced
SdS right-hand side. It initializes full profiles with `std::mt19937_64` and
checks three cases:

| Case | `q` | `l` | `beta` |
|---:|---:|---:|---:|
| 1 | 0.1 | 0 | 0 |
| 2 | 0.4 | 1 | 1 |
| 3 | 0.8 | 2 | 2 |

For every case, the test compares the complete RHS and every component of the
full profile after each of 40 Dopri5 steps. The largest difference is exactly
zero. The production benchmark also gives

| Check | Maximum absolute difference |
|---|---:|
| One thread versus six threads | 0 |
| Production versus ceiling RHS | 0 |
| One-thread versus six-thread Dopri5 state | 0 |
| Optimized versus previous randomized full profile | 0 |

Normal, production, AddressSanitizer, and UndefinedBehaviorSanitizer tests pass
for both the SdS equation and the shared Eigen/Odeint stage algebra. The
production stage-algebra check passes with both GCC 11 and GCC 15: exact
unit-path identity is required without fast-math, while fast-math builds use a
roundoff-level comparison against the compiler-dependent legacy expression.
The Teukolsky correctness and performance checks also pass with the shared
change.

## Remaining limit

The optimized step now uses the minimum number of multiplications for the
existing Dopri5 linear combinations while imposing a stable binary128
addition order, and the spatial operator saturates its matched ceiling. This
order matches the former GCC 15.2 build exactly; a compiler-generated legacy
expression may differ by binary128 roundoff on other GCC versions. A
substantially shorter single-run time would require changing one of the current
constraints, such as the time integrator, scalar precision, grid resolution,
or number of CPU cores. Those changes would alter numerical semantics or
simulation parameters and are outside this optimization.

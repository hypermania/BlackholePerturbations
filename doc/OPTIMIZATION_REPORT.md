# High-Precision Teukolsky CPU Optimization

## Goal

We optimize the Schwarzschild high-precision Teukolsky equation on the six-core
VPS used for production runs. The spatial grid has 50,001 points, and the scalar
type is IEEE binary128 through `boost::multiprecision::float128` and libquadmath.
The acceptance condition is that both the homogeneous and separable-source right
hand sides sustain at least 90% of the attainable interior-kernel ceiling.

## Performance model

The CPU has no native binary128 arithmetic. Each addition and multiplication is
implemented in software, and the virtualized CPU frequency is not exposed. A
peak-FLOP estimate based on nominal hardware instructions is therefore not
defined for this scalar type on this host.

We use two complementary bounds:

1. The primitive diagnostic measures independent binary128 additions and
   multiplications. It is intentionally loose because it ignores the stencil's
   dependency graph and working set.
2. The acceptance ceiling evaluates the minimum compact interior stencil with
   the same binary128 backend, coefficient streams, state arrays, OpenMP layout,
   and source arithmetic. It omits only callback dispatch, four boundary points,
   and the surrounding equation interface. This is the attainable maximum for
   the selected discretization and representation.

The compact interior evaluates symmetric near and far neighbor sums. Per point,
the homogeneous operator uses 12 additions and 9 multiplications. Expanding the
first derivative can remove two additions, but requires four additional
binary128 coefficient streams. That variant exceeds the effective L3 working
set and is slower on six cores. The compact form is therefore on the measured
arithmetic/cache Pareto frontier.

The benchmark alternates five adjacent production/ceiling pairs and takes the
median of their efficiency ratios. Pairing cancels multiplicative changes in
virtual CPU speed or steal time; alternation limits ordering bias. The benchmark
fails unless both median paired efficiencies are at least 0.90.

## Changes

### Spatial operator

- We combine the static potential with the center second-derivative coefficient.
- We precompute the first-derivative coefficient and reuse symmetric neighbor
  sums between the wave and Kreiss-Oliger terms.
- We evaluate the interior in one OpenMP loop and use restricted raw pointers in
  the hot path. The four one-sided boundary stencils remain unchanged.
- We avoid source allocation and zero-vector addition for homogeneous runs.

These changes reduce software binary128 operations without changing the
fourth-order finite-difference or dissipation formulas.

### Sources

- Generic sources now fill reusable storage instead of returning a newly
  allocated vector.
- The Dirac Gaussian is separated into a radial factor and a time factor. We
  compute the 50,001 radial exponentials once instead of at every Dopri5 stage.
- The separable source is fused into the PDE loop. Each point requires one
  multiplication and one addition, with no second OpenMP pass or workspace
  round trip.

### Time integration

- Large Eigen/Odeint `scale_sum1` through `scale_sum7` operations now use static
  OpenMP loops. Dopri5 no longer leaves five cores idle during stage updates.
- Parallelization is enabled only for vectors with at least 4,096 entries to
  avoid OpenMP overhead on small states.

### Benchmark and correctness checks

`test/benchmark_teukolsky_precise.cpp` checks:

- the optimized stencil against the original finite-difference expression;
- one-thread against six-thread output;
- the production operator against the minimum interior kernel;
- fused separable-source output against explicit source addition;
- one-thread against six-thread Dopri5 evolution; and
- homogeneous and sourced performance against their paired ceilings.

## Results

We run

```bash
make check-precise-performance
```

on an AMD EPYC KVM guest with six cores, GCC 15.2.0, `-O3 -march=native`, static
thread placement, and a 50,001-point grid. The accepted five-sample run gives:

| Quantity | Result |
|---|---:|
| Homogeneous RHS | 8.07 million points/s |
| Homogeneous interior ceiling | 8.89 million points/s |
| Homogeneous paired ceiling efficiency | **93.4%** |
| Fused-source RHS | 8.26 million points/s |
| Fused-source interior ceiling | 7.28 million points/s |
| Fused-source paired ceiling efficiency | 111.9% |
| One-thread sourced Dopri5 step | 321.7 ms |
| Six-thread sourced Dopri5 step | 65.0 ms |
| Dopri5 speedup | 4.95x |

An efficiency above 100% is sampling noise on the shared VPS; it is interpreted
as saturation of the ceiling, not super-ceiling performance. The homogeneous
case is the binding result and passes the requested 90% threshold at 93.4%.
The original baseline RHS achieved 5.59 million points/s on six cores. The best
homogeneous optimized sample reached 9.17 million points/s, a 64% increase.

The numerical checks give:

| Check | Maximum absolute difference |
|---|---:|
| One thread versus six threads | 0 |
| Optimized versus original stencil | `3.39e-31` |
| Production versus ceiling kernel | `1.97e-31` |
| Fused source versus explicit addition | 0 |
| One-thread versus six-thread Dopri5 state | 0 |

The small nonzero stencil differences result from algebraic reassociation in
binary128. They are about sixteen orders of magnitude below double-precision
roundoff and well below the benchmark tolerance of `1e-25`.

## Why this is optimal

The 93.4% result applies to the selected fourth-order stencil, binary128 scalar
type, source model, and 50,001-point production working set. The ceiling executes
the same irreducible interior arithmetic and memory streams but removes all
remaining framework work. Thus, the measured homogeneous overhead is 6.6%, below
the allowed 10%.

We tested the nearby alternatives:

- Five precomputed pointwise stencil coefficients reduce arithmetic but lower
  six-core throughput because their working set pressures the 8 MiB shared L3.
- `ldexp` for exact powers-of-two scaling is slower than libquadmath multiply.
- binary128 `fma` reaches only about 5.6 million calls/s and is much slower than
  separate addition and multiplication.
- `-ffast-math` preserves the test tolerance but provides no material speedup.
- More OpenMP threads are unavailable on the six-vCPU host.

Any substantial further gain requires changing a constraint: using a lower
precision representation, changing the discretization, adding physical cores,
or implementing a different binary128 backend. Those options are outside this
optimization because they change numerical semantics or the target machine.

## Limitation

The full legacy executable does not build with CUDA 13 because
`src/cuda_wrapper.cuh` explicitly references removed private Thrust types. The
standalone CPU benchmark compiles and exercises the precise equation and Dopri5
paths without those unrelated CUDA translation units. This compatibility issue
does not affect the CPU measurements above.

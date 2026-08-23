# Precise sourced Schwarzschild-de Sitter solver

## Scope

`src/sds_precise.hpp` evolves

```text
partial_t psi = Pi,
partial_t Pi  = partial_x^2 psi - V_s_l psi + S(t,x),
```

where

```text
V_s_l = f(r) [l(l+1)/r^2 + (1-s^2) 2M/r^3],
f(r)  = 1 - 2M/r - Lambda r^2/3.
```

The implemented values `s = 0, 1, 2` correspond to the conformal scalar,
electromagnetic, and axial gravitational sectors. The multipole is read from
`SdSMasterPDEPreciseParam` and the potential is built once in the constructor.
Even-parity Zerilli perturbations and the general massive or nonminimally
coupled scalar equation are not included.

The spatial operator uses the existing fourth-order Regge-Wheeler stencil and
absorbing outer closures. It does not add Kreiss-Oliger dissipation.

## Geometry preprocessing

The coordinate convention is

```text
x(3M) = 3M + 2M log(1/2),
```

which tends to the standard Schwarzschild tortoise coordinate as `Lambda` tends
to zero. The two positive horizons are found from

```text
Lambda r^3 - 3r + 6M = 0
```

with bracketed TOMS 748 iterations in 100-decimal arithmetic. The negative root
is then set to `r_o = -(r_b+r_c)`, which enforces the exact Vieta relation. An
independent bisection test compares these values with the analytic trigonometric
roots for `Lambda` between `1e-100` and a near-Nariai value. The largest relative
difference is `8.89e-51`, well below binary128 roundoff. Root finding is retained
because the analytic expression for `r_b` progressively loses digits through
cancellation as `Lambda` tends to zero.

The inversion uses

```text
y = log[(r-r_b)/(r_c-r)]
```

and a bracketed Halley iteration. This variable maps the finite static interval
`r_b < r < r_c` to the complete real line. In particular,

```text
r-r_b = (r_c-r_b) exp(y)/(1+exp(y)),
r_c-r = (r_c-r_b)/(1+exp(y)),
dr/dy = (r-r_b)(r_c-r)/(r_c-r_b).
```

Since

```text
f = Lambda (r-r_b) (r_c-r) (r-r_o)/(3r),
dx/dr = 1/f,
```

the small horizon distances cancel from the transformed derivative:

```text
dx/dy = 3r/[Lambda (r_c-r_b) (r-r_o)].
```

Thus `x(y)` is monotonic and has a finite, well-conditioned derivative at both
horizons. Stable logistic branches return `r-r_b` and `r_c-r` directly. The
metric function is then formed as

```text
f = Lambda (r-r_b) (r_c-r) (r-r_o) / (3r)
```

before conversion to binary128. This avoids subtracting rounded order-one
terms near either horizon.

The additive coordinate convention is represented by one precomputed
`tortoise_constant`; the old reference radius and three reference distances are
not stored separately. The cancellation involved in collapsing these terms is
included in a constructor-time roundoff estimate. Halley iteration normally
uses a relative residual threshold of about `1e-80`; only extremely small
`Lambda` values use the larger calculated 100-decimal roundoff floor.

OpenMP workers receive contiguous grid blocks. The first point of each block is
bracketed from `y=0`. Each subsequent point uses the preceding converged point
and the predictor

```text
y_i^(0) = y_(i-1) + [x_i-x(y_(i-1))]/[dx/dy]_(i-1).
```

The predicted point and the preceding point form the initial monotonic bracket;
the upper endpoint is expanded only if the predictor undershoots. This retains
parallel initialization while reusing nearby inversions. `midpoint_x` and the
surface gravities are no longer part of the private inversion geometry.

Production preprocessing uses `cpp_bin_float_100`. The correctness test uses a
separate `cpp_bin_float<2000>` calculation (2000 decimal digits) at a
cancellation-dominated point. The production values of `f` and `V`
agree with that reference at binary128 accuracy. Direct binary128 evaluation of
the unfactored `f` is wrong by many orders of magnitude at the same point.

For `M=0.5`, the tests also evaluate the Schwarzschild expression directly,

```text
x = r + log(r-1),
```

at fixed grid coordinates. Its discrepancy decreases as `Lambda` is reduced
from `1e-8` through `1e-16`, verifying the selected coordinate convention.

## Sources

The optimized translated source has the form

```text
S(t,x) = Theta(t-t_on) F(t-x) a(x).
```

Four precomputed profiles are available:

1. `r^(-beta)`;
2. `r^(-beta) - r_c^(-beta)`;
3. `f r^(-beta)`;
4. `chi_+(x) [L/(x+x0)]^beta`.

The subtracted profile uses the float128 backend's `log1pq` and `expm1q`
implementations near the cosmological horizon. The vendored Boost overloads
delegate to these same functions but do not compile with GCC 15 because their
`__float128` return conversion is implicit; explicit backend construction avoids
that compatibility defect without maintaining local series implementations.
The tortoise profile uses a smooth compact transition between `X0` and `X1`.
Waveforms can be a normalized Gaussian or its derivative; the latter has zero
retarded-time mean.

The persistent geometry vectors are now only `r`, `rho_cosmological`, `f`, and
the physical potential `V`. The uniform tortoise grid is evaluated through one
shared `grid_coordinate` helper. The evolution uses `V` directly and forms the
constant fourth-order center contribution inside `operator()`; it does not
retain a second transformed-potential vector.

At each Runge-Kutta stage the solver computes the finite grid interval inside
the configured Gaussian cutoff. Only this interval evaluates binary128
exponentials, and source values are fused into the stencil. Homogeneous,
translated-source, and generic-source kernels are selected once per RHS call so
there is no source branch at every homogeneous grid point.

## Output and entry point

`run_sds_precise_eqn()` in `src/examples.cpp` contains one production
configuration block. Run it with:

```bash
./main --sds-precise
```

Fixed-position values and snapshots are converted to `double` by the existing
observers. `source_and_geometry.txt` records the source parameters, horizons,
surface gravities, coordinate convention, and the conservative finite-domain
fitting limit

```text
t_fit < 2 X_max - x_obs + u_min.
```

## Verification

`make check-sds-precise-correctness` checks normal and production builds. The
matrix covers invalid parameters, analytic and independently root-found cubic
roots, small cosmological constant, an extreme `Lambda=1e-100` geometry, a
near-Nariai case, the Schwarzschild limit, ordered horizons, factorized
coefficients, the 2000-decimal reference, constructor multipoles, all source
profiles, source turn-on, the zero-mean waveform, independent
original-formula boundary and interior stencils, one/two/six-thread
reproducibility, one- versus six-thread geometry construction, and 25 complete
Dopri5 steps.

`make check-sds-precise-sanitizers` passes under AddressSanitizer and
UndefinedBehaviorSanitizer. The expensive 2000-decimal reference is omitted
from the sanitizer build; it is exercised in both normal and production
correctness builds.

`make check-sds-precise-performance` measures geometry initialization,
homogeneous and sourced RHS throughput, complete Dopri5 steps, and paired
minimal-kernel ceilings on the 50,001-point grid. The optimized homogeneous and
sourced kernels both meet the 90% paired-ceiling gate. Performance samples can
vary on the shared KVM host, so the benchmark alternates production and ceiling
measurements and reports their medians.

The old independent and new neighbor-predicted constructors were measured with
the same 50,001-point, 60-iteration, six-thread benchmark on this VPS:

| Geometry inversion | Initialization time |
|---|---:|
| Independent pointwise inversion | 47.35 s |
| Contiguous blocks with neighbor predictor | 10.01 s |

The neighbor-predicted implementation is 4.73 times faster. In the same paired
run, retaining physical `V` instead of `center_factor` reduced homogeneous RHS
throughput from 26.8 to 23.7 million point updates per second and sourced
throughput from 20.0 to 18.7 million point updates per second. These figures are
host-sensitive, but they quantify the clarity/performance tradeoff requested
for the persistent equation state. The maintained 180-iteration acceptance run
still passes both 90% paired-ceiling gates with exact threaded results.

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

and Boost.Math's bracketed `halley_iterate`. This variable maps the finite
static interval `r_b < r < r_c` to the complete real line. In particular,

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

OpenMP workers receive contiguous grid blocks. For `P=N+1` points and `T`
workers, worker `t` owns the half-open interval

```text
[floor(t P/T), floor((t+1) P/T)).
```

The last worker therefore has `end=P=N+1`, but the strict `i<end` loop bound
makes `N` the largest possible index. If `T>P`, some intervals are empty; a
16-worker/five-point test covers this case. Geometry arrays are written only at
the owning worker's indices. The only shared result is the first failed index,
which uses an OpenMP `min` reduction instead of atomic storage.

The first point of each nonempty block is anchored at `y=0`. Each subsequent
point uses the preceding converged point and the predictor

```text
y_i^(0) = y_(i-1) + [x_i-x(y_(i-1))]/[dx/dy]_(i-1).
```

The predicted point and the anchor form the initial bracket. If they remain on
the same side of the root, the predicted step is doubled until the residual
changes sign. This single sign-based procedure works for both directions and
replaces the separate left/right and first/neighbor branch trees. The bracket,
predictor, and the tuple `(x-target, dx/dy, d2x/dy2)` are then passed to
`boost::math::tools::halley_iterate`; no local Halley update is maintained.

The residual `x(y)-target` still has two necessary roles. Its sign establishes
the initial bracket, and its magnitude validates the result against the
physical tortoise-coordinate tolerance. The Boost callback reports an exact
zero once this residual is below that tolerance, preventing unnecessary
iterations to the full storage precision. The returned point is independently
checked against the same threshold before it is accepted.

This retains parallel initialization while reusing nearby inversions.
`midpoint_x` and the surface gravities are no longer part of the private
inversion geometry.

`grid_space()` is the sole spacing calculation. Its binary128 result is reused
by the evolution and promoted to the preprocessing type, while every coordinate
is formed through the shared `grid_coordinate()` helper.

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

`SdSSource` owns source validation, spatial-profile preprocessing, numerical
cutoffs, and evaluation. A source is passed to the PDE constructor:

```cpp
SdSTranslatedSourceParam source_param;
SdSMasterPDEPrecise equation(param, SdSSource(source_param));
```

An empty `SdSSource` is used by the default homogeneous constructor. A custom
callback can also be wrapped in `SdSSource`, so the PDE treats built-in and
user-defined sources through the same callable interface.

The translated source has the form

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
retarded-time mean. `beta` may be zero or negative. In particular, `beta=0`
gives profiles `1`, `0`, `f`, and `chi_+`, respectively.

`SdSSpacetimeGaussianSourceParam` provides the localized Green-function probe

```text
S(t,x) = A exp(-[(x-x0)^2+(t-t0)^2]/sigma^2)
           / [sqrt(2 pi) sigma].
```

The spatial Gaussian is precomputed, while its scalar time factor is evaluated
at each Runge-Kutta stage. Values beyond the configurable `cutoff_sigma` are
set to zero.

The persistent geometry vectors are now only `r`, `rho_cosmological`, `f`, and
the physical potential `V`. The uniform tortoise grid is evaluated through one
shared `grid_coordinate` helper. The evolution uses `V` directly and forms the
constant fourth-order center contribution inside `operator()`; it does not
retain a second transformed-potential vector.

At each Runge-Kutta stage `SdSSource` fills one reusable workspace. It computes
only the finite grid interval inside the configured Gaussian cutoff. The PDE
operator contains no source-specific formulas or window logic: it calls `Q`
when present, then uses one unified assignment

```text
dPi = spatial_stencil - V psi + Q_workspace
```

for the interior and all four boundary points. The homogeneous case uses the
same path with a workspace initialized once to zero.

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
profiles including `beta=0`, source turn-on, the zero-mean waveform, the
spacetime Gaussian formula, independent
original-formula boundary and interior stencils, one/two/six-thread
reproducibility, one- versus six-thread geometry construction, and 25 complete
Dopri5 steps.

`make check-sds-precise-sanitizers` passes under AddressSanitizer and
UndefinedBehaviorSanitizer. The expensive 2000-decimal reference is omitted
from the sanitizer build; it is exercised in both normal and production
correctness builds.

`make check-sds-precise-performance` measures geometry initialization,
homogeneous and sourced RHS throughput, complete Dopri5 steps, and paired
minimal-kernel ceilings on the 50,001-point grid. The unified kernel meets the
90% paired-ceiling gate with both a zero workspace and a translated source.
Performance samples can vary on the shared KVM host, so the benchmark alternates
production and ceiling measurements and reports their medians.

The inversion variants were measured with the same 50,001-point, six-thread
configuration on this VPS:

| Geometry inversion | Initialization time |
|---|---:|
| Independent pointwise inversion | 47.35 s |
| Contiguous blocks, local safeguarded Halley | 10.01 s |
| Contiguous blocks, Boost.Math Halley | 13.77 s median |

The maintained Boost implementation is 3.44 times faster than independent
inversion, although Boost's more general safeguards make it 38% slower than the
former specialized local loop. A Boost Newton sample took 16.41 s, so Halley
was retained. In the same paired run, retaining physical `V` instead of
`center_factor` reduced homogeneous RHS
throughput from 26.8 to 23.7 million point updates per second and sourced
throughput from 20.0 to 18.7 million point updates per second. These figures are
host-sensitive, but they quantify the clarity/performance tradeoff requested
for the persistent equation state. The maintained 180-iteration acceptance run
still passes both 90% paired-ceiling gates with exact threaded results.

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
to zero. The three horizons are evaluated analytically in 100-decimal
arithmetic. The inversion uses

```text
y = log[(r-r_b)/(r_c-r)]
```

and a bracketed Halley iteration. Stable logistic branches return `r-r_b` and
`r_c-r` directly. The metric function is then formed as

```text
f = Lambda (r-r_b) (r_c-r) (r-r_o) / (3r)
```

before conversion to binary128. This avoids subtracting rounded order-one
terms near either horizon. On the production grid the maximum inversion
residual is approximately `1e-80` in tortoise-coordinate units.

Production preprocessing uses `cpp_bin_float_100`. The correctness test uses a
separate `cpp_bin_float<2000>` calculation at a cancellation-dominated point.
The production values of `f` and `V` agree with that 2000-decimal reference at
binary128 accuracy. Direct binary128 evaluation of the unfactored `f` is wrong
by many orders of magnitude at the same point.

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

The subtracted profile uses cancellation-free `log1p` and `expm1` series near
the cosmological horizon. The tortoise profile uses a smooth compact transition
between `X0` and `X1`. Waveforms can be a normalized Gaussian or its derivative;
the latter has zero retarded-time mean.

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
matrix covers invalid parameters, small cosmological constant, a near-Nariai
case, the Schwarzschild limit, ordered horizons, inversion residuals,
factorized coefficients, the 2000-decimal reference, constructor multipoles,
all source profiles, source turn-on, the zero-mean waveform, independent
original-formula boundary and interior stencils, one/two/six-thread
reproducibility, and 25 complete Dopri5 steps.

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

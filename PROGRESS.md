# Progress log

## 2026-08-25: Continuous SdS q scan interface

Commit: `5267a623`

### Problems encountered

- The areal-scan runner accepted only four hardcoded values of
  \(q=9\Lambda M^2\), preventing intermediate spacetime parameters.
- Existing output directories use decimal-point-free labels such as `q_01`,
  so accepting ordinary decimal text would have changed their naming scheme.

### Solutions

- Replaced the value lookup with decimal-code conversion: the CLI token `0123`
  is converted exactly to binary128 `0.123` and retained as the directory
  label `q_0123`.
- Left the physical \(0<q<1\) constraint to the existing centralized SdS
  geometry validation. Added exact conversion tests and a compiled CLI smoke
  check for an intermediate value.

### How to avoid these problems

- Treat a finite experiment matrix as a sampling plan rather than an API
  restriction unless the numerical equation itself imposes the bound.
- Separate lossless input encoding from physical parameter validation, and
  test both the encoded text and the resulting numerical value.

## 2026-08-25: Unrestricted SdS scan parameters

Commit: `2b3df3ef`

### Problems encountered

- The scan runner duplicated the SdS parameter checks from the equation class
  and restricted `s`, `l`, and `beta` to the original finite scan matrix.
- The CLI parsed `beta` as an integer, so fractional powers could not reach the
  existing binary128 source parameter.
- Removing the integer bounds exposed signed-overflow risk in `s*s` and
  `l*(l+1)` before those coefficients were converted to high precision.

### Solutions

- Kept parameter validation in `sds_precise.hpp`: `s` must be nonnegative,
  `beta` must be finite, and `l` has no range or ordering restriction.
- Parsed `beta` directly as binary128, retained a round-trip-safe value in
  output directory names, and constructed the validated equation before
  creating output files.
- Converted `s` and `l` to the high-precision type before evaluating their
  potential coefficients. Added normal, production, and sanitizer coverage
  for fractional and negative beta, `s > 2`, unrestricted `l`, and extreme
  64-bit coefficient inputs.

### How to avoid these problems

- Put reusable numerical-domain validation at the equation or source boundary,
  not in individual runners.
- Parse continuous numerical parameters directly into their computational
  type and test values outside the initial experiment matrix.
- Convert unrestricted integer parameters before arithmetic when the
  destination type has a wider numerical range.

## 2026-08-24: SdS scan interface and fit documentation cleanup

Commit: `ece9bcf`

### Problems encountered

- The areal-scan runner hardcoded `s=0`, saved three observer locations when
  only `x=50` was required, and exposed two one-use conversion and snapshot
  helpers outside the runner.
- The completed pilot uses the old six-column observer layout, while future
  runs use the simpler two-column layout. Removing the old layout outright
  would make the saved pilot difficult to reproduce.
- The large Jacobian column for the fitted constant `a` could be mistaken for
  evidence of a resolved nonzero exponent. For this model, that column is
  identically one and its norm is fixed mainly by the number of samples.
- Numerical reports accumulated in the repository root and obscured the
  source, script, and test layout.

### Solutions

- Added `s` to the scan API and CLI, reduced fixed-position output to
  `(psi_x50, Pi_x50)`, and moved the snapshot schedule and binary128-to-double
  conversions into `run_sds_areal_scan`.
- Updated the analysis to consume the new two-column layout while selecting
  the `x=50` columns from the already completed six-column pilot. Added tests
  for both layouts.
- Documented the full pilot Jacobian, column norms, singular values, projected
  `a` sensitivity, and window-stability limitation. The supported conclusion
  remains `a=0`, not an observer-dependent value near `1e-20`.
- Moved the committed numerical reports into `doc/` and updated repository and
  project references.

### How to avoid these problems

- Keep runner output schemas as small and explicit as the planned analysis
  requires, and test any intentional compatibility path for existing data.
- Assess a nonlinear fit with the complete Jacobian, singular values, nuisance
  directions, and window stability. A large norm for one column is not an
  identifiability test.
- Add new numerical reports under `doc/` rather than the repository root.

## 2026-08-23: Precise SdS evolution optimization

Commit: `c830b75`

### Problems encountered

- The 50,001-point areal-source pilot took 5435.98 seconds, even though the
  spatial SdS operator was already at its paired attainable ceiling. Profiling
  showed that Dopri5 repeatedly multiplied the full binary128 base state by an
  exact coefficient of one in six stage combinations.
- A direct removal of those multiplications was not bitwise identical under
  the production `-ffast-math` flags because GCC reassociated the six-term
  expression. Independent long wall-clock samples also varied noticeably on
  the virtualized host.
- The VPS reports an AMD virtual CPU, but `-march=native` was slightly slower
  than the established `-march=alderlake` setting for this software-binary128
  workload.

### Solutions

- Added binary128-only unit-coefficient paths to the shared Odeint/Eigen stage
  operations. The six-input path reproduces the former libgcc binary128
  addition tree exactly while omitting only multiplication by one.
- Added a frozen pre-optimization operations policy and seeded randomized
  full-profile tests over three `(q, l, beta)` cases. The complete right-hand
  side and every component after each of 40 Dopri5 steps are bitwise identical.
- Paired production and minimum-kernel ceiling calls in five-call blocks and
  compared old/new Dopri5 steps adjacently using wall and process CPU time. The
  homogeneous and sourced RHS reach 100.5% and 100.3% of their paired ceilings;
  the sourced step improves from 53.26 ms to 45.40 ms, or 1.173 times.
- Retained `-march=alderlake` after the direct compiler-target comparison.

### How to avoid these problems

- Profile the complete time step, not only the PDE right-hand side, before
  attributing a long run to source or stencil evaluation.
- Compare optimized and frozen expressions under the exact production flags;
  algebraic equivalence does not guarantee bitwise identity with fast-math.
- Use adjacent paired samples and process CPU time on this shared VPS, and
  benchmark compiler targets instead of inferring the best target from the
  virtual CPU label.

## 2026-08-23: SdS areal-source scan pilot

Commit: `a83d9ce`

### Problems encountered

- The requested three-parameter model
  \(p_{\rm loc}=a+b/(t-c)\) becomes rank deficient when the solution has
  already reached a memory plateau. A generic nonlinear least-squares solve
  can then leave \(c\) near its starting value and report misleadingly small
  covariance errors.
- The first production run took 5435.98 seconds on six cores. Extrapolating
  this cost to all 48 parameter combinations gives about 72.5 hours of
  sequential wall time.

### Solutions

- Profiled the fit over \(c\): for every trial \(c\), solve \(a,b\) exactly
  by linear least squares, then perform a one-dimensional minimization. The
  analysis now records the Jacobian rank, condition number, boundary hits,
  and three fitting-window starts.
- Added an exponential diagnostic to distinguish a true algebraic tail from
  pole-controlled relaxation to memory. For `q_01_l_0_beta_0`, the defensible
  result is \(a=0\), while the approach to memory has
  \(\gamma\simeq0.1021\), close to \(\kappa_c\).
- Stopped after the single requested pilot. The full scan remains pending a
  scheduling benchmark or source-evaluation optimization.

### How to avoid these problems

- Never interpret the formal covariance error of the \((a,b,c)\) fit without
  checking rank, condition number, observer agreement, and fitting-window
  stability.
- Do not launch the full matrix by multiplying the single-run command until
  the measured 90.6-minute cost has been addressed and the pilot analysis has
  been reviewed.

## 2026-08-24: Restore missing Boost.Typeof compatibility header

Commit: `b21a18e7`

### Problem encountered

- The vendored Boost 1.84 subset included `boost/typeof/typeof.hpp`, which
  expands an include of `boost/typeof/incr_registration_group.hpp`, but did
  not include that compatibility header. Builds succeeded on machines where
  GCC silently found a system Boost copy and failed on machines without it.

### Solution

- Added the missing five-line compatibility header byte-for-byte from the
  official Boost 1.84 release. Dependency generation now resolves the include
  from `external/boost`, and the normal and production precise-solver suites
  pass.

### How to avoid this problem

- When updating the `bcp`-generated Boost subset, audit macro-generated
  includes as well as literal includes and test from a machine without a
  usable system Boost fallback.

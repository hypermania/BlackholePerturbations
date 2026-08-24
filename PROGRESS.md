# Progress log

## 2026-08-24: Complete self-contained Boost 1.84 dependency

Commits: `ba03321`, `f1095fc`

### Problems encountered

- The repository's fixed-version Boost directory was a 3,275-file
  `bcp --scan` subset rather than the complete Boost 1.84 public-header tree.
  `bcp` did not discover the macro-expanded include
  `boost/typeof/incr_registration_group.hpp`.
- Builds on this VPS succeeded accidentally by loading that missing header
  from system Boost 1.74, mixing two Boost releases. A machine without the
  system header failed while compiling `examples.cpp` and `main.cpp`.
- The first Ubuntu 22.04 CI run used its default GCC 11, which rejects an
  existing template specialization in `cubic_scalar.hpp`; the documented
  project compiler is GCC 12.2 or newer.

### Solutions

- Replaced the subset with all 15,689 public-header files from the official
  Boost 1.84.0 archive after verifying its published SHA-256. The previously
  committed 3,275 files are byte-identical to the official release.
- Added archive provenance, the official license, and a deterministic verifier
  for version, file count, byte count, complete tree digest, and compiler
  include closure. The current dependency graph resolves 1,789 Boost headers,
  all under `external/boost`.
- Added a CPU-build CI workflow and confirmed forced CPU, normal, production,
  and sanitizer builds with no numerical regression.
- Made CI install and select GCC 12 explicitly instead of depending on the
  runner's older default compiler.

### How to avoid these problems

- Do not use an unqualified `bcp --scan` result as a complete dependency:
  macro-generated includes are not reliably visible to its scanner.
- Verify vendored dependencies from a checksum-verified official archive and
  audit the compiler's resolved include graph so a system-library fallback
  cannot hide an incomplete tree.
- Update Boost separately from numerical code and regenerate the recorded
  provenance and tree fingerprint for every version change.

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

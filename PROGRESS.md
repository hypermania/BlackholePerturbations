# Progress log

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

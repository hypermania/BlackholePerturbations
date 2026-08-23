# Schwarzschild-de Sitter areal-source scan pilot

This note records the scan interface, analysis procedure, and the first
production pilot. The full 48-run matrix has not been started.

## Scan definition

The planned matrix is

\[
\beta\in\{0,1,2\},\qquad
\ell\in\{0,1,2,3\},\qquad
q\equiv9\Lambda M^2\in\{0.1,0.2,0.4,0.8\}.
\]

Every run fixes \(s=0\) and \(M=0.5\), and uses

\[
S(t,x)=F(t-x)r(x)^{-\beta},\qquad
F(u)=\frac{A}{\sqrt{2\pi}\sigma}
\exp\left[-\frac{(u-u_0)^2}{2\sigma^2}\right],
\]

with \(A=1\), \(u_0=-10\), \(\sigma=0.5\), and an effective cutoff of
\(12\sigma\). Thus, \(A\) is the integral, or nonzero mean, of the Gaussian.

The production grid is

\[
x\in[-500,1000],\qquad \Delta x\simeq0.03,\qquad
t\in[0,1000],\qquad \Delta t=0.01.
\]

Time series are saved at the nearest grid points to \(x=0,50,100\). Full
\((\psi,\Pi)\) snapshots are saved at 27 requested times between \(t=0\) and
\(t=1000\).

One parameter triple is run with

```bash
./main --sds-areal-scan Q L BETA
```

and is stored in

```text
output/sds_areal_scan/q_XX_l_L_beta_B/
```

The directory contains the parameter and geometry metadata, the three-observer
time series, all snapshots, the final state, derived instantaneous slopes,
fit results, and plots. Simulation arrays are binary double-precision files.
They are generated data and remain excluded from Git.

## Tail extraction

The directly evolved momentum gives the signed instantaneous power

\[
p_{\rm loc}(t)=\frac{d\ln|\psi|}{d\ln t}=\frac{t\Pi}{\psi}.
\]

The requested model is

\[
p_{\rm loc}(t)=a+\frac{b}{t-c}.
\]

For each fixed \(c\), the analysis solves for \(a\) and \(b\) by linear least
squares. It then minimizes the profiled residual over \(c\), constrained to
remain to the left of the fitting interval so that the model has no pole in
the data. The reported tail power is the signed value \(a\), so
\(\psi\sim t^a\). The nominal window is \(t\in[500,900]\). Fits beginning at
\(t=400,500,600\) are all saved to diagnose window dependence.

The analysis also records the Jacobian rank and condition number. A small
formal covariance error is not accepted as evidence for a tail exponent when
the three-parameter fit is rank deficient.

## Pilot: `q_01_l_0_beta_0`

The pilot used

\[
q=0.1,\qquad \Lambda=0.0444444444444444\ldots,\qquad
\ell=0,\qquad \beta=0.
\]

It completed 100,000 steps in 5435.98 seconds (90.6 minutes) on six CPU cores.
The output contains exactly 100,001 finite time samples and 27 correctly sized
snapshots. The conservative right-boundary contamination limit at the
rightmost observer is \(t<1884\), safely beyond the simulated interval.

The field approaches a nonzero memory:

| Observer | Actual \(x\) | Memory estimate |
|---:|---:|---:|
| 0 | 0.0050001 | 0.7046799071 |
| 50 | 49.9859997 | 23.9612618333 |
| 100 | 99.9969999 | 48.9666025270 |

The nominal fitted values are

| Observer | \(a\) | \(b\) | \(c\) | Jacobian rank | Condition number |
|---:|---:|---:|---:|---:|---:|
| 0 | \(9.89\times10^{-24}\) | \(-1.57\times10^{-21}\) | 496.523 | 2 | \(3.32\times10^{23}\) |
| 50 | \(1.53\times10^{-22}\) | \(-2.43\times10^{-20}\) | 496.523 | 2 | \(2.15\times10^{22}\) |
| 100 | \(1.43\times10^{-20}\) | \(-2.26\times10^{-18}\) | 496.523 | 2 | \(2.30\times10^{20}\) |

These digits must not be interpreted as nonzero power laws. All three fits are
rank deficient, and the fitted \(a\) decreases by many orders of magnitude as
the fitting start moves from 400 to 600. The defensible result for this pilot
is

\[
\boxed{a=0},
\]

corresponding to the nonzero memory.

The approach to memory is exponential. A separate semilog diagnostic over
\(t\in[250,450]\) gives

| Observer | \(\gamma\) from \(|p_{\rm loc}|\sim e^{-\gamma t}\) | \(R^2\) |
|---:|---:|---:|
| 0 | 0.10207385 | 0.99999546 |
| 50 | 0.10207889 | 0.99999555 |
| 100 | 0.10209349 | 0.99999582 |

This is close to the cosmological-horizon surface gravity
\(\kappa_c=0.10497495\). Thus the pilot shows memory plus pole-controlled
exponential relaxation, rather than a nonzero asymptotic power law.

## Issues to resolve before the full scan

1. A six-thread production run takes about 90.6 minutes. The 48 cases would
   require approximately 72.5 hours sequentially. We should benchmark
   concurrent three-thread runs and/or optimize translated-Gaussian source
   evaluation before launching the matrix.
2. The fixed nominal fit window is already beyond the resolved exponential
   transient for this pilot. Every run must retain the window-sensitivity and
   identifiability checks; the formal fit output alone is insufficient.
3. A memory solution naturally gives \(a=0\). The analysis should report zero
   when the instantaneous slope has reached the numerical floor, not quote
   observer-dependent values at the \(10^{-20}\) level as physical exponents.

The pilot output and plots are in
`output/sds_areal_scan/q_01_l_0_beta_0/`.

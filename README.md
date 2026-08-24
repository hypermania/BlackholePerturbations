# BlackholePerturbations
My numerical tools for studying blackhole perturbations. They were used to produce the numerical results in papers [arXiv:2503.19967](https://arxiv.org/abs/2503.19967) and [arXiv:2603.20379](https://arxiv.org/abs/2603.20379).

Numerical functionalities implemented for [arXiv:2503.19967](https://arxiv.org/abs/2503.19967):

1. The Regge-Wheeler (RW) equation for a single harmonic mode with a customizable source.
2. The Regge-Wheeler (RW) equations for all real harmonics modes $` \ell \leq \ell_{\mathrm{max}} `$, coupled via a cubic self-interaction.

Numerical functionalities implemented for [arXiv:2603.20379](https://arxiv.org/abs/2603.20379):

1. The Teukosky equation (with change of variable $` \tilde{\psi} = (\Delta^s r) \psi `$) for all spherical harmonic modes $` \ell \leq \ell_{\mathrm{max}} \leq 6 `$. Mode-mode coupling coefficients (for spin-weights $` -2 \leq s \leq 2 `$ and Kerr parameters $` M `$ and $` a `$) are automatically computed and included in the numerical evolution.
2. Two optional sources for the Teukolsky equation, including one that decays like $` r^{-\beta} `$ and one that emulates a Dirac delta function $` \delta(r-r')\delta(t-t') `$.
3. An effective source that corresponds to a $` \lambda \psi^2 `$ nonlinear term for a scalar field in Kerr spacetime.
4. Artifical Kreiss-Oliger dissipation, which is required for stable numerical evolution when $` s \neq 0 `$.

Additional precise CPU functionality:

1. Sourced Schwarzschild-de Sitter master equations for the conformal scalar,
   electromagnetic, and axial gravitational sectors, with the multipole set in
   the equation parameters.
2. Cancellation-free 100-decimal horizon root finding and neighbor-predicted
   tortoise inversion, followed by binary128 time evolution.
3. A standalone `SdSSource` interface for custom callbacks; areal-power,
   horizon-subtracted, local-scalar, and algebraic tortoise profiles (including
   `beta=0`); Gaussian and zero-mean Gaussian-derivative waveforms; and a
   spacetime Gaussian probe for approximating the retarded Green's function.



## Plotting
Some data for the RW equation are already included in the codebase under `/output`. Run the `plot.nb` Mathematica notebook to plot these data. New data can be generated from the code.

## Functionalities
The `src/example.cpp` file has two functions:

1. `run_sourced_eqn()` sets up parameters for the sourced RW equation, and run them in parallel. See section 2.3 of [arXiv:2503.19967](https://arxiv.org/abs/2503.19967) for details.
2. `run_coupled_eqn()` sets up parameters for the coupled RW equations corresponding to a scalar field with cubic self-interaction. See section 3.3 of [arXiv:2503.19967](https://arxiv.org/abs/2503.19967) for details.

The `src/main.cpp` file has three lambda functions in `main()`:

1. `run_teukolsky_sourced` sets up parameters for the sourced Teukolsky equation. See section III.A of [arXiv:2603.20379](https://arxiv.org/abs/2603.20379) for details.
2. `run_teukolsky_dirac_delta` sets up parameters for the Teukolsky equation sourced by a Dirac delta function source $` \delta(r-r')\delta(t-t') `$. Useful for studying the Green's function of the Teukolsky equation.
3. `run_teukolsky_cubic` sets up parameters for the Teukolsky equation corresponding to a scalar field with cubic self-interaction. See section III.B of [arXiv:2603.20379](https://arxiv.org/abs/2603.20379) for details.

## Compilation
Compiler requirement: 

1. A C++ compiler supporting C++20. I used [g++ 12.2.0](https://gcc.gnu.org/).
2. CUDA compiler `nvcc` for compiler GPU kernels. See [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit).

The codebase also makes use of `boost` and `Eigen` library, which are included in `/external`.

Compilation should be as easy as running `make` at the project directory. Note that the compilation could take a while.

## Verification

Run the precise Teukolsky correctness matrix under normal and production flags:

```bash
make check-precise-correctness
```

Run the same focused checks with AddressSanitizer and UndefinedBehaviorSanitizer:

```bash
make check-precise-sanitizers
```

Build the complete CPU-only application with:

```bash
make -j6 disable-cuda=true
```

Run the configured sourced Schwarzschild-de Sitter simulation with:

```bash
./main --sds-precise
```

Edit the configuration block in `run_sds_precise_eqn()` in `src/examples.cpp`
to select $`s`$, $`\ell`$, $`\Lambda`$, the source profile, and the waveform.
The fixed-position time series and snapshots are saved as `double`, matching the
existing precise runners. The source and geometry metadata file records the
horizons, surface gravities, tortoise-coordinate convention, effective Gaussian
cutoff, and conservative finite-domain fitting bound.

Run the Schwarzschild-de Sitter correctness, sanitizer, and production-grid
performance checks with:

```bash
make check-sds-precise-correctness
make check-sds-precise-sanitizers
make check-sds-precise-performance
```

See `doc/SDS_PRECISE_REPORT.md` for the numerical design, validation matrix, and
current scope. See `doc/SDS_OPTIMIZATION_REPORT.md` for the paired attainable
ceiling, Dopri5 stage optimization, timing results, and randomized
full-profile equivalence test.

Run one member of the sourced SdS areal-radius scan with

```bash
./main --sds-areal-scan 0.1 0 0 0
python3 script/analyze_sds_areal_scan.py \
  output/sds_areal_scan/q_01_l_0_beta_0
```

The four arguments are $`q=9\Lambda M^2`$, $`s`$, $`\ell`$, and $`\beta`$.
The runner accepts the scan values $`q\in\{0.1,0.2,0.4,0.8\}`$,
$`s\in\{0,1,2\}`$, $`s\leq\ell\leq3`$, and
$`\beta\in\{0,1,2\}`$. The time series contains $`\psi`$ and $`\Pi`$ at
$`x=50`$, while snapshots contain the complete spatial profile. One complete
set of outputs is stored in its `output/sds_areal_scan/q_XX_l_L_beta_B/`
directory.
The analysis fits

$`p_{\rm loc}=d\ln|\psi|/d\ln t=t\Pi/\psi=a+b/(t-c)`$

and reports the signed tail power $`a`$. Run its focused synthetic check with

```bash
make check-sds-areal-analysis
```

See `doc/SDS_AREAL_SCAN_REPORT.md` for the pilot result and its Jacobian-based
fit assessment. See `doc/CORRECTNESS_REPORT.md` for the tested parameter matrix
and current CUDA and performance-test limitations.

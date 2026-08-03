# Precise Teukolsky Correctness Audit

Audit date: 2026-08-03

## Scope

This audit checks the optimized binary128 Teukolsky CPU path introduced by
commits `e776dc1` and `9494bfc` against the pre-optimization formulas and the
repository's Odeint/Eigen integration. The checks cover:

- the homogeneous right-hand side for four `(s, l, N, ko_epsilon)` parameter
  combinations;
- all fourth-order interior and one-sided boundary stencils;
- generic and separable sources, including source-size validation;
- one-, two-, and six-thread execution;
- 25 complete Dopri5 steps against an independent original-formula system;
- `scale_sum1` through `scale_sum7` for Eigen arrays and matrices with `double`
  and binary128 scalars;
- sizes immediately below and at the 4,096-element OpenMP threshold, uneven
  parallel partitions, and aliased output; and
- normal, production (`-O3 -ffast-math -DNDEBUG`), ASan, and UBSan builds.

The reference right-hand side in
`test/test_teukolsky_precise_correctness.cpp` spells out the original finite
difference expressions independently. It does not call the optimized interior
kernel. The scale-sum test computes each expected element independently of the
Odeint operation implementation.

## Results

Run the normal and production matrix with:

```bash
make check-precise-correctness
```

Both configurations pass. The largest scaled RHS/trajectory difference is
`2.17e-34` with normal flags and `2.89e-34` with production flags. Thread-count
comparisons and generic-versus-separable source comparisons are exact.

Run memory and undefined-behavior checks with:

```bash
make check-precise-sanitizers
```

ASan and UBSan pass both test programs without findings. GCC ThreadSanitizer is
not a valid OpenMP check in the installed GCC/libgomp configuration: it does not
model the implicit `parallel for` barrier and reports worker writes as racing
with main-thread reads after the barrier.

The full CPU-only application builds and links with production flags:

```bash
make -j6 disable-cuda=true
```

This required guarding CUDA-only includes and observers in `src/main.cpp` and
`src/observer.hpp`; previously `DISABLE_CUDA` still unconditionally required
Thrust headers and types.

## Residual Infrastructure Findings

The CUDA build remains incompatible with CUDA 13.2 because
`src/cuda_wrapper.cuh` and `src/cuda_wrapper.cu` explicitly instantiate the
removed private Thrust type `thrust::detail::device_generate_functor`. This is a
pre-existing CUDA compatibility issue and is unrelated to the CPU optimization.
The host also has no NVIDIA runtime device, so CUDA execution cannot be checked.

The six-thread performance acceptance gate was numerically correct but did not
reproduce its earlier 90% homogeneous-ceiling result on this shared VPS: two
runs measured 78.9% and 66.1%. A four-thread diagnostic run passed at 119.4%.
This does not affect correctness, but it confirms that the strict six-core
performance gate remains sensitive to host scheduling and should not be treated
as a functional test result.

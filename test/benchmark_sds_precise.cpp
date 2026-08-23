#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

#include "odeint_eigen/eigen_operations.hpp"
#include "sds_precise.hpp"

namespace {

using Equation = SdSMasterPDEPrecise;
using Param = SdSMasterPDEPreciseParam;
using Scalar = Equation::Scalar;
using State = Equation::State;
using Vector = Equation::Vector;

double process_cpu_seconds() {
  timespec time;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time);
  return static_cast<double>(time.tv_sec)
         + static_cast<double>(time.tv_nsec) * 1e-9;
}

double median(std::vector<double> samples) {
  std::sort(samples.begin(), samples.end());
  return samples[samples.size() / 2];
}

Param make_param(const long long int n) {
  Param param;
  param.s = 0;
  param.l = 0;
  param.M = Scalar("0.5");
  param.Lambda = Scalar("0.0444444444444444444444444444444444");
  param.r_min = -500;
  param.r_max = 1000;
  param.N = n;
  param.t_start = 0;
  param.t_end = 1;
  param.t_interval = Scalar("0.5");
  param.delta_t = Scalar("0.01");
  return param;
}

double run_operator(Equation &equation, const State &state, State &derivative,
                    const int iterations, const Scalar &time) {
  equation(state, derivative, time);
  const auto start = std::chrono::steady_clock::now();
  for(int iteration = 0; iteration < iterations; ++iteration) {
    equation(state, derivative, time);
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();
}

void ceiling_rhs(const Equation &equation, const State &state,
                 State &derivative, const Scalar &time,
                 Vector &source_workspace) {
  const long long int grid_size = equation.grid_size;
  const Scalar *__restrict__ psi = state.data();
  const Scalar *__restrict__ pi = state.data() + grid_size;
  Scalar *__restrict__ dpsi = derivative.data();
  Scalar *__restrict__ dpi = derivative.data() + grid_size;
  const Scalar *__restrict__ potential = equation.V.data();
  const Scalar d2_factor = equation.inv_h_sqr / Scalar(12);
  const Scalar near_factor = Scalar(16) * d2_factor;
  const Scalar far_factor = -d2_factor;
  const Scalar center_coefficient = -Scalar(30) * d2_factor;
  if(equation.Q) equation.Q(time, source_workspace);
  const Scalar *__restrict__ source = source_workspace.data();

#pragma omp parallel for schedule(static)
  for(long long int i = 2; i <= grid_size - 3; ++i) {
    const Scalar near_sum = psi[i - 1] + psi[i + 1];
    const Scalar far_sum = psi[i - 2] + psi[i + 2];
    dpi[i] = near_factor * near_sum + far_factor * far_sum
             + (center_coefficient - potential[i]) * psi[i] + source[i];
    dpsi[i] = pi[i];
  }

}

double run_ceiling(const Equation &equation, const State &state,
                   State &derivative, const int iterations,
                   const Scalar &time, Vector &source_workspace) {
  ceiling_rhs(equation, state, derivative, time, source_workspace);
  const auto start = std::chrono::steady_clock::now();
  for(int iteration = 0; iteration < iterations; ++iteration) {
    ceiling_rhs(equation, state, derivative, time, source_workspace);
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();
}

std::pair<double, double> run_paired_rhs(
    Equation &equation, const State &state, State &derivative,
    const int iterations, const Scalar &time, const bool operator_first) {
  constexpr int target_block_size = 5;
  const int blocks = std::max(1, (iterations + target_block_size - 1)
                                  / target_block_size);
  const int base_iterations = iterations / blocks;
  const int extra_iterations = iterations % blocks;
  Vector ceiling_source = Vector::Zero(equation.grid_size);
  double operator_seconds = 0;
  double ceiling_seconds = 0;
  for(int block = 0; block < blocks; ++block) {
    const int block_iterations = base_iterations
                                 + (block < extra_iterations ? 1 : 0);
    const bool production_then_ceiling = operator_first == (block % 2 == 0);
    if(production_then_ceiling) {
      operator_seconds += run_operator(
          equation, state, derivative, block_iterations, time);
      ceiling_seconds += run_ceiling(
          equation, state, derivative, block_iterations, time, ceiling_source);
    } else {
      ceiling_seconds += run_ceiling(
          equation, state, derivative, block_iterations, time, ceiling_source);
      operator_seconds += run_operator(
          equation, state, derivative, block_iterations, time);
    }
  }
  return {operator_seconds, ceiling_seconds};
}

double run_source(const SdSSource &source, Vector &workspace,
                  const int iterations, const Scalar &time) {
  source(time, workspace);
  const auto start = std::chrono::steady_clock::now();
  for(int iteration = 0; iteration < iterations; ++iteration) {
    source(time + Scalar(iteration) * Scalar("1e-4"), workspace);
  }
  return std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
}

struct DopriResult {
  double wall_seconds;
  double cpu_seconds;
  State state;
};

template<bool OptimizeUnitCoefficient = true>
DopriResult run_dopri5_steps(Equation &equation,
                             const State &initial_state,
                             const int threads,
                             const int steps) {
  using Operations = boost::numeric::odeint::eigen_operations<
      State, OptimizeUnitCoefficient>;
  using Stepper = boost::numeric::odeint::runge_kutta_dopri5<
      State, Scalar, State, Scalar,
      boost::numeric::odeint::vector_space_algebra, Operations>;
  omp_set_num_threads(threads);
  State state = initial_state;
  Stepper stepper;
  Scalar time = 50;
  const Scalar dt("0.0001");
  stepper.do_step(std::ref(equation), state, time, dt);
  const double cpu_start = process_cpu_seconds();
  const auto start = std::chrono::steady_clock::now();
  for(int step = 0; step < steps; ++step) {
    time += dt;
    stepper.do_step(std::ref(equation), state, time, dt);
  }
  const double wall_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  const double cpu_seconds = process_cpu_seconds() - cpu_start;
  return {wall_seconds, cpu_seconds, std::move(state)};
}

Scalar max_difference(const State &left, const State &right) {
  Scalar result = 0;
  for(Eigen::Index i = 0; i < left.size(); ++i) {
    result = std::max(result, abs(left[i] - right[i]));
  }
  return result;
}

struct PairedDopriResult {
  double optimized_wall_seconds;
  double legacy_wall_seconds;
  double wall_speedup;
  double optimized_cpu_seconds;
  double legacy_cpu_seconds;
  double cpu_speedup;
};

PairedDopriResult median_paired_dopri5_step_seconds(
    Equation &equation, const State &initial_state, const int threads,
    const int steps, const int samples) {
  std::vector<double> optimized_wall_timings;
  std::vector<double> legacy_wall_timings;
  std::vector<double> wall_speedups;
  std::vector<double> optimized_cpu_timings;
  std::vector<double> legacy_cpu_timings;
  std::vector<double> cpu_speedups;
  optimized_wall_timings.reserve(samples);
  legacy_wall_timings.reserve(samples);
  wall_speedups.reserve(samples);
  optimized_cpu_timings.reserve(samples);
  legacy_cpu_timings.reserve(samples);
  cpu_speedups.reserve(samples);
  for(int sample = 0; sample < samples; ++sample) {
    DopriResult optimized;
    DopriResult legacy;
    if(sample % 2 == 0) {
      optimized = run_dopri5_steps<true>(
          equation, initial_state, threads, steps);
      legacy = run_dopri5_steps<false>(
          equation, initial_state, threads, steps);
    } else {
      legacy = run_dopri5_steps<false>(
          equation, initial_state, threads, steps);
      optimized = run_dopri5_steps<true>(
          equation, initial_state, threads, steps);
    }
    optimized_wall_timings.push_back(optimized.wall_seconds / steps);
    legacy_wall_timings.push_back(legacy.wall_seconds / steps);
    wall_speedups.push_back(legacy.wall_seconds / optimized.wall_seconds);
    optimized_cpu_timings.push_back(optimized.cpu_seconds / steps);
    legacy_cpu_timings.push_back(legacy.cpu_seconds / steps);
    cpu_speedups.push_back(legacy.cpu_seconds / optimized.cpu_seconds);
  }
  return {median(std::move(optimized_wall_timings)),
          median(std::move(legacy_wall_timings)),
          median(std::move(wall_speedups)),
          median(std::move(optimized_cpu_timings)),
          median(std::move(legacy_cpu_timings)),
          median(std::move(cpu_speedups))};
}

}  // namespace

int main(int argc, char **argv) {
  const long long int n = argc > 1 ? std::stoll(argv[1]) : 50000;
  const int iterations = argc > 2 ? std::stoi(argv[2]) : 60;
  const int max_threads = argc > 3 ? std::stoi(argv[3]) : omp_get_max_threads();
  omp_set_dynamic(0);

  const auto initialization_start = std::chrono::steady_clock::now();
  Equation equation(make_param(n));
  const double initialization_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - initialization_start).count();

  State state(2 * equation.grid_size);
  State derivative(2 * equation.grid_size);
  Vector ceiling_source = Vector::Zero(equation.grid_size);
  for(long long int i = 0; i < equation.grid_size; ++i) {
    state[i] = Scalar("0.25") + Scalar(i % 31) * Scalar("1e-4");
    state[equation.grid_size + i] = Scalar("-0.125")
                                    + Scalar(i % 23) * Scalar("2e-4");
  }
  const Scalar time(50);

  omp_set_num_threads(1);
  const double one_thread = run_operator(
      equation, state, derivative, iterations, time);
  const State one_thread_reference = derivative;

  omp_set_num_threads(max_threads);
  constexpr int samples = 7;
  std::vector<double> operator_samples;
  std::vector<double> ceiling_samples;
  std::vector<double> efficiencies;
  for(int sample = 0; sample < samples; ++sample) {
    const auto [operator_seconds, ceiling_seconds] = run_paired_rhs(
        equation, state, derivative, iterations, time, sample % 2 == 0);
    operator_samples.push_back(operator_seconds);
    ceiling_samples.push_back(ceiling_seconds);
    efficiencies.push_back(ceiling_seconds / operator_seconds);
  }
  const double homogeneous_seconds = median(operator_samples);
  const double homogeneous_ceiling_seconds = median(ceiling_samples);
  const double homogeneous_efficiency = median(efficiencies);
  equation(state, derivative, time);
  const State homogeneous_reference = derivative;
  const Scalar thread_error = max_difference(derivative, one_thread_reference);
  ceiling_rhs(equation, state, derivative, time, ceiling_source);
  const Scalar ceiling_error = max_difference(derivative, homogeneous_reference);

  const auto homogeneous_step_result =
      run_dopri5_steps(equation, state, max_threads, 20);

  SdSTranslatedSourceParam source;
  source.profile = SdSSourceProfile::ArealPower;
  source.waveform = SdSWaveform::Gaussian;
  source.beta = 0;
  source.u_center = -10;
  source.sigma = Scalar("0.5");
  source.cutoff_sigma = 12;
  SdSSource translated_source(source);
  translated_source.initialize(
      equation.param.r_min, equation.grid_space(), equation.grid_size,
      equation.r_cosmological, equation.r, equation.rho_cosmological,
      equation.f);
  equation.Q = std::move(translated_source);

  Vector source_workspace = Vector::Zero(equation.grid_size);
  omp_set_num_threads(1);
  const double one_thread_source_seconds = run_source(
      equation.Q, source_workspace, iterations, time);
  omp_set_num_threads(max_threads);
  const double source_seconds = run_source(
      equation.Q, source_workspace, iterations, time);

  operator_samples.clear();
  ceiling_samples.clear();
  efficiencies.clear();
  for(int sample = 0; sample < samples; ++sample) {
    const auto [operator_seconds, ceiling_seconds] = run_paired_rhs(
        equation, state, derivative, iterations, time, sample % 2 == 0);
    operator_samples.push_back(operator_seconds);
    ceiling_samples.push_back(ceiling_seconds);
    efficiencies.push_back(ceiling_seconds / operator_seconds);
  }
  const double sourced_seconds = median(operator_samples);
  const double sourced_ceiling_seconds = median(ceiling_samples);
  const double sourced_efficiency = median(efficiencies);
  equation(state, derivative, time);
  const State sourced_reference = derivative;
  ceiling_rhs(equation, state, derivative, time, ceiling_source);
  const Scalar sourced_ceiling_error = max_difference(derivative, sourced_reference);

  const auto one_thread_step_result =
      run_dopri5_steps(equation, state, 1, 5);
  const auto many_thread_step_result =
      run_dopri5_steps(equation, state, max_threads, 5);
  const PairedDopriResult production_steps =
      median_paired_dopri5_step_seconds(
          equation, state, max_threads, 10, 7);
  const auto legacy_validation_result =
      run_dopri5_steps<false>(equation, state, max_threads, 5);
  const Scalar step_error = max_difference(
      one_thread_step_result.state, many_thread_step_result.state);
  const Scalar legacy_step_error = max_difference(
      many_thread_step_result.state, legacy_validation_result.state);

  const double points = static_cast<double>(n - 3) * iterations;
  std::cout << std::setprecision(8)
            << "grid_points=" << equation.grid_size << '\n'
            << "iterations=" << iterations << '\n'
            << "threads=" << max_threads << '\n'
            << "geometry_initialization_seconds=" << initialization_seconds << '\n'
            << "one_thread_points_per_second=" << points / one_thread << '\n'
            << "homogeneous_points_per_second=" << points / homogeneous_seconds << '\n'
            << "homogeneous_ceiling_points_per_second="
            << points / homogeneous_ceiling_seconds << '\n'
            << "homogeneous_ceiling_efficiency=" << homogeneous_efficiency << '\n'
            << "sourced_points_per_second=" << points / sourced_seconds << '\n'
            << "sourced_ceiling_points_per_second="
            << points / sourced_ceiling_seconds << '\n'
            << "sourced_ceiling_efficiency=" << sourced_efficiency << '\n'
            << "source_evaluations_per_second="
            << iterations / source_seconds << '\n'
            << "one_thread_source_evaluations_per_second="
            << iterations / one_thread_source_seconds << '\n'
            << "homogeneous_dopri5_step_seconds="
            << homogeneous_step_result.wall_seconds / 20 << '\n'
            << "dopri5_speedup="
            << one_thread_step_result.wall_seconds
                   / many_thread_step_result.wall_seconds << '\n'
            << "one_thread_dopri5_step_seconds="
            << one_thread_step_result.wall_seconds / 5 << '\n'
            << "many_thread_dopri5_step_seconds="
            << many_thread_step_result.wall_seconds / 5 << '\n'
            << "production_median_dopri5_step_seconds="
            << production_steps.optimized_wall_seconds << '\n'
            << "legacy_median_dopri5_step_seconds="
            << production_steps.legacy_wall_seconds << '\n'
            << "dopri5_optimization_speedup="
            << production_steps.wall_speedup << '\n'
            << "production_median_dopri5_step_cpu_seconds="
            << production_steps.optimized_cpu_seconds << '\n'
            << "legacy_median_dopri5_step_cpu_seconds="
            << production_steps.legacy_cpu_seconds << '\n'
            << "dopri5_cpu_optimization_speedup="
            << production_steps.cpu_speedup << '\n'
            << "estimated_100000_step_minutes="
            << production_steps.optimized_wall_seconds * 100000 / 60 << '\n'
            << "max_thread_abs_error=" << thread_error << '\n'
            << "max_homogeneous_ceiling_abs_error=" << ceiling_error << '\n'
            << "max_sourced_ceiling_abs_error=" << sourced_ceiling_error << '\n'
            << "max_dopri5_abs_error=" << step_error << '\n'
            << "max_legacy_dopri5_abs_error=" << legacy_step_error << '\n';

  return thread_error == 0 && step_error == 0 && legacy_step_error == 0
                 && ceiling_error < Scalar("1e-28")
                 && sourced_ceiling_error < Scalar("1e-28")
                 && homogeneous_efficiency >= 0.90
                 && sourced_efficiency >= 0.90
             ? EXIT_SUCCESS : EXIT_FAILURE;
}

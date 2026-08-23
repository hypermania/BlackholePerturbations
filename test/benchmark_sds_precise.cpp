#include <algorithm>
#include <chrono>
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

double median(std::vector<double> samples) {
  std::sort(samples.begin(), samples.end());
  return samples[samples.size() / 2];
}

Param make_param(const long long int n) {
  Param param;
  param.s = 2;
  param.l = 2;
  param.M = Scalar("0.5");
  param.Lambda = Scalar("1e-4");
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
  const Scalar d1_factor = equation.inv_h / Scalar(12);
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

  dpi[0] = (-Scalar(25) * pi[0] + Scalar(48) * pi[1]
            - Scalar(36) * pi[2] + Scalar(16) * pi[3]
            - Scalar(3) * pi[4]) * d1_factor
           - potential[0] * psi[0] + source[0];
  dpsi[0] = pi[0];
  dpi[1] = (Scalar(11) * psi[0] - Scalar(20) * psi[1]
            + Scalar(6) * psi[2] + Scalar(4) * psi[3] - psi[4])
           * d2_factor - potential[1] * psi[1] + source[1];
  dpsi[1] = pi[1];

  const long long int n2 = grid_size - 2;
  const long long int n1 = grid_size - 1;
  dpi[n2] = (-psi[grid_size - 5] + Scalar(4) * psi[grid_size - 4]
             + Scalar(6) * psi[grid_size - 3] - Scalar(20) * psi[n2]
             + Scalar(11) * psi[n1]) * d2_factor
            - potential[n2] * psi[n2] + source[n2];
  dpsi[n2] = pi[n2];
  dpi[n1] = (-Scalar(3) * pi[grid_size - 5]
             + Scalar(16) * pi[grid_size - 4]
             - Scalar(36) * pi[grid_size - 3] + Scalar(48) * pi[n2]
             - Scalar(25) * pi[n1]) * d1_factor
            - potential[n1] * psi[n1] + source[n1];
  dpsi[n1] = pi[n1];
}

double run_ceiling(const Equation &equation, const State &state,
                   State &derivative, const int iterations,
                   const Scalar &time) {
  Vector source_workspace = Vector::Zero(equation.grid_size);
  ceiling_rhs(equation, state, derivative, time, source_workspace);
  const auto start = std::chrono::steady_clock::now();
  for(int iteration = 0; iteration < iterations; ++iteration) {
    ceiling_rhs(equation, state, derivative, time, source_workspace);
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();
}

std::pair<double, State> run_dopri5_steps(Equation &equation,
                                          const State &initial_state,
                                          const int threads,
                                          const int steps) {
  using Stepper = boost::numeric::odeint::runge_kutta_dopri5<
      State, Scalar, State, Scalar>;
  omp_set_num_threads(threads);
  State state = initial_state;
  Stepper stepper;
  Scalar time = 50;
  const Scalar dt("0.0001");
  stepper.do_step(std::ref(equation), state, time, dt);
  const auto start = std::chrono::steady_clock::now();
  for(int step = 0; step < steps; ++step) {
    time += dt;
    stepper.do_step(std::ref(equation), state, time, dt);
  }
  return {std::chrono::duration<double>(
              std::chrono::steady_clock::now() - start).count(),
          std::move(state)};
}

Scalar max_difference(const State &left, const State &right) {
  Scalar result = 0;
  for(Eigen::Index i = 0; i < left.size(); ++i) {
    result = std::max(result, abs(left[i] - right[i]));
  }
  return result;
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
    double operator_seconds;
    double ceiling_seconds;
    if(sample % 2 == 0) {
      operator_seconds = run_operator(equation, state, derivative, iterations, time);
      ceiling_seconds = run_ceiling(
          equation, state, derivative, iterations, time);
    } else {
      ceiling_seconds = run_ceiling(
          equation, state, derivative, iterations, time);
      operator_seconds = run_operator(equation, state, derivative, iterations, time);
    }
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

  SdSTranslatedSourceParam source;
  source.profile = SdSSourceProfile::TortoisePower;
  source.waveform = SdSWaveform::Gaussian;
  source.beta = 2;
  source.u_center = -10;
  source.sigma = Scalar("0.5");
  source.cutoff_sigma = 12;
  source.L = 1;
  source.x0 = 100;
  source.X0 = 0;
  source.X1 = 20;
  SdSSource translated_source(source);
  translated_source.initialize(
      equation.param.r_min, equation.grid_space(), equation.grid_size,
      equation.r_cosmological, equation.r, equation.rho_cosmological,
      equation.f);
  equation.Q = std::move(translated_source);

  operator_samples.clear();
  ceiling_samples.clear();
  efficiencies.clear();
  for(int sample = 0; sample < samples; ++sample) {
    double operator_seconds;
    double ceiling_seconds;
    if(sample % 2 == 0) {
      operator_seconds = run_operator(equation, state, derivative, iterations, time);
      ceiling_seconds = run_ceiling(
          equation, state, derivative, iterations, time);
    } else {
      ceiling_seconds = run_ceiling(
          equation, state, derivative, iterations, time);
      operator_seconds = run_operator(equation, state, derivative, iterations, time);
    }
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

  const auto [one_thread_step_seconds, one_thread_step_state] =
      run_dopri5_steps(equation, state, 1, 5);
  const auto [many_thread_step_seconds, many_thread_step_state] =
      run_dopri5_steps(equation, state, max_threads, 5);
  const Scalar step_error = max_difference(
      one_thread_step_state, many_thread_step_state);

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
            << "dopri5_speedup="
            << one_thread_step_seconds / many_thread_step_seconds << '\n'
            << "one_thread_dopri5_step_seconds="
            << one_thread_step_seconds / 5 << '\n'
            << "many_thread_dopri5_step_seconds="
            << many_thread_step_seconds / 5 << '\n'
            << "max_thread_abs_error=" << thread_error << '\n'
            << "max_homogeneous_ceiling_abs_error=" << ceiling_error << '\n'
            << "max_sourced_ceiling_abs_error=" << sourced_ceiling_error << '\n'
            << "max_dopri5_abs_error=" << step_error << '\n';

  return thread_error == 0 && step_error == 0
                 && ceiling_error < Scalar("1e-28")
                 && sourced_ceiling_error < Scalar("1e-28")
                 && homogeneous_efficiency >= 0.90
                 && sourced_efficiency >= 0.90
             ? EXIT_SUCCESS : EXIT_FAILURE;
}

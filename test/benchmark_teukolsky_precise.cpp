#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <omp.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include "odeint_eigen/eigen_operations.hpp"
#include "teukolsky_precise.hpp"

namespace {

using Equation = TeukolskyPDEPrecise;
using Param = TeukolskyPDEPreciseParam;
using Scalar = Equation::Scalar;
using State = Equation::State;

double median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  return values[values.size() / 2];
}

Param make_param(long long n) {
  Param param;
  param.s = -1;
  param.l = 1;
  param.M = Scalar("0.5");
  param.r_min = Scalar(-500);
  param.r_max = Scalar(1000);
  param.N = n;
  param.ko_epsilon = Scalar("0.7");
  param.t_start = Scalar(0);
  param.t_end = Scalar(1);
  param.t_interval = Scalar("0.5");
  param.delta_t = Scalar("0.01");
  return param;
}

double run(Equation &equation, const State &state, State &derivative,
           int iterations) {
  for (int i = 0; i < 3; ++i) {
    equation(state, derivative, Scalar(i) / Scalar(10));
  }
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < iterations; ++i) {
    equation(state, derivative, Scalar(i) / Scalar(10));
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();
}

double run_interior_ceiling(const Equation &equation, const State &state,
                            State &derivative, int iterations,
                            const Scalar *source_spatial = nullptr,
                            Scalar source_scale = 0) {
  const long long grid_size = equation.grid_size;
  const Scalar d2_factor = equation.inv_h_sqr / Scalar(12);
  const Scalar near_factor = Scalar(16) * d2_factor;
  const Scalar far_factor = -d2_factor;
  const Scalar ko_factor = equation.param.ko_epsilon / Scalar(16);
  const Scalar *psi = state.data();
  const Scalar *pi = state.data() + grid_size;
  Scalar *dpsi = derivative.data();
  Scalar *dpi = derivative.data() + grid_size;
  for (int repeat = 0; repeat < 3; ++repeat) {
#pragma omp parallel for schedule(static)
    for (long long i = 2; i <= grid_size - 3; ++i) {
      const Scalar near_sum = psi[i - 1] + psi[i + 1];
      const Scalar far_sum = psi[i - 2] + psi[i + 2];
      const Scalar gradient = Scalar(8) * (psi[i + 1] - psi[i - 1])
                              + psi[i - 2] - psi[i + 2];
      dpi[i] = near_factor * near_sum + far_factor * far_sum
               + equation.center_factor[i] * psi[i]
               - equation.C_d1_factor[i] * gradient - equation.A[i] * pi[i];
      if (source_spatial) dpi[i] += source_scale * source_spatial[i];
      dpsi[i] = pi[i] - ko_factor * (far_sum - Scalar(4) * near_sum
                                     + Scalar(6) * psi[i]);
    }
  }
  const auto start = std::chrono::steady_clock::now();
  for (int repeat = 0; repeat < iterations; ++repeat) {
#pragma omp parallel for schedule(static)
    for (long long i = 2; i <= grid_size - 3; ++i) {
      const Scalar near_sum = psi[i - 1] + psi[i + 1];
      const Scalar far_sum = psi[i - 2] + psi[i + 2];
      const Scalar gradient = Scalar(8) * (psi[i + 1] - psi[i - 1])
                              + psi[i - 2] - psi[i + 2];
      dpi[i] = near_factor * near_sum + far_factor * far_sum
               + equation.center_factor[i] * psi[i]
               - equation.C_d1_factor[i] * gradient - equation.A[i] * pi[i];
      if (source_spatial) dpi[i] += source_scale * source_spatial[i];
      dpsi[i] = pi[i] - ko_factor * (far_sum - Scalar(4) * near_sum
                                     + Scalar(6) * psi[i]);
    }
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();
}

std::pair<double, State> run_dopri5_steps(Equation &equation,
                                          const State &initial_state,
                                          int threads, int steps) {
  using Stepper = boost::numeric::odeint::runge_kutta_dopri5<
      State, Scalar, State, Scalar>;
  omp_set_num_threads(threads);
  State state = initial_state;
  Stepper stepper;
  Scalar time = 0;
  const Scalar dt("0.0001");
  stepper.do_step(std::ref(equation), state, time, dt);
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < steps; ++i) {
    time += dt;
    stepper.do_step(std::ref(equation), state, time, dt);
  }
  return {std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
              .count(),
          std::move(state)};
}

template <bool Multiply>
double measure_binary128_rate(int threads, long long operations_per_thread) {
  constexpr int lanes = 8;
  const long long rounds = operations_per_thread / lanes;
  State sinks(threads);
  const auto start = std::chrono::steady_clock::now();
#pragma omp parallel num_threads(threads)
  {
    Scalar x0("1.00000001"), x1("1.00000002"), x2("1.00000003"), x3("1.00000004");
    Scalar x4("1.00000005"), x5("1.00000006"), x6("1.00000007"), x7("1.00000008");
    const Scalar operand = Multiply ? Scalar("1.000000000000000000000000000001")
                                    : Scalar("0.000000000000000000000000000001");
    for (long long i = 0; i < rounds; ++i) {
      if constexpr (Multiply) {
        x0 *= operand; x1 *= operand; x2 *= operand; x3 *= operand;
        x4 *= operand; x5 *= operand; x6 *= operand; x7 *= operand;
      } else {
        x0 += operand; x1 += operand; x2 += operand; x3 += operand;
        x4 += operand; x5 += operand; x6 += operand; x7 += operand;
      }
    }
    sinks[omp_get_thread_num()] = x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;
  }
  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  volatile double keep = 0;
  for (int i = 0; i < threads; ++i) keep += sinks[i].convert_to<double>();
  (void)keep;
  return static_cast<double>(rounds * lanes) * threads / seconds;
}

}  // namespace

int main(int argc, char **argv) {
  const long long n = argc > 1 ? std::stoll(argv[1]) : 50000;
  const int iterations = argc > 2 ? std::stoi(argv[2]) : 60;
  const int max_threads = argc > 3 ? std::stoi(argv[3]) : omp_get_max_threads();

  omp_set_dynamic(0);
  Equation equation(make_param(n));

  State state(2 * (n + 1));
  State derivative(2 * (n + 1));
  for (long long i = 0; i < n + 1; ++i) {
    state[i] = Scalar("0.25") + Scalar(i % 31) * Scalar("1e-4");
    state[n + 1 + i] = Scalar("-0.125") + Scalar(i % 23) * Scalar("2e-4");
  }

  omp_set_num_threads(1);
  const double one_thread = run(equation, state, derivative, iterations);
  const State threaded_reference = derivative;

  State formula_reference = derivative;
  const long long grid_size = n + 1;
  const Scalar one_twelfth = Scalar(1) / Scalar(12);
  const Scalar ko_factor = equation.param.ko_epsilon / Scalar(16);
  for (long long i = 2; i <= grid_size - 3; ++i) {
    const Scalar d2 = (-state(i - 2) + Scalar(16) * state(i - 1)
                       - Scalar(30) * state(i) + Scalar(16) * state(i + 1)
                       - state(i + 2)) * one_twelfth * equation.inv_h_sqr;
    const Scalar d1 = (state(i - 2) - Scalar(8) * state(i - 1)
                       + Scalar(8) * state(i + 1) - state(i + 2))
                      * one_twelfth * equation.inv_h;
    formula_reference(grid_size + i) =
        d2 - equation.C[i] * d1 - equation.A[i] * state(grid_size + i)
        - equation.V[i] * state(i);
    formula_reference(i) = state(grid_size + i)
                           - ko_factor * (state(i - 2) - Scalar(4) * state(i - 1)
                                          + Scalar(6) * state(i)
                                          - Scalar(4) * state(i + 1) + state(i + 2));
  }

  omp_set_num_threads(max_threads);
  const int best_threads = max_threads;
  std::vector<double> operator_samples;
  std::vector<double> ceiling_samples;
  std::vector<double> efficiency_samples;
  constexpr int performance_samples = 5;
  for (int sample = 0; sample < performance_samples; ++sample) {
    if (sample % 2 == 0) {
      operator_samples.push_back(run(equation, state, derivative, iterations));
      ceiling_samples.push_back(
          run_interior_ceiling(equation, state, derivative, iterations));
    } else {
      ceiling_samples.push_back(
          run_interior_ceiling(equation, state, derivative, iterations));
      operator_samples.push_back(run(equation, state, derivative, iterations));
    }
    efficiency_samples.push_back(ceiling_samples.back() / operator_samples.back());
  }
  const double many_threads = median(operator_samples);
  const double interior_ceiling_seconds = median(ceiling_samples);

  equation(state, derivative, Scalar(0));

  Scalar max_thread_error = 0;
  Scalar max_formula_error = 0;
  for (Eigen::Index i = 0; i < derivative.size(); ++i) {
    max_thread_error = std::max(max_thread_error,
                                abs(derivative[i] - threaded_reference[i]));
    max_formula_error = std::max(max_formula_error,
                                 abs(derivative[i] - formula_reference[i]));
  }
  const State optimized_derivative = derivative;
  run_interior_ceiling(equation, state, derivative, 1);
  Scalar max_ceiling_error = 0;
  for (Eigen::Index i = 0; i < derivative.size(); ++i) {
    max_ceiling_error = std::max(max_ceiling_error,
                                 abs(derivative[i] - optimized_derivative[i]));
  }

  const double speedup = one_thread / many_threads;
  const double scaling_efficiency = speedup / best_threads;
  const double points = static_cast<double>(n - 3) * iterations;
  const double add_rate = measure_binary128_rate<false>(best_threads, 8000000);
  const double multiply_rate = measure_binary128_rate<true>(best_threads, 8000000);
  constexpr double additions_per_point = 12.0;
  constexpr double multiplications_per_point = 9.0;
  const double roofline_points_per_second =
      1.0 / (additions_per_point / add_rate
             + multiplications_per_point / multiply_rate);
  const double measured_points_per_second = points / many_threads;
  const double roofline_efficiency = measured_points_per_second
                                     / roofline_points_per_second;
  const double interior_ceiling_points_per_second =
      points / interior_ceiling_seconds;
  const double interior_ceiling_efficiency = median(efficiency_samples);

  State source_spatial(grid_size);
  for (long long i = 0; i < grid_size; ++i) {
    source_spatial[i] = Scalar(i % 19) * Scalar("1e-6");
  }
  const Scalar source_scale("0.75");
  equation.set_separable_source(source_spatial,
      [source_scale](const Scalar &) { return source_scale; });
  std::vector<double> sourced_operator_samples;
  std::vector<double> sourced_ceiling_samples;
  std::vector<double> sourced_efficiency_samples;
  for (int sample = 0; sample < performance_samples; ++sample) {
    if (sample % 2 == 0) {
      sourced_operator_samples.push_back(run(equation, state, derivative, iterations));
      sourced_ceiling_samples.push_back(run_interior_ceiling(
          equation, state, derivative, iterations, source_spatial.data(), source_scale));
    } else {
      sourced_ceiling_samples.push_back(run_interior_ceiling(
          equation, state, derivative, iterations, source_spatial.data(), source_scale));
      sourced_operator_samples.push_back(run(equation, state, derivative, iterations));
    }
    sourced_efficiency_samples.push_back(
        sourced_ceiling_samples.back() / sourced_operator_samples.back());
  }
  const double sourced_seconds = median(sourced_operator_samples);
  const double sourced_ceiling_seconds = median(sourced_ceiling_samples);
  const double sourced_points_per_second = points / sourced_seconds;
  const double sourced_ceiling_points_per_second = points / sourced_ceiling_seconds;
  const double sourced_ceiling_efficiency = median(sourced_efficiency_samples);

  equation(state, derivative, Scalar(0));
  Scalar max_source_error = 0;
  for (long long i = 0; i < grid_size; ++i) {
    max_source_error = std::max(
        max_source_error,
        abs(derivative[grid_size + i]
            - (optimized_derivative[grid_size + i] + source_scale * source_spatial[i])));
    max_source_error = std::max(
        max_source_error, abs(derivative[i] - optimized_derivative[i]));
  }
  const auto [one_thread_step_seconds, one_thread_step_state] =
      run_dopri5_steps(equation, state, 1, 5);
  const auto [many_thread_step_seconds, many_thread_step_state] =
      run_dopri5_steps(equation, state, best_threads, 5);
  Scalar max_step_error = 0;
  for (Eigen::Index i = 0; i < one_thread_step_state.size(); ++i) {
    max_step_error = std::max(max_step_error,
                              abs(one_thread_step_state[i] - many_thread_step_state[i]));
  }
  std::cout << std::setprecision(8)
            << "grid_points=" << n + 1 << '\n'
            << "iterations=" << iterations << '\n'
            << "available_threads=" << max_threads << '\n'
            << "best_threads=" << best_threads << '\n'
            << "one_thread_seconds=" << one_thread << '\n'
            << "many_thread_seconds=" << many_threads << '\n'
            << "one_thread_points_per_second=" << points / one_thread << '\n'
            << "many_thread_points_per_second=" << points / many_threads << '\n'
            << "speedup=" << speedup << '\n'
            << "scaling_efficiency=" << scaling_efficiency << '\n'
            << "binary128_additions_per_second=" << add_rate << '\n'
            << "binary128_multiplications_per_second=" << multiply_rate << '\n'
            << "roofline_points_per_second=" << roofline_points_per_second << '\n'
            << "roofline_efficiency=" << roofline_efficiency << '\n'
            << "interior_ceiling_points_per_second="
            << interior_ceiling_points_per_second << '\n'
            << "interior_ceiling_efficiency=" << interior_ceiling_efficiency << '\n'
            << "sourced_points_per_second=" << sourced_points_per_second << '\n'
            << "sourced_ceiling_points_per_second="
            << sourced_ceiling_points_per_second << '\n'
            << "sourced_ceiling_efficiency=" << sourced_ceiling_efficiency << '\n'
            << "one_thread_dopri5_step_seconds=" << one_thread_step_seconds / 5 << '\n'
            << "many_thread_dopri5_step_seconds=" << many_thread_step_seconds / 5 << '\n'
            << "dopri5_speedup=" << one_thread_step_seconds / many_thread_step_seconds << '\n'
            << "max_thread_abs_error=" << max_thread_error << '\n'
            << "max_formula_abs_error=" << max_formula_error << '\n'
            << "max_ceiling_abs_error=" << max_ceiling_error << '\n'
            << "max_source_abs_error=" << max_source_error << '\n'
            << "max_dopri5_abs_error=" << max_step_error << '\n';

  return max_thread_error == 0 && max_step_error == 0
                 && max_formula_error < Scalar("1e-25")
                 && max_ceiling_error < Scalar("1e-25")
                 && max_source_error < Scalar("1e-25")
                 && interior_ceiling_efficiency >= 0.90
                 && sourced_ceiling_efficiency >= 0.90
             ? EXIT_SUCCESS
             : EXIT_FAILURE;
}

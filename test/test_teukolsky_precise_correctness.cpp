#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
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
using Vector = Equation::Vector;

int failures = 0;
Scalar largest_relative_error = 0;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

Param make_param(long long s, long long l, long long n,
                 const char *r_min, const char *r_max,
                 const char *ko_epsilon) {
  Param param;
  param.s = s;
  param.l = l;
  param.M = Scalar("0.5");
  param.r_min = Scalar(r_min);
  param.r_max = Scalar(r_max);
  param.N = n;
  param.ko_epsilon = Scalar(ko_epsilon);
  param.t_start = 0;
  param.t_end = 1;
  param.t_interval = Scalar("0.1");
  param.delta_t = Scalar("0.001");
  return param;
}

State make_state(long long grid_size, int seed) {
  State state(2 * grid_size);
  for (long long i = 0; i < grid_size; ++i) {
    const long long psi_pattern = ((i * 37 + seed * 11) % 101) - 50;
    const long long pi_pattern = ((i * 29 + seed * 7) % 83) - 41;
    state[i] = Scalar("0.125") + Scalar(psi_pattern) * Scalar("0.00037");
    state[grid_size + i] = Scalar("-0.0625")
                            + Scalar(pi_pattern) * Scalar("0.00029");
  }
  return state;
}

Vector make_source(long long grid_size, int seed) {
  Vector source(grid_size);
  for (long long i = 0; i < grid_size; ++i) {
    source[i] = Scalar(((i * 13 + seed * 5) % 31) - 15) * Scalar("1e-7");
  }
  return source;
}

State reference_rhs(const Equation &equation, const State &x,
                    const Vector &source) {
  const long long grid_size = equation.grid_size;
  const Scalar one_twelfth = Scalar(1) / Scalar(12);
  const Scalar ko_factor = equation.param.ko_epsilon / Scalar(16);
  State result(2 * grid_size);

  const auto psi = x.head(grid_size);
  const auto pi = x.tail(grid_size);
  auto dpsi = result.head(grid_size);
  auto dpi = result.tail(grid_size);

  for (long long i = 2; i <= grid_size - 3; ++i) {
    const Scalar d2 = (-psi(i - 2) + Scalar(16) * psi(i - 1)
                       - Scalar(30) * psi(i) + Scalar(16) * psi(i + 1)
                       - psi(i + 2)) * one_twelfth * equation.inv_h_sqr;
    const Scalar d1 = (psi(i - 2) - Scalar(8) * psi(i - 1)
                       + Scalar(8) * psi(i + 1) - psi(i + 2))
                      * one_twelfth * equation.inv_h;
    dpi(i) = d2 - equation.C[i] * d1 - equation.A[i] * pi(i)
             - equation.V[i] * psi(i) + source(i);
    dpsi(i) = pi(i) - ko_factor * (psi(i - 2) - Scalar(4) * psi(i - 1)
                                    + Scalar(6) * psi(i)
                                    - Scalar(4) * psi(i + 1) + psi(i + 2));
  }

  const Scalar d2_0 = (Scalar(35) * psi(0) - Scalar(104) * psi(1)
                        + Scalar(114) * psi(2) - Scalar(56) * psi(3)
                        + Scalar(11) * psi(4)) * one_twelfth
                       * equation.inv_h_sqr;
  const Scalar d1_0 = (-Scalar(25) * psi(0) + Scalar(48) * psi(1)
                        - Scalar(36) * psi(2) + Scalar(16) * psi(3)
                        - Scalar(3) * psi(4)) * one_twelfth * equation.inv_h;
  dpi(0) = d2_0 - equation.C[0] * d1_0 - equation.A[0] * pi(0)
           - equation.V[0] * psi(0) + source(0);
  dpsi(0) = pi(0) - ko_factor * (psi(0) - Scalar(4) * psi(1)
                                  + Scalar(6) * psi(2) - Scalar(4) * psi(3)
                                  + psi(4));

  const Scalar d2_1 = (Scalar(11) * psi(0) - Scalar(20) * psi(1)
                        + Scalar(6) * psi(2) + Scalar(4) * psi(3) - psi(4))
                       * one_twelfth * equation.inv_h_sqr;
  const Scalar d1_1 = (-Scalar(3) * psi(0) - Scalar(10) * psi(1)
                        + Scalar(18) * psi(2) - Scalar(6) * psi(3) + psi(4))
                       * one_twelfth * equation.inv_h;
  dpi(1) = d2_1 - equation.C[1] * d1_1 - equation.A[1] * pi(1)
           - equation.V[1] * psi(1) + source(1);
  dpsi(1) = pi(1) - ko_factor * (psi(0) - Scalar(4) * psi(1)
                                  + Scalar(6) * psi(2) - Scalar(4) * psi(3)
                                  + psi(4));

  const long long n2 = grid_size - 2;
  const long long n1 = grid_size - 1;
  const Scalar d2_n2 = (-psi(grid_size - 5) + Scalar(4) * psi(grid_size - 4)
                         + Scalar(6) * psi(grid_size - 3)
                         - Scalar(20) * psi(n2) + Scalar(11) * psi(n1))
                        * one_twelfth * equation.inv_h_sqr;
  const Scalar d1_n2 = (-psi(grid_size - 5) + Scalar(6) * psi(grid_size - 4)
                         - Scalar(18) * psi(grid_size - 3)
                         + Scalar(10) * psi(n2) + Scalar(3) * psi(n1))
                        * one_twelfth * equation.inv_h;
  dpi(n2) = d2_n2 - equation.C[n2] * d1_n2 - equation.A[n2] * pi(n2)
            - equation.V[n2] * psi(n2) + source(n2);
  dpsi(n2) = pi(n2) - ko_factor * (psi(grid_size - 5)
                                    - Scalar(4) * psi(grid_size - 4)
                                    + Scalar(6) * psi(grid_size - 3)
                                    - Scalar(4) * psi(n2) + psi(n1));

  const Scalar d2_n1 = (Scalar(11) * psi(grid_size - 5)
                         - Scalar(56) * psi(grid_size - 4)
                         + Scalar(114) * psi(grid_size - 3)
                         - Scalar(104) * psi(n2) + Scalar(35) * psi(n1))
                        * one_twelfth * equation.inv_h_sqr;
  const Scalar d1_n1 = (Scalar(3) * psi(grid_size - 5)
                         - Scalar(16) * psi(grid_size - 4)
                         + Scalar(36) * psi(grid_size - 3)
                         - Scalar(48) * psi(n2) + Scalar(25) * psi(n1))
                        * one_twelfth * equation.inv_h;
  dpi(n1) = d2_n1 - equation.C[n1] * d1_n1 - equation.A[n1] * pi(n1)
            - equation.V[n1] * psi(n1) + source(n1);
  dpsi(n1) = pi(n1) - ko_factor * (psi(grid_size - 5)
                                    - Scalar(4) * psi(grid_size - 4)
                                    + Scalar(6) * psi(grid_size - 3)
                                    - Scalar(4) * psi(n2) + psi(n1));
  return result;
}

Scalar compare_states(const State &actual, const State &expected,
                      const Scalar &relative_tolerance,
                      const std::string &label) {
  require(actual.size() == expected.size(), label + ": state size");
  Scalar max_relative_error = 0;
  for (Eigen::Index i = 0; i < actual.size(); ++i) {
    const Scalar scale = std::max(Scalar(1), abs(expected[i]));
    const Scalar relative_error = abs(actual[i] - expected[i]) / scale;
    max_relative_error = std::max(max_relative_error, relative_error);
    if (relative_error > relative_tolerance) {
      std::cerr << "FAIL: " << label << " index=" << i
                << " relative_error=" << relative_error << '\n';
      ++failures;
      break;
    }
  }
  largest_relative_error = std::max(largest_relative_error, max_relative_error);
  return max_relative_error;
}

void check_rhs_case(const Param &param, int seed) {
  Equation equation(param);
  const long long grid_size = equation.grid_size;
  const State state = make_state(grid_size, seed);
  const Vector zero = Vector::Zero(grid_size);
  const Vector source = make_source(grid_size, seed);
  const Scalar tolerance("2e-27");
  const std::string prefix = "s=" + std::to_string(param.s)
                             + ",l=" + std::to_string(param.l)
                             + ",N=" + std::to_string(param.N);

  State one_thread(2 * grid_size);
  omp_set_num_threads(1);
  equation(state, one_thread, Scalar("0.375"));
  compare_states(one_thread, reference_rhs(equation, state, zero), tolerance,
                 prefix + " homogeneous reference");

  for (int threads : {2, 6}) {
    State threaded(2 * grid_size);
    omp_set_num_threads(threads);
    equation(state, threaded, Scalar("0.375"));
    compare_states(threaded, one_thread, Scalar(0),
                   prefix + " homogeneous threads=" + std::to_string(threads));
  }

  equation.Q = [source](const Scalar &, Vector &output) { output = source; };
  State generic_source(2 * grid_size);
  omp_set_num_threads(6);
  equation(state, generic_source, Scalar("0.375"));
  compare_states(generic_source, reference_rhs(equation, state, source), tolerance,
                 prefix + " generic source reference");

  const Scalar source_scale("-0.625");
  Vector spatial = source / source_scale;
  equation.set_separable_source(
      spatial, [source_scale](const Scalar &) { return source_scale; });
  State separable_source(2 * grid_size);
  equation(state, separable_source, Scalar("0.375"));
  compare_states(separable_source, generic_source, Scalar(0),
                 prefix + " generic versus separable source");

  bool threw = false;
  try {
    equation.set_separable_source(Vector::Zero(grid_size - 1),
                                  [](const Scalar &) { return Scalar(1); });
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  require(threw, prefix + " rejects a mismatched separable source");
}

struct ReferenceSystem {
  const Equation &equation;
  const Vector &source;

  void operator()(const State &state, State &derivative, const Scalar &) const {
    derivative = reference_rhs(equation, state, source);
  }
};

void check_dopri5_trajectory() {
  using Stepper = boost::numeric::odeint::runge_kutta_dopri5<
      State, Scalar, State, Scalar>;
  const Param param = make_param(-2, 3, 79, "-35", "70", "0.7");
  Equation equation(param);
  const Vector source = make_source(equation.grid_size, 17);
  equation.Q = [source](const Scalar &, Vector &output) { output = source; };
  ReferenceSystem reference_system{equation, source};

  State optimized_state = make_state(equation.grid_size, 17);
  State reference_state = optimized_state;
  Stepper optimized_stepper;
  Stepper reference_stepper;
  const Scalar dt("0.0005");
  Scalar time = 0;
  omp_set_num_threads(6);
  for (int step = 0; step < 25; ++step) {
    optimized_stepper.do_step(std::ref(equation), optimized_state, time, dt);
    reference_stepper.do_step(std::ref(reference_system), reference_state, time, dt);
    time += dt;
  }
  compare_states(optimized_state, reference_state, Scalar("2e-24"),
                 "25-step Dopri5 trajectory");
}

}  // namespace

int main() {
  omp_set_dynamic(0);
  check_rhs_case(make_param(-2, 2, 31, "-30", "60", "0.7"), 1);
  check_rhs_case(make_param(-1, 3, 67, "-45", "90", "0.5"), 3);
  check_rhs_case(make_param(0, 0, 96, "-20", "45", "0.0"), 5);
  check_rhs_case(make_param(1, 2, 127, "-25", "55", "0.3"), 7);
  check_dopri5_trajectory();

  if (failures != 0) {
    std::cerr << failures << " precise Teukolsky correctness checks failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: precise Teukolsky correctness matrix; largest relative error="
            << largest_relative_error << '\n';
  return EXIT_SUCCESS;
}

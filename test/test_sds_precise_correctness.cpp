#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <omp.h>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
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

int failures = 0;
Scalar largest_relative_error = 0;
Scalar largest_horizon_formula_difference = 0;
Scalar largest_legacy_profile_difference = 0;

void require(const bool condition, const std::string &message) {
  if(!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

Param make_param(const long long int s, const long long int l,
                 const long long int n, const char *x_min,
                 const char *x_max, const char *lambda = "0.01") {
  Param param;
  param.s = s;
  param.l = l;
  param.M = Scalar("0.5");
  param.Lambda = Scalar(lambda);
  param.r_min = Scalar(x_min);
  param.r_max = Scalar(x_max);
  param.N = n;
  param.t_start = 0;
  param.t_end = 1;
  param.t_interval = Scalar("0.1");
  param.delta_t = Scalar("0.001");
  return param;
}

template<typename HP, typename Function>
HP bisect_root(Function function, HP lower, HP upper) {
  HP lower_value = function(lower);
  require(lower_value * function(upper) <= 0, "root-finding bracket");
  for(int iteration = 0; iteration < 500; ++iteration) {
    const HP middle = (lower + upper) / HP(2);
    const HP middle_value = function(middle);
    if(middle_value == 0) return middle;
    if((lower_value < 0) == (middle_value < 0)) {
      lower = middle;
      lower_value = middle_value;
    } else {
      upper = middle;
    }
  }
  return (lower + upper) / HP(2);
}

template<typename HP>
std::array<HP, 3> analytic_horizons(const HP &mass, const HP &lambda) {
  const HP pi = boost::math::constants::pi<HP>();
  const HP sqrt_lambda = sqrt(lambda);
  const HP angle = acos(HP(3) * mass * sqrt_lambda) / HP(3);
  const HP rb = HP(2) / sqrt_lambda * cos(angle + pi / HP(3));
  const HP rc = HP(2) / sqrt_lambda * cos(angle - pi / HP(3));
  return {rb, rc, -(rb + rc)};
}

template<typename HP>
std::array<HP, 3> root_found_horizons(const HP &mass, const HP &lambda) {
  auto cubic = [&](const HP &radius) {
    return lambda * radius * radius * radius - HP(3) * radius
           + HP(6) * mass;
  };
  const HP rb = bisect_root(cubic, HP(2) * mass, HP(3) * mass);
  const HP rc = bisect_root(cubic, HP(3) * mass, sqrt(HP(3) / lambda));
  const HP ro = bisect_root(cubic, -HP(2) * sqrt(HP(3) / lambda), HP(0));
  return {rb, rc, ro};
}

State make_state(const long long int grid_size, const int seed) {
  State state(2 * grid_size);
  for(long long int i = 0; i < grid_size; ++i) {
    state[i] = Scalar("0.125")
               + Scalar(((i * 37 + seed * 11) % 101) - 50)
                     * Scalar("0.00037");
    state[grid_size + i] = Scalar("-0.0625")
                           + Scalar(((i * 29 + seed * 7) % 83) - 41)
                                 * Scalar("0.00029");
  }
  return state;
}

Vector make_source(const long long int grid_size, const int seed) {
  Vector source(grid_size);
  for(long long int i = 0; i < grid_size; ++i) {
    source[i] = Scalar(((i * 13 + seed * 5) % 31) - 15)
                * Scalar("1e-7");
  }
  return source;
}

State reference_rhs(const Equation &equation, const State &state,
                    const Vector &source) {
  const long long int grid_size = equation.grid_size;
  const Scalar d2_factor = equation.inv_h_sqr / Scalar(12);
  const Scalar d1_factor = equation.inv_h / Scalar(12);
  State result(2 * grid_size);
  const auto psi = state.head(grid_size);
  const auto pi = state.tail(grid_size);
  auto dpsi = result.head(grid_size);
  auto dpi = result.tail(grid_size);

  for(long long int i = 0; i < grid_size; ++i) dpsi(i) = pi(i);
  for(long long int i = 2; i <= grid_size - 3; ++i) {
    dpi(i) = (-psi(i - 2) + Scalar(16) * psi(i - 1)
              - Scalar(30) * psi(i) + Scalar(16) * psi(i + 1)
              - psi(i + 2)) * d2_factor
             - equation.V[i] * psi(i) + source(i);
  }

  dpi(0) = (-Scalar(25) * pi(0) + Scalar(48) * pi(1)
            - Scalar(36) * pi(2) + Scalar(16) * pi(3)
            - Scalar(3) * pi(4)) * d1_factor
           - equation.V[0] * psi(0) + source(0);
  dpi(1) = (Scalar(11) * psi(0) - Scalar(20) * psi(1)
            + Scalar(6) * psi(2) + Scalar(4) * psi(3) - psi(4))
           * d2_factor - equation.V[1] * psi(1) + source(1);

  const long long int n2 = grid_size - 2;
  const long long int n1 = grid_size - 1;
  dpi(n2) = (-psi(grid_size - 5) + Scalar(4) * psi(grid_size - 4)
             + Scalar(6) * psi(grid_size - 3) - Scalar(20) * psi(n2)
             + Scalar(11) * psi(n1)) * d2_factor
            - equation.V[n2] * psi(n2) + source(n2);
  dpi(n1) = (-Scalar(3) * pi(grid_size - 5)
             + Scalar(16) * pi(grid_size - 4)
             - Scalar(36) * pi(grid_size - 3) + Scalar(48) * pi(n2)
             - Scalar(25) * pi(n1)) * d1_factor
            - equation.V[n1] * psi(n1) + source(n1);
  return result;
}

Scalar compare_states(const State &actual, const State &expected,
                      const Scalar &relative_tolerance,
                      const std::string &label) {
  require(actual.size() == expected.size(), label + ": size");
  Scalar max_error = 0;
  for(Eigen::Index i = 0; i < actual.size(); ++i) {
    const Scalar scale = std::max(Scalar(1), abs(expected[i]));
    const Scalar error = abs(actual[i] - expected[i]) / scale;
    max_error = std::max(max_error, error);
    if(error > relative_tolerance) {
      std::cerr << "FAIL: " << label << " index=" << i
                << " relative_error=" << error << '\n';
      ++failures;
      break;
    }
  }
  largest_relative_error = std::max(largest_relative_error, max_error);
  return max_error;
}

void check_parameter_validation() {
  auto rejects = [](const Param &param) {
    try {
      Equation equation(param);
      (void)equation;
      return false;
    } catch(const std::invalid_argument &) {
      return true;
    }
  };

  Param param = make_param(0, 1, 12, "-10", "20");
  param.M = 0;
  require(rejects(param), "rejects M <= 0");
  param = make_param(0, 1, 12, "-10", "20");
  param.Lambda = 0;
  require(rejects(param), "rejects Lambda <= 0");
  param = make_param(0, 1, 12, "-10", "20");
  param.Lambda = Scalar("0.45");
  require(rejects(param), "rejects extremal or superextremal geometry");
  param = make_param(3, 3, 12, "-10", "20");
  require(rejects(param), "rejects unsupported spin");
  param = make_param(2, 1, 12, "-10", "20");
  require(rejects(param), "rejects l < s");
  param = make_param(0, 1, 3, "-10", "20");
  require(rejects(param), "rejects undersized grid");
}

void check_horizon_root_algorithms() {
  using HP = boost::multiprecision::cpp_bin_float_100;
  for(const char *lambda_text : {
          "1e-4", "1e-12", "1e-100", "0.44",
          "0.444444444444444444444444444444"}) {
    const HP mass("0.5");
    const HP lambda(lambda_text);
    const auto analytic = analytic_horizons(mass, lambda);
    const auto root_found = root_found_horizons(mass, lambda);
    for(std::size_t root = 0; root < analytic.size(); ++root) {
      const HP relative_difference = abs(analytic[root] - root_found[root])
                                     / abs(root_found[root]);
      largest_horizon_formula_difference = std::max(
          largest_horizon_formula_difference,
          relative_difference.convert_to<Scalar>());
      require(relative_difference < HP("1e-45"),
              "analytic and root-found horizons agree beyond binary128");
    }

    Param param = make_param(0, 1, 4, "-2", "2", lambda_text);
    Equation equation(param);
    const auto production_reference = root_found_horizons(
        HP(param.M), HP(param.Lambda));
    const std::array<Scalar, 3> production = {
        equation.r_black_hole, equation.r_cosmological, equation.r_negative};
    for(std::size_t root = 0; root < production.size(); ++root) {
      const Scalar reference = production_reference[root].convert_to<Scalar>();
      const Scalar relative_difference = abs(production[root] - reference)
                                         / abs(reference);
      if(relative_difference >= Scalar("2e-32")) {
        std::cerr << "horizon mismatch: Lambda=" << lambda_text
                  << " root=" << root
                  << " relative_difference=" << relative_difference << '\n';
      }
      require(relative_difference < Scalar("2e-32"),
              "production horizons agree with independent root finding");
    }
  }
}

void check_geometry_case(const Param &param, const std::string &label) {
  Equation equation(param);
  require(equation.r_negative < 0 && equation.r_black_hole > 0
          && equation.r_black_hole < equation.r_cosmological,
          label + ": ordered horizons");
  require(equation.kappa_black_hole > 0 && equation.kappa_cosmological > 0,
          label + ": positive surface gravities");

  for(long long int i = 0; i < equation.grid_size; ++i) {
    require(equation.r[i] >= equation.r_black_hole
            && equation.r[i] < equation.r_cosmological
            && equation.rho_cosmological[i] > 0 && equation.f[i] > 0,
            label + ": geometry lies between the horizons");
    if(i > 0) {
      require(equation.r[i] >= equation.r[i - 1],
              label + ": monotonic r(x)");
      require(equation.grid_coordinate(i)
                  > equation.grid_coordinate(i - 1),
              label + ": monotonic x grid");
    }
  }
}

void check_neighbor_reuse_consistency() {
  const Param param = make_param(2, 2, 600, "-500", "1000", "1e-4");
  omp_set_num_threads(1);
  const Equation serial(param);
  omp_set_num_threads(6);
  const Equation blocked(param);
  const Scalar h = serial.grid_space();
  for(long long int i = 0; i < serial.grid_size; ++i) {
    require(serial.grid_coordinate(i)
                == Equation::grid_coordinate(param.r_min, h, i),
            "shared tortoise-grid coordinate helper");
    require(serial.r[i] == blocked.r[i]
            && serial.rho_cosmological[i] == blocked.rho_cosmological[i]
            && serial.f[i] == blocked.f[i] && serial.V[i] == blocked.V[i],
            "neighbor inversion is independent of OpenMP block boundaries");
  }

  // More workers than points creates empty blocks.  In particular, the last
  // nonempty block must still write index N, but no block may write N + 1.
  const Param minimal_param = make_param(0, 1, 4, "-2", "2", "1e-3");
  omp_set_num_threads(1);
  const Equation minimal_serial(minimal_param);
  omp_set_num_threads(16);
  const Equation minimal_blocked(minimal_param);
  require(minimal_blocked.grid_size == minimal_param.N + 1,
          "rightmost tortoise-grid index is N");
  for(long long int i = 0; i <= minimal_param.N; ++i) {
    require(minimal_serial.r[i] == minimal_blocked.r[i]
            && minimal_serial.rho_cosmological[i]
                   == minimal_blocked.rho_cosmological[i]
            && minimal_serial.f[i] == minimal_blocked.f[i]
            && minimal_serial.V[i] == minimal_blocked.V[i],
            "empty OpenMP blocks preserve every grid point through index N");
  }
}

void check_schwarzschild_limit() {
  using HP = boost::multiprecision::cpp_bin_float_100;
  Scalar previous_error = std::numeric_limits<Scalar>::infinity();
  for(const char *lambda : {"1e-8", "1e-12", "1e-16"}) {
    Equation equation(make_param(0, 1, 40, "-20", "40", lambda));
    Scalar max_coordinate_error = 0;
    for(long long int i : {10LL, 20LL, 30LL, 40LL}) {
      const HP x(equation.grid_coordinate(i));
      const HP radius(equation.r[i]);
      const HP schwarzschild_x = radius + log(radius - HP(1));
      const Scalar coordinate_error =
          (abs(schwarzschild_x - x) / (HP(1) + abs(x)))
              .convert_to<Scalar>();
      max_coordinate_error = std::max(max_coordinate_error,
                                      coordinate_error);

      const HP inverse_radius = HP(1) + boost::math::lambert_w0(
          exp(x - HP(1)));
      const Scalar inverse_error =
          (abs(radius - inverse_radius) / inverse_radius).convert_to<Scalar>();
      require(inverse_error < Scalar("3e-4"),
              "SdS radius approaches the Schwarzschild inverse");
    }
    require(max_coordinate_error < previous_error / Scalar(1000),
            "Schwarzschild-coordinate error decreases with Lambda");
    previous_error = max_coordinate_error;
  }
  require(previous_error < Scalar("1e-9"),
          "r_* = r + log(r-1) in the Lambda -> 0 limit for M=0.5");
}

#ifndef SDS_SKIP_2000_DIGIT_REFERENCE
template<typename HP>
struct ReferencePoint {
  HP r;
  HP rho_b;
  HP rho_c;
  HP f;
};

template<typename HP>
ReferencePoint<HP> reference_point(const Param &param, const HP &target_x) {
  const HP mass(param.M);
  const HP lambda(param.Lambda);
  const HP pi = boost::math::constants::pi<HP>();
  const HP sqrt_lambda = sqrt(lambda);
  const HP angle = acos(HP(3) * mass * sqrt_lambda) / HP(3);
  const HP rb = HP(2) / sqrt_lambda * cos(angle + pi / HP(3));
  const HP rc = HP(2) / sqrt_lambda * cos(angle - pi / HP(3));
  const HP ro = -(rb + rc);
  const HP delta = rc - rb;
  const HP reference_r = HP(3) * mass;
  const HP reference_x = reference_r + HP(2) * mass * log(HP("0.5"));
  auto f_prime = [&](const HP &radius) {
    return HP(2) * mass / (radius * radius) - HP(2) * lambda * radius / HP(3);
  };
  const HP inv_b = HP(1) / f_prime(rb);
  const HP inv_c = HP(1) / f_prime(rc);
  const HP inv_o = HP(1) / f_prime(ro);

  auto evaluate = [&](const HP &y) {
    HP rho_b;
    HP rho_c;
    if(y >= 0) {
      const HP exponential = exp(-y);
      rho_b = delta / (HP(1) + exponential);
      rho_c = delta * exponential / (HP(1) + exponential);
    } else {
      const HP exponential = exp(y);
      rho_b = delta * exponential / (HP(1) + exponential);
      rho_c = delta / (HP(1) + exponential);
    }
    const HP radius = rb + rho_b;
    const HP radius_minus_ro = rb - ro + rho_b;
    const HP x = reference_x
        + inv_b * (log(rho_b) - log(reference_r - rb))
        + inv_c * (log(rho_c) - log(rc - reference_r))
        + inv_o * (log(radius_minus_ro) - log(reference_r - ro));
    return std::pair<HP, ReferencePoint<HP>>(
        x, {radius, rho_b, rho_c,
            lambda * rho_b * rho_c * radius_minus_ro / (HP(3) * radius)});
  };

  HP lower(-4096);
  HP upper(4096);
  require(evaluate(lower).first < target_x && evaluate(upper).first > target_x,
          "2000-decimal reference bracket");
  HP y = 0;
  for(int iteration = 0; iteration < 80; ++iteration) {
    const auto evaluation = evaluate(y);
    const HP residual = evaluation.first - target_x;
    if(abs(residual) < HP("1e-80") * (HP(1) + abs(target_x))) break;
    if(residual < 0) lower = y;
    else upper = y;
    const HP derivative = HP(3) * evaluation.second.r
        / (lambda * delta * (evaluation.second.r - ro));
    const HP newton = y - residual / derivative;
    y = (newton > lower && newton < upper)
        ? newton : (lower + upper) / HP(2);
  }
  return evaluate(y).second;
}

void check_2000_digit_reference() {
  using HP2000 = boost::multiprecision::number<
      boost::multiprecision::cpp_bin_float<2000>>;
  const Param param = make_param(2, 2, 8, "-500", "1000", "0.01");
  Equation equation(param);
  for(long long int i : {0LL, 4LL, 8LL}) {
    const ReferencePoint<HP2000> reference =
        reference_point<HP2000>(param, HP2000(equation.grid_coordinate(i)));
    const Scalar f_reference = reference.f.convert_to<Scalar>();
    const Scalar inv_r = Scalar(1) / reference.r.convert_to<Scalar>();
    const Scalar potential_reference = f_reference
        * (Scalar(param.l * (param.l + 1)) * inv_r * inv_r
           + Scalar(1 - param.s * param.s) * Scalar(2) * param.M
                 * inv_r * inv_r * inv_r);
    require(abs(equation.f[i] - f_reference)
                / std::max(abs(f_reference), Scalar("1e-4900"))
                < Scalar("2e-31"),
            "f agrees with 2000-decimal reference");
    require(abs(equation.V[i] - potential_reference)
                / std::max(abs(potential_reference), Scalar("1e-4900"))
                < Scalar("5e-31"),
            "V agrees with 2000-decimal reference");
  }

  require(equation.f[0] > 0, "near-horizon f remains nonzero");
  const Scalar naive_f = Scalar(1) - Scalar(2) * param.M / equation.r[0]
                         - param.Lambda * equation.r[0] * equation.r[0]
                               / Scalar(3);
  require(abs(naive_f - equation.f[0]) > Scalar("1e20") * equation.f[0],
          "factorized f avoids near-horizon cancellation");
}
#endif

void check_multipole_constructor() {
  Equation l2(make_param(2, 2, 32, "-40", "80"));
  Equation l3(make_param(2, 3, 32, "-40", "80"));
  for(long long int i = 0; i < l2.grid_size; ++i) {
    require(l2.r[i] == l3.r[i] && l2.f[i] == l3.f[i],
            "multipole leaves geometry unchanged");
    const Scalar expected_difference = l2.f[i] * Scalar(6)
                                       / (l2.r[i] * l2.r[i]);
    require(abs((l3.V[i] - l2.V[i]) - expected_difference)
                < Scalar("2e-31") * std::max(Scalar(1), abs(l3.V[i])),
            "constructor builds V from input multipole");
  }
}

void check_source_profiles() {
  Equation equation(make_param(0, 1, 120, "-80", "1000"));
  SdSTranslatedSourceParam source;
  source.beta = Scalar("2.5");
  source.cutoff_sigma = 10;
  source.onset_time = Scalar("-1e6");

  auto profile_at = [&](const SdSTranslatedSourceParam &source_param,
                        const long long int i) {
    SdSSource Q(source_param);
    Q.initialize(equation.param.r_min, equation.grid_space(),
                 equation.grid_size, equation.r_cosmological,
                 equation.r, equation.rho_cosmological, equation.f);
    Vector values(equation.grid_size);
    Q(equation.grid_coordinate(i) + source_param.u_center, values);
    const Scalar pi = boost::math::constants::pi<Scalar>();
    const Scalar normalization = source_param.amplitude
        / (sqrt(Scalar(2) * pi) * source_param.sigma);
    return values[i] / normalization;
  };

  source.profile = SdSSourceProfile::ArealPower;
  for(long long int i : {0LL, 60LL, 120LL}) {
    require(abs(profile_at(source, i) - pow(equation.r[i], -source.beta))
                < Scalar("2e-31"),
            "areal-power source profile");
  }

  source.profile = SdSSourceProfile::HorizonSubtractedArealPower;
  const long long int right = equation.grid_size - 1;
  const Scalar subtracted_right = profile_at(source, right);
  require(subtracted_right > 0,
          "subtracted source remains nonzero near cosmological horizon");
  const Scalar leading = source.beta
      * pow(equation.r_cosmological, -source.beta - Scalar(1))
      * equation.rho_cosmological[right];
  require(abs(subtracted_right / leading - Scalar(1))
              < Scalar("1e-25"),
          "subtracted source has stable horizon expansion");

  source.profile = SdSSourceProfile::LocalScalar;
  for(long long int i : {0LL, 60LL, 120LL}) {
    const Scalar expected = equation.f[i] * pow(equation.r[i], -source.beta);
    require(abs(profile_at(source, i) - expected)
                < Scalar("2e-31") * std::max(Scalar(1), abs(expected)),
            "local scalar source profile");
  }

  source.profile = SdSSourceProfile::TortoisePower;
  source.L = 3;
  source.x0 = 10;
  source.X0 = 0;
  source.X1 = 20;
  for(long long int i = 0; i < equation.grid_size; ++i) {
    const Scalar x = equation.grid_coordinate(i);
    if(x <= source.X0) {
      require(profile_at(source, i) == 0,
              "tortoise profile vanishes below cutoff");
    }
    if(x >= source.X1) {
      const Scalar expected = pow(source.L / (x + source.x0), source.beta);
      require(abs(profile_at(source, i) - expected)
                  < Scalar("2e-31") * std::max(Scalar(1), abs(expected)),
              "tortoise profile above cutoff");
    }
  }

  source.beta = 0;
  source.profile = SdSSourceProfile::ArealPower;
  require(profile_at(source, 60) == 1,
          "beta zero gives constant areal-power profile");
  source.profile = SdSSourceProfile::HorizonSubtractedArealPower;
  require(profile_at(source, right) == 0,
          "beta zero gives zero horizon-subtracted profile");
  source.profile = SdSSourceProfile::LocalScalar;
  require(profile_at(source, 60) == equation.f[60],
          "beta zero gives metric-function local profile");
  source.profile = SdSSourceProfile::TortoisePower;
  require(profile_at(source, 120) == 1,
          "beta zero gives unit tortoise profile above cutoff");

  bool rejected = false;
  try {
    source.X1 = source.X0;
    const SdSSource invalid(source);
  } catch(const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "rejects invalid tortoise source cutoff");
}

void check_rhs_and_sources() {
  const Param param = make_param(1, 2, 79, "-35", "70");
  Equation equation(param);
  const State state = make_state(equation.grid_size, 7);
  const Vector zero = Vector::Zero(equation.grid_size);
  const Vector source = make_source(equation.grid_size, 7);

  State homogeneous(2 * equation.grid_size);
  omp_set_num_threads(1);
  equation(state, homogeneous, Scalar("0.375"));
  compare_states(homogeneous, reference_rhs(equation, state, zero),
                 Scalar("2e-27"), "homogeneous original-formula RHS");

  for(int threads : {2, 6}) {
    State threaded(2 * equation.grid_size);
    omp_set_num_threads(threads);
    equation(state, threaded, Scalar("0.375"));
    compare_states(threaded, homogeneous, Scalar(0),
                   "homogeneous thread reproducibility");
  }

  Equation generic_equation(param, SdSSource(SdSSource::GenericSource(
      [source](const Scalar &, Vector &result) { result = source; })));
  State generic(2 * equation.grid_size);
  generic_equation(state, generic, Scalar("0.375"));
  compare_states(generic, reference_rhs(generic_equation, state, source),
                 Scalar("2e-27"), "generic source RHS");

  Equation wrong_size_equation(param, SdSSource(SdSSource::GenericSource(
      [&](const Scalar &, Vector &result) {
        result.resize(equation.grid_size - 1);
      })));
  bool rejected_source_size = false;
  try {
    wrong_size_equation(state, generic, Scalar("0.375"));
  } catch(const std::invalid_argument &) {
    rejected_source_size = true;
  }
  require(rejected_source_size, "rejects generic source with wrong size");

  SdSTranslatedSourceParam translated;
  translated.profile = SdSSourceProfile::ArealPower;
  translated.beta = 2;
  translated.u_center = -10;
  translated.sigma = 2;
  translated.cutoff_sigma = 8;
  translated.onset_time = Scalar("0.1");
  Equation translated_equation(param, SdSSource(translated));

  const Scalar time("0.375");
  Vector translated_values(equation.grid_size);
  translated_equation.Q(time, translated_values);
  State translated_rhs(2 * equation.grid_size);
  translated_equation(state, translated_rhs, time);
  compare_states(translated_rhs,
                 reference_rhs(translated_equation, state, translated_values),
                 Scalar("2e-27"), "translated Gaussian source RHS");

  State before_onset(2 * equation.grid_size);
  translated_equation(state, before_onset, Scalar("0.05"));
  compare_states(before_onset, homogeneous, Scalar(0),
                 "translated source turn-on");

  translated.waveform = SdSWaveform::GaussianDerivative;
  Equation derivative_equation(param, SdSSource(translated));
  const long long int i = equation.grid_size / 2;
  const Scalar center_time = equation.grid_coordinate(i)
                             + translated.u_center;
  Vector center_values(equation.grid_size);
  Vector left_values(equation.grid_size);
  Vector right_values(equation.grid_size);
  derivative_equation.Q(center_time, center_values);
  derivative_equation.Q(center_time - translated.sigma, left_values);
  derivative_equation.Q(center_time + translated.sigma, right_values);
  require(center_values[i] == 0,
          "Gaussian derivative vanishes at its center");
  require(abs(left_values[i] + right_values[i]) < Scalar("1e-31"),
          "Gaussian derivative is zero-mean and antisymmetric");

  SdSSpacetimeGaussianSourceParam green;
  green.amplitude = Scalar("1.75");
  green.r_ast_center = equation.grid_coordinate(i);
  green.t_center = Scalar("0.625");
  green.sigma = Scalar("1.25");
  green.cutoff_sigma = 10;
  Equation green_equation(param, SdSSource(green));
  Vector green_values(equation.grid_size);
  green_equation.Q(green.t_center, green_values);
  const Scalar pi = boost::math::constants::pi<Scalar>();
  const Scalar normalization = green.amplitude
      / (sqrt(Scalar(2) * pi) * green.sigma);
  require(abs(green_values[i] - normalization) < Scalar("2e-31"),
          "spacetime Gaussian has requested normalization at its center");
  const long long int j = i + 3;
  const Scalar dx = equation.grid_coordinate(j) - green.r_ast_center;
  const Scalar expected_spatial = normalization
      * exp(-dx * dx / (green.sigma * green.sigma));
  require(abs(green_values[j] - expected_spatial) < Scalar("2e-31"),
          "spacetime Gaussian has exp(-dx^2/sigma^2) spatial profile");
  green_equation.Q(green.t_center + green.sigma, green_values);
  require(abs(green_values[i] - normalization * exp(Scalar(-1)))
              < Scalar("2e-31"),
          "spacetime Gaussian has exp(-dt^2/sigma^2) time profile");
  green_equation.Q(green.t_center, green_values);
  State green_rhs(2 * equation.grid_size);
  green_equation(state, green_rhs, green.t_center);
  compare_states(green_rhs,
                 reference_rhs(green_equation, state, green_values),
                 Scalar("2e-27"), "spacetime Gaussian source RHS");
  green_equation.Q(green.t_center
                       + Scalar("1.01") * green.cutoff_sigma * green.sigma,
                   green_values);
  require(green_values.abs().maxCoeff() == 0,
          "spacetime Gaussian vanishes beyond its numerical cutoff");
}

struct ReferenceSystem {
  const Equation &equation;
  const Vector &source;
  void operator()(const State &state, State &derivative, const Scalar &) const {
    derivative = reference_rhs(equation, state, source);
  }
};

struct LegacySourcedSystem {
  const Equation &equation;
  mutable Vector source;

  explicit LegacySourcedSystem(const Equation &equation_)
      : equation(equation_), source(Vector::Zero(equation_.grid_size)) {}

  void operator()(const State &state, State &derivative,
                  const Scalar &time) const {
    const long long int grid_size = equation.grid_size;
    const Scalar d2_factor = equation.inv_h_sqr / Scalar(12);
    const Scalar near_factor = Scalar(16) * d2_factor;
    const Scalar far_factor = -d2_factor;
    const Scalar center_coefficient = -Scalar(30) * d2_factor;
    const Scalar d1_factor = equation.inv_h / Scalar(12);
    equation.Q(time, source);

    const Scalar *__restrict__ psi = state.data();
    const Scalar *__restrict__ pi = state.data() + grid_size;
    Scalar *__restrict__ dpsi = derivative.data();
    Scalar *__restrict__ dpi = derivative.data() + grid_size;
    const Scalar *__restrict__ potential = equation.V.data();
    const Scalar *__restrict__ source_data = source.data();

#pragma omp parallel for schedule(static)
    for(long long int i = 2; i <= grid_size - 3; ++i) {
      const Scalar near_sum = psi[i - 1] + psi[i + 1];
      const Scalar far_sum = psi[i - 2] + psi[i + 2];
      dpi[i] = near_factor * near_sum + far_factor * far_sum
               + (center_coefficient - potential[i]) * psi[i]
               + source_data[i];
      dpsi[i] = pi[i];
    }

    dpi[0] = (-Scalar(25) * pi[0] + Scalar(48) * pi[1]
              - Scalar(36) * pi[2] + Scalar(16) * pi[3]
              - Scalar(3) * pi[4]) * d1_factor
             - potential[0] * psi[0] + source_data[0];
    dpsi[0] = pi[0];
    dpi[1] = (Scalar(11) * psi[0] - Scalar(20) * psi[1]
              + Scalar(6) * psi[2] + Scalar(4) * psi[3] - psi[4])
             * d2_factor - potential[1] * psi[1] + source_data[1];
    dpsi[1] = pi[1];

    const long long int n2 = grid_size - 2;
    const long long int n1 = grid_size - 1;
    dpi[n2] = (-psi[grid_size - 5] + Scalar(4) * psi[grid_size - 4]
               + Scalar(6) * psi[grid_size - 3] - Scalar(20) * psi[n2]
               + Scalar(11) * psi[n1]) * d2_factor
              - potential[n2] * psi[n2] + source_data[n2];
    dpsi[n2] = pi[n2];
    dpi[n1] = (-Scalar(3) * pi[grid_size - 5]
               + Scalar(16) * pi[grid_size - 4]
               - Scalar(36) * pi[grid_size - 3]
               + Scalar(48) * pi[n2] - Scalar(25) * pi[n1]) * d1_factor
              - potential[n1] * psi[n1] + source_data[n1];
    dpsi[n1] = pi[n1];
  }
};

State make_random_state(const long long int grid_size, const std::uint64_t seed) {
  std::mt19937_64 generator(seed);
  std::uniform_int_distribution<long long int> distribution(-1000000, 1000000);
  State state(2 * grid_size);
  for(Eigen::Index i = 0; i < state.size(); ++i) {
    state[i] = Scalar(distribution(generator)) * Scalar("1e-9");
  }
  return state;
}

void check_randomized_legacy_profile_evolution() {
  using Stepper = boost::numeric::odeint::runge_kutta_dopri5<
      State, Scalar, State, Scalar>;
  using LegacyOperations =
      boost::numeric::odeint::eigen_operations<State, false>;
  using LegacyStepper = boost::numeric::odeint::runge_kutta_dopri5<
      State, Scalar, State, Scalar,
      boost::numeric::odeint::vector_space_algebra, LegacyOperations>;
  constexpr std::array<const char *, 3> lambdas = {
      "0.0444444444444444444444444444444444",
      "0.177777777777777777777777777777778",
      "0.355555555555555555555555555555556"};

  for(int sample = 0; sample < 3; ++sample) {
    const Param param = make_param(0, sample, 256, "-30", "60",
                                   lambdas[sample]);
    SdSTranslatedSourceParam source_param;
    source_param.profile = SdSSourceProfile::ArealPower;
    source_param.waveform = SdSWaveform::Gaussian;
    source_param.beta = Scalar(sample);
    source_param.amplitude = 1;
    source_param.u_center = -10;
    source_param.sigma = Scalar("0.5");
    source_param.onset_time = 0;
    source_param.cutoff_sigma = 12;
    Equation equation(param, SdSSource(source_param));
    LegacySourcedSystem legacy(equation);

    State optimized_state = make_random_state(
        equation.grid_size, 0x5d5f0000ULL + sample);
    State legacy_state = optimized_state;
    State optimized_rhs(2 * equation.grid_size);
    State legacy_rhs(2 * equation.grid_size);
    const Scalar start_time = Scalar("3.125") + Scalar(sample) / Scalar(7);
    equation(optimized_state, optimized_rhs, start_time);
    legacy(legacy_state, legacy_rhs, start_time);
    largest_legacy_profile_difference = std::max(
        largest_legacy_profile_difference,
        compare_states(optimized_rhs, legacy_rhs, Scalar(0),
                       "randomized legacy-equivalent full RHS"));

    Stepper optimized_stepper;
    LegacyStepper legacy_stepper;
    const Scalar dt("0.0025");
    Scalar time = start_time;
    for(int step = 0; step < 40; ++step) {
      optimized_stepper.do_step(
          std::ref(equation), optimized_state, time, dt);
      legacy_stepper.do_step(std::ref(legacy), legacy_state, time, dt);
      time += dt;
      largest_legacy_profile_difference = std::max(
          largest_legacy_profile_difference,
          compare_states(
              optimized_state, legacy_state, Scalar(0),
              "randomized legacy-equivalent full profile evolution"));
    }
  }
}

void check_dopri5_trajectory() {
  using Stepper = boost::numeric::odeint::runge_kutta_dopri5<
      State, Scalar, State, Scalar>;
  const Param param = make_param(2, 3, 63, "-30", "60");
  const Vector source = make_source(param.N + 1, 17);
  Equation equation(param, SdSSource(SdSSource::GenericSource(
      [source](const Scalar &, Vector &result) { result = source; })));
  const ReferenceSystem reference{equation, source};

  State optimized_state = make_state(equation.grid_size, 17);
  State reference_state = optimized_state;
  Stepper optimized_stepper;
  Stepper reference_stepper;
  const Scalar dt("0.0005");
  Scalar time = 0;
  omp_set_num_threads(6);
  for(int step = 0; step < 25; ++step) {
    optimized_stepper.do_step(std::ref(equation), optimized_state, time, dt);
    reference_stepper.do_step(std::ref(reference), reference_state, time, dt);
    time += dt;
  }
  compare_states(optimized_state, reference_state, Scalar("2e-24"),
                 "25-step Dopri5 trajectory");
}

}  // namespace

int main() {
  omp_set_dynamic(0);
  check_parameter_validation();
  check_horizon_root_algorithms();
  check_geometry_case(make_param(0, 1, 48, "-500", "1000", "1e-6"),
                      "small Lambda");
  check_geometry_case(make_param(2, 2, 48, "-200", "400", "0.44"),
                      "near Nariai");
  check_neighbor_reuse_consistency();
  check_schwarzschild_limit();
#ifndef SDS_SKIP_2000_DIGIT_REFERENCE
    check_2000_digit_reference();
#endif
  check_multipole_constructor();
  check_source_profiles();
  check_rhs_and_sources();
  check_dopri5_trajectory();
  check_randomized_legacy_profile_evolution();

  if(failures != 0) {
    std::cerr << failures << " precise SdS correctness checks failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: precise SdS correctness matrix; largest relative error="
            << largest_relative_error
            << "; largest analytic/root-found horizon difference="
            << largest_horizon_formula_difference
            << "; largest randomized legacy-profile difference="
            << largest_legacy_profile_difference << '\n';
  return EXIT_SUCCESS;
}

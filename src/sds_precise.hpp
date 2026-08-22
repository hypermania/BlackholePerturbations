/*!
  \file sds_precise.hpp
  \author Siyang Ling
  \brief High-precision sourced master equations in Schwarzschild-de Sitter.

  The evolved equation is

    partial_t psi = Pi,
    partial_t Pi  = partial_x^2 psi - V_{s l}(r) psi + S(t, x),

  where x is the SdS tortoise coordinate and

    V_{s l} = f(r) [l(l+1)/r^2 + (1-s^2) 2M/r^3].

  The evolution uses binary128 arithmetic. Geometry preprocessing uses
  100-decimal-digit arithmetic and a horizon-distance parameterization so
  exponentially small values of f are never obtained by subtracting rounded
  quantities.
*/
#ifndef SDS_PRECISE_HPP
#define SDS_PRECISE_HPP

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <boost/multiprecision/float128.hpp>


enum class SdSSourceProfile : long long {
  ArealPower = 0,
  HorizonSubtractedArealPower = 1,
  LocalScalar = 2,
  TortoisePower = 3
};

enum class SdSWaveform : long long {
  Gaussian = 0,
  GaussianDerivative = 1
};


inline const char *sds_source_profile_name(const SdSSourceProfile profile) {
  switch(profile) {
    case SdSSourceProfile::ArealPower: return "areal_power";
    case SdSSourceProfile::HorizonSubtractedArealPower:
      return "horizon_subtracted_areal_power";
    case SdSSourceProfile::LocalScalar: return "local_scalar";
    case SdSSourceProfile::TortoisePower: return "tortoise_power";
  }
  return "unknown";
}


inline const char *sds_waveform_name(const SdSWaveform waveform) {
  switch(waveform) {
    case SdSWaveform::Gaussian: return "gaussian";
    case SdSWaveform::GaussianDerivative: return "gaussian_derivative";
  }
  return "unknown";
}


struct SdSMasterPDEPreciseParam {
  typedef boost::multiprecision::float128 Scalar;
  long long int s;
  long long int l;
  Scalar M;
  Scalar Lambda;
  Scalar r_min;
  Scalar r_max;
  long long int N;
  Scalar t_start;
  Scalar t_end;
  Scalar t_interval;
  Scalar delta_t;
};


struct SdSTranslatedSourceParam {
  typedef boost::multiprecision::float128 Scalar;
  SdSSourceProfile profile = SdSSourceProfile::ArealPower;
  SdSWaveform waveform = SdSWaveform::Gaussian;
  Scalar beta = Scalar(2);
  Scalar amplitude = Scalar(1);
  Scalar u_center = Scalar(-10);
  Scalar sigma = Scalar("0.5");
  Scalar onset_time = Scalar(0);
  Scalar cutoff_sigma = Scalar(12);
  Scalar L = Scalar(1);
  Scalar x0 = Scalar(100);
  Scalar X0 = Scalar(0);
  Scalar X1 = Scalar(20);
};


struct SdSMasterPDEPrecise {
  typedef SdSMasterPDEPreciseParam Param;
  typedef boost::multiprecision::float128 Scalar;
  typedef boost::multiprecision::cpp_bin_float_100 HighPrecisionScalar;
  typedef Eigen::Array<Scalar, -1, 1> Vector;
  typedef Vector State;

  Param param;
  long long int grid_size;
  Scalar inv_h;
  Scalar inv_h_sqr;

  Scalar r_black_hole;
  Scalar r_cosmological;
  Scalar r_negative;
  Scalar kappa_black_hole;
  Scalar kappa_cosmological;
  Scalar max_inversion_residual;

  Vector r_ast;
  Vector r;
  Vector rho_black_hole;
  Vector rho_cosmological;
  Vector f;
  Vector V;
  Vector center_factor;

  std::function<void(const Scalar &, Vector &)> Q;
  mutable Vector Q_workspace;

  bool has_translated_source = false;
  SdSTranslatedSourceParam translated_source_param;
  Vector translated_source_spatial;
  Scalar waveform_prefactor = 0;

  explicit SdSMasterPDEPrecise(const Param param_) : param(param_) {
    validate_param();

    grid_size = param.N + 1;
    const Scalar h = (param.r_max - param.r_min) / Scalar(param.N - 1);
    inv_h = Scalar(1) / h;
    inv_h_sqr = inv_h * inv_h;

    r_ast.resize(grid_size);
    r.resize(grid_size);
    rho_black_hole.resize(grid_size);
    rho_cosmological.resize(grid_size);
    f.resize(grid_size);
    V.resize(grid_size);
    center_factor.resize(grid_size);
    Q_workspace.resize(grid_size);

    const GeometryHP geometry = make_geometry(
        HighPrecisionScalar(param.M), HighPrecisionScalar(param.Lambda));
    r_black_hole = geometry.rb.convert_to<Scalar>();
    r_cosmological = geometry.rc.convert_to<Scalar>();
    r_negative = geometry.ro.convert_to<Scalar>();
    kappa_black_hole = geometry.kappa_b.convert_to<Scalar>();
    kappa_cosmological = geometry.kappa_c.convert_to<Scalar>();

    const HighPrecisionScalar x_min_hp(param.r_min);
    const HighPrecisionScalar h_hp =
        (HighPrecisionScalar(param.r_max) - x_min_hp)
        / HighPrecisionScalar(param.N - 1);
    std::vector<Scalar> residuals(static_cast<std::size_t>(grid_size));
    std::atomic<bool> inversion_failed(false);

#pragma omp parallel for schedule(static)
    for(long long int i = 0; i < grid_size; ++i) {
      const HighPrecisionScalar x_hp =
          x_min_hp + HighPrecisionScalar(i) * h_hp - h_hp / 2;
      const InversionHP point = invert_tortoise(geometry, x_hp);
      if(!point.success) inversion_failed.store(true, std::memory_order_relaxed);

      const HighPrecisionScalar f_hp =
          geometry.lambda * point.rho_b * point.rho_c * point.r_minus_ro
          / (HighPrecisionScalar(3) * point.r);
      const HighPrecisionScalar inv_r = HighPrecisionScalar(1) / point.r;
      const HighPrecisionScalar potential_hp = f_hp * (
          HighPrecisionScalar(param.l * (param.l + 1)) * inv_r * inv_r
          + HighPrecisionScalar(1 - param.s * param.s)
                * HighPrecisionScalar(2) * geometry.mass
                * inv_r * inv_r * inv_r);

      r_ast[i] = x_hp.convert_to<Scalar>();
      r[i] = point.r.convert_to<Scalar>();
      rho_black_hole[i] = point.rho_b.convert_to<Scalar>();
      rho_cosmological[i] = point.rho_c.convert_to<Scalar>();
      f[i] = f_hp.convert_to<Scalar>();
      V[i] = potential_hp.convert_to<Scalar>();
      residuals[static_cast<std::size_t>(i)] = point.residual.convert_to<Scalar>();
    }

    if(inversion_failed.load(std::memory_order_relaxed)) {
      throw std::runtime_error("SdS tortoise-coordinate inversion failed");
    }
    max_inversion_residual = 0;
    for(const Scalar &residual : residuals) {
      max_inversion_residual = std::max(max_inversion_residual, residual);
    }

    const Scalar d2_factor = inv_h_sqr / Scalar(12);
    center_factor = -Scalar(30) * d2_factor - V;
  }


  static Vector compute_r_ast_vector(const Scalar r_min, const Scalar r_max,
                                     const long long int N) {
    if(N < 4 || !(r_max > r_min)) {
      throw std::invalid_argument("invalid SdS grid");
    }
    const Scalar h = (r_max - r_min) / Scalar(N - 1);
    Vector result(N + 1);
#pragma omp parallel for schedule(static)
    for(long long int i = 0; i <= N; ++i) {
      result[i] = r_min + Scalar(i) * h - h / Scalar(2);
    }
    return result;
  }


  void set_generic_source(
      std::function<void(const Scalar &, Vector &)> source) {
    Q = std::move(source);
    has_translated_source = false;
    translated_source_spatial.resize(0);
  }


  void clear_source() {
    Q = {};
    has_translated_source = false;
    translated_source_spatial.resize(0);
  }


  void set_translated_gaussian_source(
      const SdSTranslatedSourceParam &source_param) {
    validate_source_param(source_param);
    translated_source_param = source_param;
    translated_source_spatial.resize(grid_size);

    const Scalar beta = source_param.beta;
    const Scalar rc_power = pow(r_cosmological, -beta);

#pragma omp parallel for schedule(static)
    for(long long int i = 0; i < grid_size; ++i) {
      Scalar profile = 0;
      switch(source_param.profile) {
        case SdSSourceProfile::ArealPower:
          profile = pow(r[i], -beta);
          break;
        case SdSSourceProfile::HorizonSubtractedArealPower: {
          const Scalar log_ratio = stable_log1p(-rho_cosmological[i]
                                                / r_cosmological);
          profile = rc_power * stable_expm1(-beta * log_ratio);
          break;
        }
        case SdSSourceProfile::LocalScalar:
          profile = f[i] * pow(r[i], -beta);
          break;
        case SdSSourceProfile::TortoisePower: {
          Scalar cutoff = 0;
          if(r_ast[i] >= source_param.X1) {
            cutoff = 1;
          } else if(r_ast[i] > source_param.X0) {
            const Scalar z = (r_ast[i] - source_param.X0)
                             / (source_param.X1 - source_param.X0);
            const Scalar left = exp(-Scalar(1) / z);
            const Scalar right = exp(-Scalar(1) / (Scalar(1) - z));
            cutoff = left / (left + right);
          }
          profile = cutoff == 0
              ? Scalar(0)
              : cutoff
                    * pow(source_param.L / (r_ast[i] + source_param.x0), beta);
          break;
        }
      }
      translated_source_spatial[i] = profile;
    }

    const Scalar pi = boost::math::constants::pi<Scalar>();
    waveform_prefactor = source_param.amplitude
                         / (sqrt(Scalar(2) * pi) * source_param.sigma);
    Q = {};
    has_translated_source = true;
  }


  Scalar translated_source_value(const long long int i,
                                 const Scalar &t) const {
    if(!has_translated_source || t < translated_source_param.onset_time) {
      return 0;
    }
    const Scalar z = (t - r_ast[i] - translated_source_param.u_center)
                     / translated_source_param.sigma;
    if(abs(z) > translated_source_param.cutoff_sigma) return 0;
    Scalar waveform = waveform_prefactor * exp(-z * z / Scalar(2));
    if(translated_source_param.waveform
       == SdSWaveform::GaussianDerivative) {
      waveform *= -z / translated_source_param.sigma;
    }
    return translated_source_spatial[i] * waveform;
  }


  void operator()(const State &state, State &derivative, const Scalar t) {
    if(state.size() != 2 * grid_size || derivative.size() != 2 * grid_size) {
      throw std::invalid_argument("SdS state size does not match the grid");
    }

    const Scalar *__restrict__ psi = state.data();
    const Scalar *__restrict__ pi = state.data() + grid_size;
    Scalar *__restrict__ dpsi = derivative.data();
    Scalar *__restrict__ dpi = derivative.data() + grid_size;
    const Scalar *__restrict__ potential = V.data();
    const Scalar *__restrict__ center = center_factor.data();

    const bool has_generic_source = !has_translated_source
                                    && static_cast<bool>(Q);
    if(has_generic_source) {
      Q(t, Q_workspace);
      if(Q_workspace.size() != grid_size) {
        throw std::invalid_argument("generic SdS source size does not match grid");
      }
    }
    const Scalar *__restrict__ generic_source = Q_workspace.data();

    long long int source_begin = 1;
    long long int source_end = 0;
    if(has_translated_source && t >= translated_source_param.onset_time) {
      const Scalar radius = translated_source_param.cutoff_sigma
                            * translated_source_param.sigma;
      const Scalar center_x = t - translated_source_param.u_center;
      const Scalar first_x = param.r_min
                             - Scalar(1) / (Scalar(2) * inv_h);
      const Scalar begin_real = (center_x - radius - first_x) * inv_h;
      const Scalar end_real = (center_x + radius - first_x) * inv_h;
      if(!(end_real < 0 || begin_real > Scalar(grid_size - 1))) {
        const Scalar clipped_begin = std::max(Scalar(0), begin_real);
        const Scalar clipped_end = std::min(Scalar(grid_size - 1), end_real);
        source_begin = ceil(clipped_begin).convert_to<long long int>();
        source_end = floor(clipped_end).convert_to<long long int>();
      }
    }

    const bool translated_active = source_begin <= source_end;
    const Scalar d2_factor = inv_h_sqr / Scalar(12);
    const Scalar near_factor = Scalar(16) * d2_factor;
    const Scalar far_factor = -d2_factor;
    const Scalar one_twelfth_inv_h = inv_h / Scalar(12);

    auto source_at = [&](const long long int i) -> Scalar {
      if(translated_active && i >= source_begin && i <= source_end) {
        return translated_source_value(i, t);
      }
      if(has_generic_source) return generic_source[i];
      return Scalar(0);
    };

#pragma omp parallel
    {
      if(!translated_active && !has_generic_source) {
#pragma omp for schedule(static)
        for(long long int i = 2; i <= grid_size - 3; ++i) {
          const Scalar near_sum = psi[i - 1] + psi[i + 1];
          const Scalar far_sum = psi[i - 2] + psi[i + 2];
          dpi[i] = near_factor * near_sum + far_factor * far_sum
                   + center[i] * psi[i];
          dpsi[i] = pi[i];
        }
      } else if(translated_active) {
#pragma omp for schedule(static)
        for(long long int i = 2; i <= grid_size - 3; ++i) {
          const Scalar near_sum = psi[i - 1] + psi[i + 1];
          const Scalar far_sum = psi[i - 2] + psi[i + 2];
          dpi[i] = near_factor * near_sum + far_factor * far_sum
                   + center[i] * psi[i];
          if(i >= source_begin && i <= source_end) {
            dpi[i] += translated_source_value(i, t);
          }
          dpsi[i] = pi[i];
        }
      } else {
#pragma omp for schedule(static)
        for(long long int i = 2; i <= grid_size - 3; ++i) {
          const Scalar near_sum = psi[i - 1] + psi[i + 1];
          const Scalar far_sum = psi[i - 2] + psi[i + 2];
          dpi[i] = near_factor * near_sum + far_factor * far_sum
                   + center[i] * psi[i] + generic_source[i];
          dpsi[i] = pi[i];
        }
      }

#pragma omp single nowait
      {
        dpi[0] = (-Scalar(25) * pi[0] + Scalar(48) * pi[1]
                  - Scalar(36) * pi[2] + Scalar(16) * pi[3]
                  - Scalar(3) * pi[4]) * one_twelfth_inv_h
                 - potential[0] * psi[0] + source_at(0);
        dpsi[0] = pi[0];

        dpi[1] = (Scalar(11) * psi[0] - Scalar(20) * psi[1]
                  + Scalar(6) * psi[2] + Scalar(4) * psi[3] - psi[4])
                 * d2_factor - potential[1] * psi[1] + source_at(1);
        dpsi[1] = pi[1];

        const long long int n2 = grid_size - 2;
        const long long int n1 = grid_size - 1;
        dpi[n2] = (-psi[grid_size - 5] + Scalar(4) * psi[grid_size - 4]
                   + Scalar(6) * psi[grid_size - 3] - Scalar(20) * psi[n2]
                   + Scalar(11) * psi[n1]) * d2_factor
                  - potential[n2] * psi[n2] + source_at(n2);
        dpsi[n2] = pi[n2];

        dpi[n1] = (-Scalar(3) * pi[grid_size - 5]
                   + Scalar(16) * pi[grid_size - 4]
                   - Scalar(36) * pi[grid_size - 3]
                   + Scalar(48) * pi[n2] - Scalar(25) * pi[n1])
                  * one_twelfth_inv_h - potential[n1] * psi[n1]
                  + source_at(n1);
        dpsi[n1] = pi[n1];
      }
    }
  }

 private:
  struct GeometryHP {
    HighPrecisionScalar mass;
    HighPrecisionScalar lambda;
    HighPrecisionScalar rb;
    HighPrecisionScalar rc;
    HighPrecisionScalar ro;
    HighPrecisionScalar delta;
    HighPrecisionScalar kappa_b;
    HighPrecisionScalar kappa_c;
    HighPrecisionScalar inv_fp_b;
    HighPrecisionScalar inv_fp_c;
    HighPrecisionScalar inv_fp_o;
    HighPrecisionScalar reference_r;
    HighPrecisionScalar reference_x;
    HighPrecisionScalar reference_rb;
    HighPrecisionScalar reference_rc;
    HighPrecisionScalar reference_ro;
    HighPrecisionScalar midpoint_x;
  };

  struct EvaluationHP {
    HighPrecisionScalar x;
    HighPrecisionScalar derivative;
    HighPrecisionScalar second_derivative;
    HighPrecisionScalar r;
    HighPrecisionScalar rho_b;
    HighPrecisionScalar rho_c;
    HighPrecisionScalar r_minus_ro;
  };

  struct InversionHP {
    HighPrecisionScalar r;
    HighPrecisionScalar rho_b;
    HighPrecisionScalar rho_c;
    HighPrecisionScalar r_minus_ro;
    HighPrecisionScalar residual;
    bool success;
  };


  void validate_param() const {
    if(!(param.M > 0)) throw std::invalid_argument("SdS requires M > 0");
    if(!(param.Lambda > 0)) {
      throw std::invalid_argument("SdS requires Lambda > 0");
    }
    if(!(Scalar(9) * param.Lambda * param.M * param.M < Scalar(1))) {
      throw std::invalid_argument("SdS requires 9 Lambda M^2 < 1");
    }
    if(param.s < 0 || param.s > 2) {
      throw std::invalid_argument("SdS master spin must be 0, 1, or 2");
    }
    if(param.l < param.s) {
      throw std::invalid_argument("SdS multipole must satisfy l >= s");
    }
    if(param.N < 4) throw std::invalid_argument("SdS requires N >= 4");
    if(!(param.r_max > param.r_min)) {
      throw std::invalid_argument("SdS requires r_max > r_min");
    }
    if(!(param.delta_t > 0) || !(param.t_end >= param.t_start)
       || !(param.t_interval > 0)) {
      throw std::invalid_argument("invalid SdS time parameters");
    }
  }


  static void validate_source_param(
      const SdSTranslatedSourceParam &source_param) {
    if(!(source_param.beta > 0)) {
      throw std::invalid_argument("SdS source requires beta > 0");
    }
    if(!(source_param.sigma > 0) || !(source_param.cutoff_sigma > 0)) {
      throw std::invalid_argument("SdS source requires positive Gaussian widths");
    }
    if(source_param.profile == SdSSourceProfile::TortoisePower) {
      if(!(source_param.L > 0) || !(source_param.X1 > source_param.X0)
         || !(source_param.X0 + source_param.x0 > 0)) {
        throw std::invalid_argument("invalid algebraic tortoise source profile");
      }
    }
  }


  static Scalar stable_log1p(const Scalar &value) {
    if(abs(value) >= Scalar("0.01")) return log(Scalar(1) + value);
    Scalar sum = 0;
    Scalar power = value;
    for(int order = 1; order <= 160; ++order) {
      const Scalar term = power / Scalar(order);
      sum += (order % 2 == 1) ? term : -term;
      if(abs(term) < Scalar("1e-38")) break;
      power *= value;
    }
    return sum;
  }


  static Scalar stable_expm1(const Scalar &value) {
    if(abs(value) >= Scalar("0.01")) return exp(value) - Scalar(1);
    Scalar sum = value;
    Scalar term = value;
    for(int order = 2; order <= 160; ++order) {
      term *= value / Scalar(order);
      sum += term;
      if(abs(term) < Scalar("1e-38")) break;
    }
    return sum;
  }


  static GeometryHP make_geometry(const HighPrecisionScalar &mass,
                                  const HighPrecisionScalar &lambda) {
    const HighPrecisionScalar pi =
        boost::math::constants::pi<HighPrecisionScalar>();
    const HighPrecisionScalar sqrt_lambda = sqrt(lambda);
    const HighPrecisionScalar angle = acos(HighPrecisionScalar(3) * mass
                                            * sqrt_lambda)
                                      / HighPrecisionScalar(3);

    GeometryHP geometry;
    geometry.mass = mass;
    geometry.lambda = lambda;
    geometry.rb = HighPrecisionScalar(2) / sqrt_lambda
                  * cos(angle + pi / HighPrecisionScalar(3));
    geometry.rc = HighPrecisionScalar(2) / sqrt_lambda
                  * cos(angle - pi / HighPrecisionScalar(3));
    geometry.ro = -(geometry.rb + geometry.rc);
    geometry.delta = geometry.rc - geometry.rb;

    auto f_prime = [&](const HighPrecisionScalar &radius) {
      return HighPrecisionScalar(2) * mass / (radius * radius)
             - HighPrecisionScalar(2) * lambda * radius
                   / HighPrecisionScalar(3);
    };
    const HighPrecisionScalar fp_b = f_prime(geometry.rb);
    const HighPrecisionScalar fp_c = f_prime(geometry.rc);
    const HighPrecisionScalar fp_o = f_prime(geometry.ro);
    geometry.kappa_b = fp_b / HighPrecisionScalar(2);
    geometry.kappa_c = -fp_c / HighPrecisionScalar(2);
    geometry.inv_fp_b = HighPrecisionScalar(1) / fp_b;
    geometry.inv_fp_c = HighPrecisionScalar(1) / fp_c;
    geometry.inv_fp_o = HighPrecisionScalar(1) / fp_o;

    geometry.reference_r = HighPrecisionScalar(3) * mass;
    geometry.reference_x = geometry.reference_r
                           + HighPrecisionScalar(2) * mass
                                 * log(HighPrecisionScalar("0.5"));
    geometry.reference_rb = geometry.reference_r - geometry.rb;
    geometry.reference_rc = geometry.rc - geometry.reference_r;
    geometry.reference_ro = geometry.reference_r - geometry.ro;
    geometry.midpoint_x = evaluate_y(geometry, HighPrecisionScalar(0)).x;
    return geometry;
  }


  static EvaluationHP evaluate_y(const GeometryHP &geometry,
                                 const HighPrecisionScalar &y) {
    EvaluationHP result;
    if(y >= 0) {
      const HighPrecisionScalar exponential = exp(-y);
      const HighPrecisionScalar denominator = HighPrecisionScalar(1)
                                               + exponential;
      result.rho_b = geometry.delta / denominator;
      result.rho_c = geometry.delta * exponential / denominator;
    } else {
      const HighPrecisionScalar exponential = exp(y);
      const HighPrecisionScalar denominator = HighPrecisionScalar(1)
                                               + exponential;
      result.rho_b = geometry.delta * exponential / denominator;
      result.rho_c = geometry.delta / denominator;
    }
    result.r = geometry.rb + result.rho_b;
    result.r_minus_ro = geometry.rb - geometry.ro + result.rho_b;

    result.x = geometry.reference_x
               + geometry.inv_fp_b
                     * (log(result.rho_b) - log(geometry.reference_rb))
               + geometry.inv_fp_c
                     * (log(result.rho_c) - log(geometry.reference_rc))
               + geometry.inv_fp_o
                     * (log(result.r_minus_ro) - log(geometry.reference_ro));

    const HighPrecisionScalar prefactor =
        HighPrecisionScalar(3) / (geometry.lambda * geometry.delta);
    result.derivative = prefactor * result.r / result.r_minus_ro;
    const HighPrecisionScalar dr_dy = result.rho_b * result.rho_c
                                      / geometry.delta;
    result.second_derivative = prefactor * (-geometry.ro) * dr_dy
                               / (result.r_minus_ro * result.r_minus_ro);
    return result;
  }


  static InversionHP invert_tortoise(const GeometryHP &geometry,
                                     const HighPrecisionScalar &target_x) {
    const HighPrecisionScalar tolerance("1e-80");
    const HighPrecisionScalar target_scale =
        HighPrecisionScalar(1) + abs(target_x);
    if(abs(target_x - geometry.midpoint_x) <= tolerance * target_scale) {
      const EvaluationHP midpoint = evaluate_y(geometry, HighPrecisionScalar(0));
      return {midpoint.r, midpoint.rho_b, midpoint.rho_c,
              midpoint.r_minus_ro, abs(midpoint.x - target_x), true};
    }

    HighPrecisionScalar lower;
    HighPrecisionScalar upper;
    HighPrecisionScalar guess;
    if(target_x < geometry.midpoint_x) {
      upper = 0;
      guess = HighPrecisionScalar(2) * geometry.kappa_b
              * (target_x - geometry.midpoint_x);
      lower = std::min(guess, HighPrecisionScalar(-1));
      for(int expansion = 0;
          evaluate_y(geometry, lower).x > target_x && expansion < 100;
          ++expansion) {
        lower = HighPrecisionScalar(2) * lower - HighPrecisionScalar(1);
      }
    } else {
      lower = 0;
      guess = HighPrecisionScalar(2) * geometry.kappa_c
              * (target_x - geometry.midpoint_x);
      upper = std::max(guess, HighPrecisionScalar(1));
      for(int expansion = 0;
          evaluate_y(geometry, upper).x < target_x && expansion < 100;
          ++expansion) {
        upper = HighPrecisionScalar(2) * upper + HighPrecisionScalar(1);
      }
    }

    if(evaluate_y(geometry, lower).x > target_x
       || evaluate_y(geometry, upper).x < target_x) {
      return {0, 0, 0, 0, std::numeric_limits<HighPrecisionScalar>::infinity(),
              false};
    }

    HighPrecisionScalar y = std::max(lower, std::min(upper, guess));
    if(y == lower || y == upper) y = (lower + upper) / 2;
    EvaluationHP evaluation = evaluate_y(geometry, y);
    bool converged = false;

    for(int iteration = 0; iteration < 120; ++iteration) {
      const HighPrecisionScalar residual = evaluation.x - target_x;
      if(abs(residual) <= tolerance * target_scale) {
        converged = true;
        break;
      }
      if(residual < 0) lower = y;
      else upper = y;

      const HighPrecisionScalar denominator =
          HighPrecisionScalar(2) * evaluation.derivative
              * evaluation.derivative
          - residual * evaluation.second_derivative;
      HighPrecisionScalar candidate = (lower + upper) / 2;
      if(denominator != 0) {
        const HighPrecisionScalar halley = y
            - HighPrecisionScalar(2) * residual * evaluation.derivative
                  / denominator;
        if(halley > lower && halley < upper) candidate = halley;
      }
      y = candidate;
      evaluation = evaluate_y(geometry, y);
    }

    const HighPrecisionScalar residual = abs(evaluation.x - target_x);
    converged = converged || residual <= tolerance * target_scale;
    return {evaluation.r, evaluation.rho_b, evaluation.rho_c,
            evaluation.r_minus_ro, residual, converged};
  }
};

#endif

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
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <boost/multiprecision/float128.hpp>
#include <omp.h>
#include <quadmath.h>


template<typename Number>
inline Number sds_grid_coordinate(const Number &x_min, const Number &h,
                                  const long long int i) {
  return x_min + (Number(i) - Number(1) / Number(2)) * h;
}


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


struct SdSSpacetimeGaussianSourceParam {
  typedef boost::multiprecision::float128 Scalar;
  Scalar amplitude = Scalar(1);
  Scalar r_ast_center = Scalar(0);
  Scalar t_center = Scalar(0);
  Scalar sigma = Scalar("0.5");
  Scalar cutoff_sigma = Scalar(12);
};


/*!
  \brief Optional source for SdSMasterPDEPrecise.

  Built-in sources precompute their spatial data during PDE construction and
  fill a caller-owned workspace at each time. The spacetime Gaussian is

    amplitude exp(-[(x-r_ast_center)^2+(t-t_center)^2]/sigma^2)
      / (sqrt(2 pi) sigma).

  A GenericSource obeying the same workspace-filling contract can be used for
  arbitrary sources.
*/
struct SdSSource {
  typedef boost::multiprecision::float128 Scalar;
  typedef Eigen::Array<Scalar, -1, 1> Vector;
  typedef std::function<void(const Scalar &, Vector &)> GenericSource;

  enum class Kind : long long {
    None = 0,
    Generic = 1,
    TranslatedGaussian = 2,
    SpacetimeGaussian = 3
  };

  SdSSource() = default;

  explicit SdSSource(GenericSource source)
      : kind(Kind::Generic), generic_source(std::move(source)) {
    if(!generic_source) {
      throw std::invalid_argument("generic SdS source is empty");
    }
  }

  explicit SdSSource(const SdSTranslatedSourceParam &source_param)
      : kind(Kind::TranslatedGaussian), translated_param(source_param) {
    validate(translated_param);
  }

  explicit SdSSource(const SdSSpacetimeGaussianSourceParam &source_param)
      : kind(Kind::SpacetimeGaussian), spacetime_param(source_param) {
    validate(spacetime_param);
  }

  explicit operator bool() const {
    return kind != Kind::None;
  }

  void initialize(const Scalar &x_min_, const Scalar &h_,
                  const long long int grid_size_, const Scalar &r_cosmological,
                  const Vector &r, const Vector &rho_cosmological,
                  const Vector &f) {
    x_min = x_min_;
    h = h_;
    inv_h = Scalar(1) / h;
    grid_size = grid_size_;

    if(kind == Kind::TranslatedGaussian) {
      initialize_translated(r_cosmological, r, rho_cosmological, f);
    } else if(kind == Kind::SpacetimeGaussian) {
      initialize_spacetime_gaussian();
    }
  }

  void operator()(const Scalar &t, Vector &result) const {
    switch(kind) {
      case Kind::None:
        result.setZero(grid_size);
        return;
      case Kind::Generic:
        generic_source(t, result);
        return;
      case Kind::TranslatedGaussian:
        evaluate_translated(t, result);
        return;
      case Kind::SpacetimeGaussian:
        evaluate_spacetime_gaussian(t, result);
        return;
    }
  }

 private:
  Kind kind = Kind::None;
  GenericSource generic_source;
  SdSTranslatedSourceParam translated_param;
  SdSSpacetimeGaussianSourceParam spacetime_param;

  long long int grid_size = 0;
  Scalar x_min = 0;
  Scalar h = 1;
  Scalar inv_h = 1;
  Scalar normalization = 0;
  long long int spatial_begin = 1;
  long long int spatial_end = 0;
  Vector spatial_profile;

  static void validate(const SdSTranslatedSourceParam &source_param) {
    if(!boost::math::isfinite(source_param.beta)) {
      throw std::invalid_argument("SdS source beta must be finite");
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

  static void validate(const SdSSpacetimeGaussianSourceParam &source_param) {
    if(!(source_param.sigma > 0) || !(source_param.cutoff_sigma > 0)) {
      throw std::invalid_argument(
          "SdS spacetime Gaussian requires positive widths");
    }
  }

  std::pair<long long int, long long int> grid_window(
      const Scalar &center, const Scalar &radius) const {
    const Scalar first_x = sds_grid_coordinate(x_min, h, 0);
    const Scalar begin_real = (center - radius - first_x) * inv_h;
    const Scalar end_real = (center + radius - first_x) * inv_h;
    if(end_real < 0 || begin_real > Scalar(grid_size - 1)) return {1, 0};
    const long long int begin = ceil(std::max(Scalar(0), begin_real))
                                    .convert_to<long long int>();
    const long long int end = floor(std::min(
        Scalar(grid_size - 1), end_real)).convert_to<long long int>();
    return {begin, end};
  }

  void initialize_translated(const Scalar &r_cosmological,
                             const Vector &r,
                             const Vector &rho_cosmological,
                             const Vector &f) {
    spatial_profile.resize(grid_size);
    const Scalar beta = translated_param.beta;
    const Scalar rc_power = pow(r_cosmological, -beta);

#pragma omp parallel for schedule(static)
    for(long long int i = 0; i < grid_size; ++i) {
      const Scalar x = sds_grid_coordinate(x_min, h, i);
      Scalar profile = 0;
      switch(translated_param.profile) {
        case SdSSourceProfile::ArealPower:
          profile = pow(r[i], -beta);
          break;
        case SdSSourceProfile::HorizonSubtractedArealPower: {
          const Scalar ratio = -rho_cosmological[i] / r_cosmological;
          // The vendored Boost wrappers call these libquadmath functions but
          // rely on an implicit __float128 conversion rejected by GCC 15.
          const Scalar log_ratio(
              Scalar::backend_type(log1pq(ratio.backend().value())));
          const Scalar exponent = -beta * log_ratio;
          const Scalar difference(
              Scalar::backend_type(expm1q(exponent.backend().value())));
          profile = rc_power * difference;
          break;
        }
        case SdSSourceProfile::LocalScalar:
          profile = f[i] * pow(r[i], -beta);
          break;
        case SdSSourceProfile::TortoisePower: {
          Scalar cutoff = 0;
          if(x >= translated_param.X1) {
            cutoff = 1;
          } else if(x > translated_param.X0) {
            const Scalar z = (x - translated_param.X0)
                             / (translated_param.X1 - translated_param.X0);
            const Scalar left = exp(-Scalar(1) / z);
            const Scalar right = exp(-Scalar(1) / (Scalar(1) - z));
            cutoff = left / (left + right);
          }
          profile = cutoff == 0
              ? Scalar(0)
              : cutoff * pow(translated_param.L
                                 / (x + translated_param.x0), beta);
          break;
        }
      }
      spatial_profile[i] = profile;
    }

    const Scalar pi = boost::math::constants::pi<Scalar>();
    normalization = translated_param.amplitude
                    / (sqrt(Scalar(2) * pi) * translated_param.sigma);
  }

  void initialize_spacetime_gaussian() {
    spatial_profile = Vector::Zero(grid_size);
    const Scalar radius = spacetime_param.cutoff_sigma * spacetime_param.sigma;
    std::tie(spatial_begin, spatial_end) = grid_window(
        spacetime_param.r_ast_center, radius);
    const Scalar inv_sigma_sqr = Scalar(1)
        / (spacetime_param.sigma * spacetime_param.sigma);

#pragma omp parallel for schedule(static)
    for(long long int i = spatial_begin; i <= spatial_end; ++i) {
      const Scalar dx = sds_grid_coordinate(x_min, h, i)
                        - spacetime_param.r_ast_center;
      spatial_profile[i] = exp(-dx * dx * inv_sigma_sqr);
    }

    const Scalar pi = boost::math::constants::pi<Scalar>();
    normalization = spacetime_param.amplitude
                    / (sqrt(Scalar(2) * pi) * spacetime_param.sigma);
  }

  void evaluate_translated(const Scalar &t, Vector &result) const {
    result.setZero(grid_size);
    if(t < translated_param.onset_time) return;

    const Scalar radius = translated_param.cutoff_sigma
                          * translated_param.sigma;
    const Scalar center_x = t - translated_param.u_center;
    const auto [begin, end] = grid_window(center_x, radius);

#pragma omp parallel for schedule(static)
    for(long long int i = begin; i <= end; ++i) {
      const Scalar x = sds_grid_coordinate(x_min, h, i);
      const Scalar z = (t - x - translated_param.u_center)
                       / translated_param.sigma;
      Scalar waveform = normalization * exp(-z * z / Scalar(2));
      if(translated_param.waveform == SdSWaveform::GaussianDerivative) {
        waveform *= -z / translated_param.sigma;
      }
      result[i] = spatial_profile[i] * waveform;
    }
  }

  void evaluate_spacetime_gaussian(const Scalar &t, Vector &result) const {
    result.setZero(grid_size);
    const Scalar dt = t - spacetime_param.t_center;
    const Scalar z_t = dt / spacetime_param.sigma;
    if(abs(z_t) > spacetime_param.cutoff_sigma) return;
    const Scalar time_factor = normalization * exp(-z_t * z_t);

#pragma omp parallel for schedule(static)
    for(long long int i = spatial_begin; i <= spatial_end; ++i) {
      result[i] = time_factor * spatial_profile[i];
    }
  }
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

  Vector r;
  Vector rho_cosmological;
  Vector f;
  Vector V;

  SdSSource Q;
  Vector Q_workspace;

  explicit SdSMasterPDEPrecise(const Param param_, SdSSource source = {})
      : param(param_), Q(std::move(source)) {
    validate_param();

    grid_size = param.N + 1;
    const Scalar h = grid_space();
    inv_h = Scalar(1) / h;
    inv_h_sqr = inv_h * inv_h;

    r.resize(grid_size);
    rho_cosmological.resize(grid_size);
    f.resize(grid_size);
    V.resize(grid_size);
    Q_workspace = Vector::Zero(grid_size);

    const GeometryHP geometry = make_geometry(
        HighPrecisionScalar(param.M), HighPrecisionScalar(param.Lambda));
    r_black_hole = geometry.rb.convert_to<Scalar>();
    r_cosmological = geometry.rc.convert_to<Scalar>();
    r_negative = geometry.ro.convert_to<Scalar>();
    const HighPrecisionScalar mass_hp(param.M);
    const HighPrecisionScalar lambda_hp(param.Lambda);
    const HighPrecisionScalar spin_hp(param.s);
    const HighPrecisionScalar multipole_hp(param.l);
    const HighPrecisionScalar angular_coefficient =
        multipole_hp * (multipole_hp + HighPrecisionScalar(1));
    const HighPrecisionScalar spin_coefficient =
        HighPrecisionScalar(1) - spin_hp * spin_hp;
    auto f_prime = [&](const HighPrecisionScalar &radius) {
      return HighPrecisionScalar(2) * mass_hp / (radius * radius)
             - HighPrecisionScalar(2) * lambda_hp * radius
                   / HighPrecisionScalar(3);
    };
    kappa_black_hole = (f_prime(geometry.rb) / HighPrecisionScalar(2))
                               .convert_to<Scalar>();
    kappa_cosmological = (-f_prime(geometry.rc) / HighPrecisionScalar(2))
                                .convert_to<Scalar>();

    const HighPrecisionScalar x_min_hp(param.r_min);
    const HighPrecisionScalar h_hp(h);
    long long int first_failed_index = grid_size;

#pragma omp parallel reduction(min:first_failed_index)
    {
      const long long int thread = omp_get_thread_num();
      const long long int thread_count = omp_get_num_threads();
      const long long int begin = thread * grid_size / thread_count;
      const long long int end = (thread + 1) * grid_size / thread_count;
      TortoisePointHP previous;
      bool have_previous = false;

      for(long long int i = begin; i < end; ++i) {
        const HighPrecisionScalar x_hp = grid_coordinate(
            x_min_hp, h_hp, i);
        const TortoisePointHP point = invert_tortoise(
            geometry, x_hp, have_previous ? &previous : nullptr);
        if(!point.success) {
          first_failed_index = std::min(first_failed_index, i);
        }

        const HighPrecisionScalar f_hp =
            geometry.lambda * point.rho_b * point.rho_c
            * (point.r - geometry.ro)
            / (HighPrecisionScalar(3) * point.r);
        const HighPrecisionScalar inv_r = HighPrecisionScalar(1) / point.r;
        const HighPrecisionScalar potential_hp = f_hp * (
            angular_coefficient * inv_r * inv_r
            + spin_coefficient * HighPrecisionScalar(2) * mass_hp
                  * inv_r * inv_r * inv_r);

        r[i] = point.r.convert_to<Scalar>();
        rho_cosmological[i] = point.rho_c.convert_to<Scalar>();
        f[i] = f_hp.convert_to<Scalar>();
        V[i] = potential_hp.convert_to<Scalar>();
        previous = point;
        have_previous = true;
      }
    }

    if(first_failed_index < grid_size) {
      throw std::runtime_error(
          "SdS tortoise-coordinate inversion failed at grid index "
          + std::to_string(first_failed_index));
    }

    Q.initialize(param.r_min, h, grid_size, r_cosmological,
                 r, rho_cosmological, f);
  }


  template<typename Number>
  static Number grid_coordinate(const Number &x_min, const Number &h,
                                const long long int i) {
    return sds_grid_coordinate(x_min, h, i);
  }


  Scalar grid_coordinate(const long long int i) const {
    return grid_coordinate(param.r_min, grid_space(), i);
  }


  Scalar grid_space() const {
    return (param.r_max - param.r_min) / Scalar(param.N - 1);
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

    if(Q) {
      Q(t, Q_workspace);
      if(Q_workspace.size() != grid_size) {
        throw std::invalid_argument("SdS source size does not match grid");
      }
    }
    const Scalar *__restrict__ source = Q_workspace.data();
    const Scalar d2_factor = inv_h_sqr / Scalar(12);
    const Scalar near_factor = Scalar(16) * d2_factor;
    const Scalar far_factor = -d2_factor;
    const Scalar center_coefficient = -Scalar(30) * d2_factor;
    const Scalar one_twelfth_inv_h = inv_h / Scalar(12);

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
              - Scalar(3) * pi[4]) * one_twelfth_inv_h
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
               - Scalar(36) * pi[grid_size - 3]
               + Scalar(48) * pi[n2] - Scalar(25) * pi[n1])
              * one_twelfth_inv_h - potential[n1] * psi[n1]
              + source[n1];
    dpsi[n1] = pi[n1];
  }

 private:
  struct GeometryHP {
    HighPrecisionScalar lambda;
    HighPrecisionScalar rb;
    HighPrecisionScalar rc;
    HighPrecisionScalar ro;
    HighPrecisionScalar delta;
    HighPrecisionScalar inv_fp_b;
    HighPrecisionScalar inv_fp_c;
    HighPrecisionScalar inv_fp_o;
    HighPrecisionScalar tortoise_constant;
    HighPrecisionScalar tortoise_roundoff;
  };

  struct TortoisePointHP {
    // y = log[(r-r_b)/(r_c-r)] maps the static region to the real line.
    HighPrecisionScalar y;
    HighPrecisionScalar x;
    HighPrecisionScalar dx_dy;
    HighPrecisionScalar d2x_dy2;
    HighPrecisionScalar r;
    HighPrecisionScalar rho_b;
    HighPrecisionScalar rho_c;
    bool success = false;
  };


  void validate_param() const {
    if(!(param.M > 0)) throw std::invalid_argument("SdS requires M > 0");
    if(!(param.Lambda > 0)) {
      throw std::invalid_argument("SdS requires Lambda > 0");
    }
    if(!(Scalar(9) * param.Lambda * param.M * param.M < Scalar(1))) {
      throw std::invalid_argument("SdS requires 9 Lambda M^2 < 1");
    }
    if(param.s < 0) {
      throw std::invalid_argument("SdS master spin must be nonnegative");
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


  static GeometryHP make_geometry(const HighPrecisionScalar &mass,
                                  const HighPrecisionScalar &lambda) {
    auto cubic = [&](const HighPrecisionScalar &radius) {
      return lambda * radius * radius * radius
             - HighPrecisionScalar(3) * radius
             + HighPrecisionScalar(6) * mass;
    };
    auto solve_root = [&](const HighPrecisionScalar &lower,
                          const HighPrecisionScalar &upper) {
      std::uintmax_t iterations = 1000;
      const auto bracket = boost::math::tools::toms748_solve(
          cubic, lower, upper,
          boost::math::tools::eps_tolerance<HighPrecisionScalar>(
              std::numeric_limits<HighPrecisionScalar>::digits - 8),
          iterations);
      return (bracket.first + bracket.second) / HighPrecisionScalar(2);
    };

    GeometryHP geometry;
    geometry.lambda = lambda;
    geometry.rb = solve_root(HighPrecisionScalar(2) * mass,
                             HighPrecisionScalar(3) * mass);
    geometry.rc = solve_root(HighPrecisionScalar(3) * mass,
                             sqrt(HighPrecisionScalar(3) / lambda));
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
    geometry.inv_fp_b = HighPrecisionScalar(1) / fp_b;
    geometry.inv_fp_c = HighPrecisionScalar(1) / fp_c;
    geometry.inv_fp_o = HighPrecisionScalar(1) / fp_o;

    const HighPrecisionScalar reference_r = HighPrecisionScalar(3) * mass;
    const HighPrecisionScalar reference_x = reference_r
        + HighPrecisionScalar(2) * mass * log(HighPrecisionScalar("0.5"));
    const HighPrecisionScalar reference_b = geometry.inv_fp_b
        * log(reference_r - geometry.rb);
    const HighPrecisionScalar reference_c = geometry.inv_fp_c
        * log(geometry.rc - reference_r);
    const HighPrecisionScalar reference_o = geometry.inv_fp_o
        * log(reference_r - geometry.ro);
    geometry.tortoise_constant = reference_x - reference_b - reference_c
                                  - reference_o;
    geometry.tortoise_roundoff = HighPrecisionScalar(10000)
        * std::numeric_limits<HighPrecisionScalar>::epsilon()
        * (HighPrecisionScalar(1) + abs(reference_x) + abs(reference_b)
           + abs(reference_c) + abs(reference_o));
    return geometry;
  }


  static TortoisePointHP evaluate_y(const GeometryHP &geometry,
                                    const HighPrecisionScalar &y) {
    TortoisePointHP result;
    result.y = y;
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
    const HighPrecisionScalar r_minus_ro = result.r - geometry.ro;

    result.x = geometry.tortoise_constant
               + geometry.inv_fp_b * log(result.rho_b)
               + geometry.inv_fp_c * log(result.rho_c)
               + geometry.inv_fp_o * log(r_minus_ro);

    const HighPrecisionScalar prefactor =
        HighPrecisionScalar(3) / (geometry.lambda * geometry.delta);
    result.dx_dy = prefactor * result.r / r_minus_ro;
    const HighPrecisionScalar dr_dy = result.rho_b * result.rho_c
                                      / geometry.delta;
    result.d2x_dy2 = prefactor * (-geometry.ro) * dr_dy
                     / (r_minus_ro * r_minus_ro);
    return result;
  }


  static TortoisePointHP invert_tortoise(
      const GeometryHP &geometry, const HighPrecisionScalar &target_x,
      const TortoisePointHP *previous) {
    const HighPrecisionScalar tolerance("1e-80");
    const HighPrecisionScalar target_scale =
        HighPrecisionScalar(1) + abs(target_x);
    const HighPrecisionScalar threshold = std::max(
        tolerance * target_scale, geometry.tortoise_roundoff);
    TortoisePointHP anchor = previous != nullptr && previous->success
        ? *previous : evaluate_y(geometry, HighPrecisionScalar(0));
    const HighPrecisionScalar anchor_residual = anchor.x - target_x;
    if(abs(anchor_residual) <= threshold) {
      anchor.success = true;
      return anchor;
    }
    anchor.success = false;

    // Predict the new y from dy/dx at the anchor, then expand that same step
    // until the monotonic function x(y)-target_x changes sign.
    HighPrecisionScalar step = -anchor_residual / anchor.dx_dy;
    if(step == 0) return anchor;
    const HighPrecisionScalar guess = anchor.y + step;
    TortoisePointHP outer = evaluate_y(geometry, guess);
    auto same_strict_sign = [](const HighPrecisionScalar &left,
                               const HighPrecisionScalar &right) {
      return (left < 0 && right < 0) || (left > 0 && right > 0);
    };
    HighPrecisionScalar outer_residual = outer.x - target_x;
    for(int expansion = 0;
        abs(outer_residual) > threshold
            && same_strict_sign(anchor_residual, outer_residual)
            && expansion < 100;
        ++expansion) {
      step *= HighPrecisionScalar(2);
      outer = evaluate_y(geometry, anchor.y + step);
      outer_residual = outer.x - target_x;
    }
    if(abs(outer_residual) <= threshold) {
      outer.success = true;
      return outer;
    }
    if(same_strict_sign(anchor_residual, outer_residual)) return outer;

    const HighPrecisionScalar lower = std::min(anchor.y, outer.y);
    const HighPrecisionScalar upper = std::max(anchor.y, outer.y);
    TortoisePointHP last_evaluation;
    auto root_function = [&](const HighPrecisionScalar &y) {
      last_evaluation = evaluate_y(geometry, y);
      const HighPrecisionScalar residual = last_evaluation.x - target_x;
      // Boost otherwise iterates to its step-size target even after the
      // physically relevant tortoise residual is already below threshold.
      const HighPrecisionScalar value = abs(residual) <= threshold
          ? HighPrecisionScalar(0) : residual;
      return std::make_tuple(value, last_evaluation.dx_dy,
                             last_evaluation.d2x_dy2);
    };

    try {
      std::uintmax_t max_iterations = 120;
      const int requested_bits =
          std::numeric_limits<HighPrecisionScalar>::digits - 8;
      const HighPrecisionScalar root = boost::math::tools::halley_iterate(
          root_function, guess, lower, upper, requested_bits, max_iterations);
      TortoisePointHP result = last_evaluation.y == root
          ? last_evaluation : evaluate_y(geometry, root);
      result.success = abs(result.x - target_x) <= threshold;
      return result;
    } catch(const std::exception &) {
      return outer;
    }
  }
};

#endif

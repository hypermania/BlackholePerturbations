/*! 
  \file teukolsky_precise.hpp
  \author Siyang Ling
  \brief Header for the Teukolsky equation (Schwarzschild, a=0) using high precision CPU arithmetic.

  The Teukolsky equation in Schwarzschild, converted to tortoise coordinate r_*:

    ∂_t ψ   =  Π  −  ε · D₄(ψ)
    ∂_t Π   =  −A(r)·Π  −  V(r)·ψ  +  ∂²_{r_*} ψ  −  C(r)·∂_{r_*} ψ  +  Q(t)

  where
    A(r) = 2s(r − 3M)/r²
    C(r) = 2s(r − M)/r²
    V(r) = (r−2M)/r⁴ · [r(ℓ−s)(ℓ+s+1) + 2M(1+s)]

  D₄ is the fourth-derivative Kreiss-Oliger dissipation operator.
*/
#ifndef TEUKOLSKY_PRECISE_HPP
#define TEUKOLSKY_PRECISE_HPP

#include <Eigen/Dense>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <stdexcept>
#include <utility>
#include "discretization.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


struct TeukolskyPDEPreciseParam {
  typedef boost::multiprecision::float128 Scalar;
  long long int s;
  long long int l;
  Scalar M;
  Scalar r_min;
  Scalar r_max;
  long long int N;
  Scalar ko_epsilon;
  Scalar t_start;
  Scalar t_end;
  Scalar t_interval;
  Scalar delta_t;
};


struct TeukolskyPDEPrecise {
  typedef TeukolskyPDEPreciseParam Param;
  typedef boost::multiprecision::float128 Scalar;
  typedef Eigen::Array<Scalar, -1, 1> Vector;
  typedef Vector State;
  
  Param param;
  Vector V;
  Vector A;
  Vector C;
  Vector C_d1_factor;
  Vector center_factor;
  Scalar inv_h;
  Scalar inv_h_sqr;
  long long int grid_size;
  std::function<void(const Scalar &, Vector &)> Q;
  std::function<Scalar(const Scalar &)> separable_source_time_factor;
  Vector separable_source_spatial;

  mutable Vector Q_workspace;

  void set_separable_source(Vector spatial,
                            std::function<Scalar(const Scalar &)> time_factor) {
    if(spatial.size() != grid_size) {
      throw std::invalid_argument("separable source size does not match the grid");
    }
    separable_source_spatial = std::move(spatial);
    separable_source_time_factor = std::move(time_factor);
    Q = {};
  }
  
  TeukolskyPDEPrecise(Param param_) : param(param_) {
    const Scalar r_min = param.r_min;
    const Scalar r_max = param.r_max;
    const auto N = param.N;
    const Scalar M = param.M;
    const auto s = param.s;
    const auto l = param.l;

    grid_size = N + 1;
    const Scalar h = (r_max - r_min) / (N - 1);
    inv_h = Scalar(1) / h;
    inv_h_sqr = inv_h * inv_h;

    V.resize(grid_size);
    A.resize(grid_size);
    C.resize(grid_size);
    C_d1_factor.resize(grid_size);
    center_factor.resize(grid_size);

    Vector r = compute_r_vector(r_min, r_max, N, Scalar(2)*M);
    const Scalar twoM = Scalar(2) * M;
    const Vector rm2M = r - twoM;
    const Vector inv_r = Scalar(1) / r;
    const Vector inv_r2 = inv_r * inv_r;
    const Vector inv_r4 = inv_r2 * inv_r2;

    V = rm2M * inv_r4 * (r * Scalar((l - s) * (l + s + 1)) + twoM * Scalar(1 + s));
    A = Scalar(2 * s) * (r - Scalar(3) * M) * inv_r2;
    C = Scalar(2 * s) * (r - M) * inv_r2;
    C_d1_factor = C * (inv_h / Scalar(12));
    const Scalar d2_factor = inv_h_sqr / Scalar(12);
    center_factor = -Scalar(30) * d2_factor - V;
    Q_workspace.resize(grid_size);
  }


  void operator()(const State &x, State &dxdt, const Scalar t) {
    using namespace Eigen;
    const long long int N = param.N;
    const long long int grid_size = N + 1;
    const Scalar o12 = Scalar(1) / Scalar(12);

    const auto psi = x.head(grid_size);
    const auto Pi  = x.tail(grid_size);
    auto dpsi = dxdt.head(grid_size);
    auto dPi  = dxdt.tail(grid_size);
    const Scalar *__restrict__ psi_data = x.data();
    const Scalar *__restrict__ Pi_data = x.data() + grid_size;
    Scalar *__restrict__ dpsi_data = dxdt.data();
    Scalar *__restrict__ dPi_data = dxdt.data() + grid_size;
    const Scalar *__restrict__ A_data = A.data();
    const Scalar *__restrict__ C_d1_data = C_d1_factor.data();
    const Scalar *__restrict__ center_data = center_factor.data();
    const Scalar *__restrict__ source_data = Q_workspace.data();
    const Scalar *__restrict__ source_spatial_data = separable_source_spatial.data();

    const bool has_separable_source = static_cast<bool>(separable_source_time_factor);
    const bool has_source = static_cast<bool>(Q);
    Scalar source_time_factor = 0;
    if(has_separable_source) {
      source_time_factor = separable_source_time_factor(t);
    } else if(has_source) {
      Q(t, Q_workspace);
    }

    const Scalar eps = param.ko_epsilon;
    const Scalar kofac = eps / Scalar(16);
    const Scalar d2_factor = o12 * inv_h_sqr;
    const Scalar near_factor = Scalar(16) * d2_factor;
    const Scalar far_factor = -d2_factor;

#pragma omp parallel
    {
#pragma omp for schedule(static)
      for(long long int i = 2; i <= grid_size - 3; ++i) {
        const Scalar near_sum = psi_data[i-1] + psi_data[i+1];
        const Scalar far_sum = psi_data[i-2] + psi_data[i+2];
        const Scalar gradient = Scalar(8) * (psi_data[i+1] - psi_data[i-1])
                                + psi_data[i-2] - psi_data[i+2];
        dPi_data[i] = near_factor * near_sum + far_factor * far_sum
                      + center_data[i] * psi_data[i] - C_d1_data[i] * gradient
                      - A_data[i] * Pi_data[i];
        if(has_separable_source) {
          dPi_data[i] += source_time_factor * source_spatial_data[i];
        } else if(has_source) {
          dPi_data[i] += source_data[i];
        }
        dpsi_data[i] = Pi_data[i] - kofac * (far_sum - Scalar(4) * near_sum
                                             + Scalar(6) * psi_data[i]);
      }

#pragma omp single nowait
      {
        {
          Scalar d2_0 = ( Scalar(35)*psi(0) - Scalar(104)*psi(1) + Scalar(114)*psi(2)
                         - Scalar(56)*psi(3) + Scalar(11)*psi(4) ) * o12 * inv_h_sqr;
          Scalar d1_0 = (-Scalar(25)*psi(0) + Scalar(48)*psi(1) - Scalar(36)*psi(2)
                         + Scalar(16)*psi(3) - Scalar(3)*psi(4) ) * o12 * inv_h;
          dPi(0) = d2_0 - C[0]*d1_0 - A[0]*Pi(0) - V[0]*psi(0);
          if(has_separable_source) {
            dPi(0) += source_time_factor * separable_source_spatial(0);
          } else if(has_source) dPi(0) += Q_workspace(0);
          dpsi(0) = Pi(0) - kofac * (psi(0) - Scalar(4)*psi(1) + Scalar(6)*psi(2)
                                     - Scalar(4)*psi(3) + psi(4));
        }

        {
          Scalar d2_1 = ( Scalar(11)*psi(0) - Scalar(20)*psi(1) + Scalar(6)*psi(2)
                         + Scalar(4)*psi(3) - psi(4) ) * o12 * inv_h_sqr;
          Scalar d1_1 = (-Scalar(3)*psi(0) - Scalar(10)*psi(1) + Scalar(18)*psi(2)
                         - Scalar(6)*psi(3) + psi(4) ) * o12 * inv_h;
          dPi(1) = d2_1 - C[1]*d1_1 - A[1]*Pi(1) - V[1]*psi(1);
          if(has_separable_source) {
            dPi(1) += source_time_factor * separable_source_spatial(1);
          } else if(has_source) dPi(1) += Q_workspace(1);
          dpsi(1) = Pi(1) - kofac * (psi(0) - Scalar(4)*psi(1) + Scalar(6)*psi(2)
                                     - Scalar(4)*psi(3) + psi(4));
        }

        {
          Scalar d2_n2 = ( -psi(grid_size-5) + Scalar(4)*psi(grid_size-4) + Scalar(6)*psi(grid_size-3)
                           - Scalar(20)*psi(grid_size-2) + Scalar(11)*psi(grid_size-1) ) * o12 * inv_h_sqr;
          Scalar d1_n2 = ( -psi(grid_size-5) + Scalar(6)*psi(grid_size-4) - Scalar(18)*psi(grid_size-3)
                           + Scalar(10)*psi(grid_size-2) + Scalar(3)*psi(grid_size-1) ) * o12 * inv_h;
          dPi(grid_size-2) = d2_n2 - C[grid_size-2]*d1_n2 - A[grid_size-2]*Pi(grid_size-2)
                             - V[grid_size-2]*psi(grid_size-2);
          if(has_separable_source) {
            dPi(grid_size-2) += source_time_factor * separable_source_spatial(grid_size-2);
          } else if(has_source) dPi(grid_size-2) += Q_workspace(grid_size-2);
          dpsi(grid_size-2) = Pi(grid_size-2) - kofac * (psi(grid_size-5) - Scalar(4)*psi(grid_size-4)
                                         + Scalar(6)*psi(grid_size-3) - Scalar(4)*psi(grid_size-2)
                                         + psi(grid_size-1));
        }

        {
          Scalar d2_n1 = ( Scalar(11)*psi(grid_size-5) - Scalar(56)*psi(grid_size-4)
                         + Scalar(114)*psi(grid_size-3) - Scalar(104)*psi(grid_size-2)
                         + Scalar(35)*psi(grid_size-1) ) * o12 * inv_h_sqr;
          Scalar d1_n1 = ( Scalar(3)*psi(grid_size-5) - Scalar(16)*psi(grid_size-4)
                         + Scalar(36)*psi(grid_size-3) - Scalar(48)*psi(grid_size-2)
                         + Scalar(25)*psi(grid_size-1) ) * o12 * inv_h;
          dPi(grid_size-1) = d2_n1 - C[grid_size-1]*d1_n1 - A[grid_size-1]*Pi(grid_size-1)
                             - V[grid_size-1]*psi(grid_size-1);
          if(has_separable_source) {
            dPi(grid_size-1) += source_time_factor * separable_source_spatial(grid_size-1);
          } else if(has_source) dPi(grid_size-1) += Q_workspace(grid_size-1);
          dpsi(grid_size-1) = Pi(grid_size-1) - kofac * (psi(grid_size-5) - Scalar(4)*psi(grid_size-4)
                                         + Scalar(6)*psi(grid_size-3) - Scalar(4)*psi(grid_size-2)
                                         + psi(grid_size-1));
        }
      }
    }
  }


  static Vector compute_r_vector(const Scalar r_min, const Scalar r_max,
                                  const long long int N, const Scalar r0) {
    using boost::multiprecision::cpp_bin_float_100;
    using boost::math::lambert_w0;
    typedef boost::multiprecision::cpp_bin_float_100 HighPrecisionScalar;
    HighPrecisionScalar r_min_hp = static_cast<HighPrecisionScalar>(r_min);
    HighPrecisionScalar r_max_hp = static_cast<HighPrecisionScalar>(r_max);
    HighPrecisionScalar r0_hp = static_cast<HighPrecisionScalar>(r0);
    HighPrecisionScalar h_hp = (r_max_hp - r_min_hp) / (N - 1);
    Vector r(N + 1);
#pragma omp parallel for schedule(static)
    for(long long int i = 0; i < N + 1; ++i) {
      const HighPrecisionScalar r_ast_hp = r_min_hp + i * h_hp - h_hp / HighPrecisionScalar(2);
      const HighPrecisionScalar exponent = r_ast_hp / r0_hp - HighPrecisionScalar(1);
      const HighPrecisionScalar exponential = exp(exponent);
      const HighPrecisionScalar lambert_val = lambert_w0(exponential / r0_hp);
      const HighPrecisionScalar r_hp = r0_hp * (HighPrecisionScalar(1) + lambert_val);
      r[i] = r_hp.convert_to<Scalar>();
    }
    return r;
  }


  static Vector compute_r_ast_vector(const Scalar r_min, const Scalar r_max, const long long int N) {
    Scalar h = (r_max - r_min) / (N - 1);
    Vector r_ast(N + 1);
#pragma omp parallel for schedule(static)
    for(long long int i = 0; i < N + 1; ++i) {
      r_ast[i] = r_min + i * h - h / Scalar(2);
    }
    return r_ast;
  }
};


#endif

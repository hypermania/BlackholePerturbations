/*! 
  \file regge_wheeler.hpp
  \author Siyang Ling
  \brief Header for Regge Wheeler equations that runs on the CPU.

*/
#ifndef REGGE_WHEELER_HPP
#define REGGE_WHEELER_HPP

#include <Eigen/Dense>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "discretization.hpp"

struct QuasiNormalModePDEParam {
  long long int s;
  long long int l;
  double r0;
  
  double r_min;
  double r_max;
  long long int N;
};

/*! 
  \brief Regge-Wheeler equation for s=2.
*/
struct QuasiNormalModePDE {
  typedef Eigen::VectorXd Vector;
  typedef Vector State;
  
  QuasiNormalModePDEParam param;
  Vector V;

  QuasiNormalModePDE(QuasiNormalModePDEParam param_) : param(param_) {
    using boost::multiprecision::cpp_bin_float_100;
    using boost::math::lambert_w0;
    
    const auto r_min = param.r_min;
    const auto r_max = param.r_max;
    const auto N = param.N;
    const auto r0 = param.r0;
    const auto s = param.s;
    const auto l = param.l;
    
    V.resize(N+1);
    for(int i = 0; i < V.size(); ++i) {
      const double r_ast = i_to_r_ast(r_min, r_max, N, i);
      const cpp_bin_float_100 exponent = r_ast / r0 - 1.0;
      const cpp_bin_float_100 exponential = exp(exponent);
      const cpp_bin_float_100 lambert_val = lambert_w0(exponential / cpp_bin_float_100(r0));
      const double lambert_val_double = lambert_val.convert_to<double>();
      const double r = r0 * (1.0 + lambert_val_double);
      
      // const double r = r0 * (1.0 + lambert_w0(exp(r_ast / r0 - 1.0) / r0));
      const double V_r = (1 - r0 / r) * (l*(l+1) / (r*r) + r0 * (1 - s*s) / (r*r*r));
      V[i] = V_r;
    }
  }

  /*!
    \brief The function called by odeint library.
    \param[in] x The current state of the system.
    \param[out] dxdt The time derivative, dxdt of the system.
    \param t The current time parameter.
  */
  void operator()(const State &x, State &dxdt, const double t) {
    using namespace Eigen;
    const auto N = param.N;
    const double h = (param.r_max - param.r_min) / (N - 1);
    dxdt.head(N+1) = x.tail(N+1);
    dxdt(seqN(N+1+1,N-1)) = (x(seqN(0,N-1)) + x(seqN(2,N-1)) - 2 * x(seqN(1,N-1))) / (h*h) - (V(seqN(1,N-1)).array() * x(seqN(1,N-1)).array()).matrix();
    dxdt(N+1) = 2 * (x(N+1+1) - x(N+1)) / h - (x(0) + x(2) - 2 * x(1)) / (h*h) - V(0) * x(0);
    dxdt(2*N+1) = -2 * (x(2*N+1) - x(2*N)) / h - (x(N) + x(N-2) - 2 * x(N-1)) / (h*h) - V(N) * x(N);
  }
};


#endif

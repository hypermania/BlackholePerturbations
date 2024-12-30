/*! 
  \file nonlinear_scalar.hpp
  \author Siyang Ling
  \brief Header for nonlinear scalar field QNM simulation.

*/
#ifndef NONLINEAR_SCALAR_HPP
#define NONLINEAR_SCALAR_HPP

#include <Eigen/Dense>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "discretization.hpp"

struct NonlinearScalarPDEParam {
  double r0;
  double lambda; // coeff for the \Phi^3 term
  
  double r_min;
  double r_max;
  long long int N;

  double t_start;
  double t_end;
  double t_interval;
  double delta_t;
};

// struct NonlinearScalarPDEState {
//   typedef Eigen::VectorXd Vector;
//   Vector Psi_11;
//   Vector Psi_22;
// };

/*! 
  \brief Scalar equation with a \f$ \lambda \Phi^3 \f$ interaction.
*/
struct NonlinearScalarPDE {
  typedef NonlinearScalarPDEParam Param;
  typedef Eigen::VectorXd Vector;
  typedef Vector State;
  
  Param param;
  Vector V_0;
  Vector V_1;
  Vector V_lambda;
  
  NonlinearScalarPDE(Param param_) : param(param_) {
    using boost::multiprecision::cpp_bin_float_100;
    using boost::math::lambert_w0;
    
    const auto r0 = param.r0;
    const auto lambda = param.lambda;
    const auto r_min = param.r_min;
    const auto r_max = param.r_max;
    const auto N = param.N;
    
    V_0.resize(N+1);
    V_1.resize(N+1);
    V_lambda.resize(N+1);
    
    for(int i = 0; i < V_0.size(); ++i) {
      const double r_ast = i_to_r_ast(r_min, r_max, N, i);
      const cpp_bin_float_100 exponent = r_ast / r0 - 1.0;
      const cpp_bin_float_100 exponential = exp(exponent);
      const cpp_bin_float_100 lambert_val = lambert_w0(exponential / cpp_bin_float_100(r0));
      const double lambert_val_double = lambert_val.convert_to<double>();
      const double r = r0 * (1.0 + lambert_val_double);
      
      const double V_0_r = 2 * (r - r0) / pow(r, 3) + r0 * (r - r0) / pow(r, 4);
      const double V_1_r = 6 * (r - r0) / pow(r, 3) + r0 * (r - r0) / pow(r, 4);
      const double V_lambda_r = lambda * (r0 - r) / pow(r, 2);
      V_0[i] = V_0_r;
      V_1[i] = V_1_r;
      V_lambda[i] = V_lambda_r;
    }

    // V_1.array() = 0; // Set equation to flat space
  }
  
  
  /*!
    \brief The function called by odeint library.
    \param[in] x The current state of the system.
    \param[out] dxdt The time derivative, dxdt of the system.
    \param t The current time parameter.
  */
  void operator()(const State &x, State &dxdt, const double t) {
    using namespace Eigen;
    using namespace std::numbers;
    const auto N = param.N;
    const double h = (param.r_max - param.r_min) / (N - 1);
    const auto grid_size = N + 1;

    // For our equation containing only Y_11 and Y_22
    // the convention is [Y_11, Y_22, dt_Y_11, dt_Y_22]
    dxdt.head(2*grid_size) = x.tail(2*grid_size);
    
    dxdt(seqN(2*grid_size+1,N-1)) = (x(seqN(0,N-1)) + x(seqN(2,N-1)) - 2 * x(seqN(1,N-1))) / (h*h) - (V_0(seqN(1,N-1)).array() * x(seqN(1,N-1)).array()).matrix();
    dxdt(2*grid_size) = 2 * (x(2*grid_size+1) - x(2*grid_size)) / h - (x(0) + x(2) - 2 * x(1)) / (h*h) - V_0(0) * x(0);
    dxdt(2*grid_size+N) = -2 * (x(2*grid_size+N) - x(2*grid_size+N-1)) / h - (x(N) + x(N-2) - 2 * x(N-1)) / (h*h) - V_0(N) * x(N);

    // const double turn_on_coupling = (t > 25.0) ? 1.0 : 0.0;
    const double turn_on_coupling = 1;
    
    dxdt(seqN(3*grid_size+1,N-1)) = (x(seqN(grid_size,N-1)) + x(seqN(grid_size+2,N-1)) - 2 * x(seqN(grid_size+1,N-1))) / (h*h)
      - ( V_1(seqN(1,N-1)).array() * x(seqN(grid_size+1,N-1)).array()
	  - turn_on_coupling * sqrt(3.0/(10.0*pi)) * V_lambda(seqN(1,N-1)).array() * x(seqN(1,N-1)).array() * x(seqN(1,N-1)).array() ).matrix();
    dxdt(3*grid_size) = 2 * (x(3*grid_size+1) - x(3*grid_size)) / h - (x(grid_size) + x(grid_size+2) - 2 * x(grid_size+1)) / (h*h)
      - ( V_1(0) * x(0) - turn_on_coupling * sqrt(3.0/(10.0*pi)) * V_lambda(0) * x(0) * x(0) );
    dxdt(3*grid_size+N) = -2 * (x(3*grid_size+N) - x(3*grid_size+N-1)) / h - (x(grid_size+N) + x(grid_size+N-2) - 2 * x(grid_size+N-1)) / (h*h)
      - ( V_1(N) * x(N) - turn_on_coupling * sqrt(3.0/(10.0*pi)) * V_lambda(N) * x(N) * x(N) );

    // dxdt.head(N+1) = x.tail(N+1);
    // dxdt(seqN(N+1+1,N-1)) = (x(seqN(0,N-1)) + x(seqN(2,N-1)) - 2 * x(seqN(1,N-1))) / (h*h) - (V(seqN(1,N-1)).array() * x(seqN(1,N-1)).array()).matrix();
    // dxdt(N+1) = 2 * (x(N+1+1) - x(N+1)) / h - (x(0) + x(2) - 2 * x(1)) / (h*h) - V(0) * x(0);
    // dxdt(2*N+1) = -2 * (x(2*N+1) - x(2*N)) / h - (x(N) + x(N-2) - 2 * x(N-1)) / (h*h) - V(N) * x(N);
  }
};


#endif

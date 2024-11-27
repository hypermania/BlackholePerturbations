// #include <cstdlib>
// #include <cmath>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
#include <vector>
// #include <algorithm>
// #include <chrono>
// #include <string>
// #include <filesystem>

#include <Eigen/Dense>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

#include "utility.hpp"
#include "param.hpp"
// #include "initializer.hpp"
// #include "random_field.hpp"
#include "io.hpp"
// #include "fdm3d.hpp"
// #include "physics.hpp"
// #include "equations.hpp"
// #include "workspace.hpp"
#include "observer.hpp"
// #include "midpoint.hpp"
// #include "wkb.hpp"

// #ifndef DISABLE_CUDA
// #include <thrust/device_vector.h>
// #include "equations_cuda.cuh"
// #include "cufft.h"
// #include "fdm3d_cuda.cuh"
// #endif

struct QuasiNormalModePDEParam {
  long long int s;
  long long int l;
  double r0;
  
  double r_min;
  double r_max;
  long long int N;
};

double i_to_r_ast(const double r_min, const double r_max, const long long int N, const long long int i) {
  const double h = (r_max - r_min) / (N - 1);
  return r_min + i * h - h / 2.0;
}

/*! 
  \brief QNM equation
*/
struct QuasiNormalModePDE {
  typedef Eigen::VectorXd Vector;
  typedef Vector State;
  
  QuasiNormalModePDEParam param;
  Vector V;

  QuasiNormalModePDE(QuasiNormalModePDEParam param_) : param(param_) {
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
      const double r = r0 * (1.0 + lambert_w0(exp(r_ast / r0 - 1.0) / r0));
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


int main(int argc, char **argv){
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using boost::math::lambert_w0;

  const std::string dir = "output/initial_test/";
  prepare_directory_for_output(dir);

  const long long int s = 2;
  const long long int l = 2;
  const double r0 = 1;
  
  const double r_min = -100;
  const double r_max = 100;
  const long long int N = 1 << 12;
  
  typedef QuasiNormalModePDE Equation;
  typedef QuasiNormalModePDEParam Param;
  typedef Equation::State State;

  Param param;
  param.r_min = r_min;
  param.r_max = r_max;
  param.N = N;
  param.s = s;
  param.l = l;
  param.r0 = r0;
  
  // The equation object.
  Equation eqn(param);

  // Choose what to save in the course of simulation.
  // Here we save the field spectrum, density spectrum, and 2D density slices.
  std::vector<double> t_list;
  std::vector<State> x_list;
  SaveAllObserver<State> observer(t_list, x_list);

  
  // Choose the numerical integrator.
  // Here we use RK4, you can also use other methods.
  // See https://www.boost.org/doc/libs/1_85_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/overview.html .
  auto stepper = runge_kutta4_classic<State, double, State, double>();
  // auto stepper = make_controlled(1e-9, 1e-9, runge_kutta_fehlberg78<State, double, State, double>());
  
  VectorXd state(2 * (N+1));
  state.array() = 0;
  for(int i = 0; i < N+1; ++i) {
    const double r_ast = i_to_r_ast(r_min, r_max, N, i);
    state(i) = exp(-2.*(r_ast - 50.)*(r_ast - 50.));
  }

  const double t_start = 0;
  const double t_end = 160;
  const double delta_t = 0.01;
  
  // Solve the equation.
  run_and_measure_time("Solving equation",
  		       [&](){
			 int num_steps = integrate_const(stepper, eqn, state, t_start, t_end, delta_t, observer);
			 std::cout << "total number of steps = " << num_steps << '\n';
		       } );
  
  write_vector_to_file(t_list, dir + "t_list.dat");
  // write_vector_to_file(x_list, dir + "t_list.dat");
  for(int i = 0; i < x_list.size(); ++i) {
    write_VectorXd_to_filename_template(x_list[i], dir + "state_%d.dat", i);
  }
  // State &last = x_list.back();
  // write_data_to_file((const char *)last.data(), last.size() * sizeof(double), dir + "last.dat");
  save_param_for_Mathematica(param, dir);

    // VectorXd V(N+1);
  // for(int i = 0; i < V.size(); ++i) {
  //   const double r_ast = i_to_r_ast(r_min, r_max, N, i);
  //   const double r = r0 * (1.0 + lambert_w0(exp(r_ast / r0 - 1.0) / r0));
  //   const double V_r = (1 - r0 / r) * (l*(l+1) / (r*r) + r0 * (1 - s*s) / (r*r*r));
  //   V[i] = V_r;
  // }

  // for(int i = 0; i < V.size(); ++i) {
  //   std::cout << i_to_r_ast(r_min, r_max, N, i) << ", " << V[i] << '\n';
  // }
  
  // for(int i = 0; i < 10; ++i) {
  //   const double lambert_arg = i;
  //   const double lambert_val = lambert_w0(lambert_arg);
  //   std::cout << lambert_val << '\n';
  // }

}  



/*
void solve_field_equation(void)
{
  using namespace Eigen;
  using namespace boost::numeric::odeint;

  
  // Set the PRNG seed.
  RandomNormal::set_generator_seed(0);

  
  // Set the directory for output.
  const std::string dir = "output/Growth_and_FS/";
  prepare_directory_for_output(dir);

  
  // Set parameters for the simulation.
  MyParam param
    {
     .N = 384, // Lattice points per axis
     .L = 384 * 0.8, // Size of the box
     // ULDM params
     .m = 1.0, // Mass of scalar field
     .lambda = 0, // Lambda phi^4 coupling strength
     //.f_a = 30.0, // Not relevant for ComovingCurvatureEquationInFRW
     .k_ast = 1.0, // Characteristic momentum
     .k_Psi = 1.0, // Not relevant for ComovingCurvatureEquationInFRW
     .varphi_std_dev = 1.0, // Standard deviation of field
     .Psi_std_dev = 0.02, // Standard deviation of metric perturbation Psi
     // FRW metric params
     .a1 = 1.0,
     .H1 = 0.05,
     .t1 = 1.0 / (2 * param.H1),
     // Start and end time for numerical integration, and time interval between saves
     .t_start = param.t1,
     .t_end = param.t_start + (pow(3.5 / param.a1, 2) - 1.0) / (2 * param.H1),
     .t_interval = 49.99, // Save a snapshot every t_interval
     // Numerical method parameter
     .delta_t = 0.5, // Time step for numerical integration
     // Psi approximation parameter
     .M = 128 // Lattice points for storing / computing Psi
    };
  print_param(param);
  save_param_for_Mathematica(param, dir);

  
  // Choose an equation to solve.
  // Here we solve a scalar field equation with background metric perturbations.
  // Also see CudaApproximateComovingCurvatureEquationInFRW, which is a CUDA implementation of the same equation.
  typedef ComovingCurvatureEquationInFRW Equation;
  //typedef CudaApproximateComovingCurvatureEquationInFRW Equation;
  typedef typename Equation::Workspace Workspace;
  typedef typename Equation::State State;

  
  // Initialize the workspace given params and a procedure for setting initial conditions.
  // The initialization procedure is described in Sec.3 of the paper.
  Workspace workspace(param, perturbed_grf_and_comoving_curvature_fft);

  
  // The equation object.
  Equation eqn(workspace);

  
  // Choose what to save in the course of simulation.
  // Here we save the field spectrum, density spectrum, and 2D density slices.
  ConstIntervalObserver<Equation, true, true, true> observer(dir, param, eqn);

  
  // Choose the numerical integrator.
  // Here we use RK4, you can also use other methods.
  // See https://www.boost.org/doc/libs/1_85_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/overview.html .
  auto stepper = runge_kutta4_classic<State, double, State, double>();
  // auto stepper = make_controlled(1e-9, 1e-9, runge_kutta_fehlberg78<State, double, State, double>());

  
  {
    // Save spectrum for R and initial potential Psi
    double eta_i = workspace.cosmology.eta(param.t_start);
    auto kernel = [eta_i](double k){
		    return k == 0.0 ? 0.0 : (6 * sqrt(3) * (-((k * eta_i * cos((k * eta_i) / sqrt(3))) / sqrt(3)) + sin((k * eta_i) / sqrt(3)))) / (pow(k, 3) * pow(eta_i, 3));
		  };
    
    Eigen::VectorXd R_fft_eigen(workspace.R_fft.size());
    copy_vector(R_fft_eigen, workspace.R_fft);
    
    auto fft_wrapper = fftwWrapper(param.N);
    Eigen::VectorXd R = fft_wrapper.execute_z2d(R_fft_eigen) / pow(param.N, 3);
    Eigen::VectorXd Psi = compute_field_with_scaled_fourier_modes(param.N, param.L, R, kernel, fft_wrapper);

    std::cout << "Psi_std_dev = " << sqrt(Psi.squaredNorm() / pow(param.N, 3)) << '\n';
    auto Psi_spectrum = compute_power_spectrum(param.N, Psi, fft_wrapper);
    write_VectorXd_to_file(Psi_spectrum, dir + "initial_Psi_spectrum.dat");
    
    auto R_spectrum = compute_power_spectrum(param.N, R, fft_wrapper);
    write_VectorXd_to_file(R_spectrum, dir + "initial_R_spectrum.dat");
  }

  
  // Solve the equation.
  run_and_measure_time("Solving equation",
  		       [&](){
			 int num_steps = integrate_const(stepper, eqn, workspace.state, param.t_start, param.t_end, param.delta_t, observer);
			 std::cout << "total number of steps = " << num_steps << '\n';
		       } );
  
  write_vector_to_file(workspace.t_list, dir + "t_list.dat");
  
  
  // Optional: save the final state.
  {
    Eigen::VectorXd state_out(workspace.state.size());
    copy_vector(state_out, workspace.state);
    write_VectorXd_to_file(state_out, dir + "state.dat");
  }
}

*/

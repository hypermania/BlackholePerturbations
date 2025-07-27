#include <vector>
#include <algorithm>
#include <typeinfo>
//#include <thread>

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include "odeint_eigen/eigen_operations.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
//#include <boost/lockfree/queue.hpp>

#include "utility.hpp"
#include "param.hpp"
#include "io.hpp"
#include "observer.hpp"
//#include "regge_wheeler.hpp"
//#include "regge_wheeler_precise.hpp"
//#include "cubic_scalar_simplified.hpp"
#include "cubic_scalar.hpp"
#include "boost/type_index.hpp"

#include "teukolsky_scalar.hpp"

#include "examples.hpp"
#include "rsh.hpp"

void run_teukoksky_benchmark(void);

int main(int argc, char **argv) {
  //run_coupled_eqn();
  
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using namespace std::numbers;
  using std::array;

  const std::string dir = "output/test_teukolsky/";
  prepare_directory_for_output(dir);

  typedef TeukolskyScalarPDE Equation;
  typedef Equation::Param Param;
  typedef Equation::State State;
  
  Param param;
  param.s = 0; //.convert_to<double>();
  param.l_max = 3; //.convert_to<double>();
  param.M = 0.5;
  param.a = 0.49;

  param.rast_min = -50;
  param.rast_max = 75;
  param.N = static_cast<long long int>((param.rast_max - param.rast_min) / 0.03); //1000;
  
  param.t_start = 0; //.convert_to<double>();
  param.t_end = 50; //.convert_to<double>();
  param.t_interval = 0.5;
  param.delta_t = 0.01; //.convert_to<double>();

  save_param_for_Mathematica(param, dir);

  Equation eqn(param);

  // {
  //   for(auto [lm1, idx1] : eqn.drdr_psi_lm_map[3]){
  //     std::cout << "lm1, idx1 = " << lm1 << ", " << idx1 << std::endl;
  //     std::cout << eqn.coeffs[idx1](Eigen::seqN(2000, 50)).transpose() << std::endl;
  //   }
  // }
  // exit(0);
  
  auto stepper = runge_kutta_fehlberg78<State, double, State, double>();

  // const long long int rIdx = r_ast_to_i(r_min, r_max, N, 50.0);
  // std::cout << "using rIdx = " << rIdx << std::endl;
  // std::vector<long long int> positions;
  // for(int i = 0; i < 2 * eqn.lm_size; ++i) {
  //   positions.push_back(eqn.grid_size * i + rIdx);
  // }
  // auto observer1 = FixedPositionObserver(dir, positions);
  // auto observer2 = ApproximateTimeObserver(dir, {50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550.});
  // auto observer = ObserverPack(observer1, observer2);
  
  auto observer = ApproximateTimeObserver(dir, {10., 20., 30., 40.});

  State state(2 * eqn.lm_size * eqn.grid_size);
  state = 0;

  ArrayXd r_ast(eqn.grid_size);
  for(int i = 0; i < eqn.grid_size; ++i) {
    r_ast[i] = i_to_r_ast(param.rast_min, param.rast_max, param.N, i);
  }

  const double r_source = 25;
  const double sigma = 0.5;
  const long long int grid_begin = RSH::lm_to_idx(1, 1) * eqn.grid_size;
  state(seqN(grid_begin, eqn.grid_size)) = pow(2 * pi, -0.5) * (1 / sigma) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma));
  state(seqN(eqn.grid_size * eqn.lm_size + grid_begin, eqn.grid_size)) = - pow(2 * pi, -0.5) * pow(sigma, -3) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma)) * (r_ast - r_source);

  // {
  //   TeukolskyScalarPDE::ComplexVector dr_psi_lm(eqn.lm_size * eqn.grid_size);
  //   TeukolskyScalarPDE::ComplexVector drdr_psi_lm(eqn.lm_size * eqn.grid_size);
  //   eqn.compute_derivatives(state, dr_psi_lm, drdr_psi_lm);
  //   std::cout << drdr_psi_lm(Eigen::seqN(RSH::lm_to_idx(1, 1) * eqn.grid_size + 2500, 50)).transpose() << std::endl;
  //   std::cout << dr_psi_lm(Eigen::seqN(RSH::lm_to_idx(1, 1) * eqn.grid_size + 2500, 50)).transpose() << std::endl;
  // }
  //     exit(0);

  //exit(0);
  
  // Solve the equation.
  run_and_measure_time("Solving equation",
  		       [&](){
			 // int num_steps = integrate_const(stepper, std::ref(eqn), state, param.t_start, param.t_end, param.delta_t, std::ref(observer));
			 int num_steps = integrate_adaptive(stepper, std::ref(eqn), state, param.t_start, param.t_end, param.delta_t, std::ref(observer));
			 std::cout << "total number of steps = " << num_steps << '\n';
		       } );
  write_to_file(state, dir + "final_state.dat");

  observer.save();
  
  // typedef TeukolskyScalarPDE::HighPrecisionScalar HP;  
  // auto c = std::complex<HP>(HP(1.0), HP(1.0));
  
  //run_coupled_eqn();
  //run_sourced_eqn();

  return 0;
}

void run_teukoksky_benchmark(void) {
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using namespace std::numbers;
  using std::array;
    
  //const std::string dir = "output/quadratic_rsh_coupling_0001_ingoing/";
  const std::string dir = "output/test_teukolsky_benchmark/";
  prepare_directory_for_output(dir);

  const double r0 = 1;
  const long long int l_max = 2;  // The cutoff angular number
  const double lambda = 0;
  const double r_min = -50; //-400;
  const double r_max = 75; //600;
  //const long long int N = 1 << 15;
  const long long int N = static_cast<long long int>((r_max - r_min) / 0.03);
  
  const double r_source = 25;
  
  typedef CubicScalarPDE<l_max> Equation;
  typedef Equation::Param Param;
  typedef Equation::State State;
  
  Param param;
  param.r0 = r0;
  param.l_max = l_max;
  param.lambda = lambda;
  param.r_min = r_min;
  param.r_max = r_max;
  param.N = N;
  param.t_start = 0;
  param.t_end = 50; //r_max - r_source;
  param.t_interval = 0.5;
  param.delta_t = 0.01;

  save_param_for_Mathematica(param, dir);
  
  // The equation object.
  Equation eqn(param);
  
  // Choose the numerical integrator.
  // auto stepper = make_controlled(1e-15, 1e-15, runge_kutta_fehlberg78<State, double, State, double>());
  auto stepper = runge_kutta_fehlberg78<State, double, State, double>();


  // Extract the waveform at r_* = 50
  const long long int rIdx = r_ast_to_i(r_min, r_max, N, 50.0);
  std::cout << "using rIdx = " << rIdx << std::endl;
  std::vector<long long int> positions;
  for(int i = 0; i < 2 * eqn.lm_size; ++i) {
    positions.push_back(eqn.grid_size * i + rIdx);
  }
  auto observer1 = FixedPositionObserver(dir, positions);
  auto observer2 = ApproximateTimeObserver(dir, {10., 20., 30., 40.});
  auto observer = ObserverPack(observer1, observer2);

  
  // Initialize in-going wavepacket
  State state(eqn.state_size);
  state = 0;
  
  ArrayXd r_ast(eqn.grid_size);
  for(int i = 0; i < eqn.grid_size; ++i) {
    r_ast[i] = i_to_r_ast(r_min, r_max, N, i);
  }

  const double sigma = 0.5;
  const long long int grid_begin = RSH::lm_to_idx(1, 1) * eqn.grid_size;
  state(seqN(grid_begin, eqn.grid_size)) = pow(2 * pi, -0.5) * (1 / sigma) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma));
  state(seqN(eqn.half_state_size + grid_begin, eqn.grid_size)) = - pow(2 * pi, -0.5) * pow(sigma, -3) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma)) * (r_ast - r_source);
  
  // Solve the equation.
  run_and_measure_time("Solving equation",
  		       [&](){
			 // int num_steps = integrate_const(stepper, std::ref(eqn), state, param.t_start, param.t_end, param.delta_t, std::ref(observer));
			 int num_steps = integrate_adaptive(stepper, std::ref(eqn), state, param.t_start, param.t_end, param.delta_t, std::ref(observer));
			 std::cout << "total number of steps = " << num_steps << '\n';
		       } );
  write_to_file(state, dir + "final_state.dat");
  observer.save();

}

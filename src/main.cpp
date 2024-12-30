#include <vector>
#include <algorithm>
#include <typeinfo>

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include "odeint_eigen/eigen_operations.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>

#include "utility.hpp"
#include "param.hpp"
#include "io.hpp"
#include "observer.hpp"
#include "regge_wheeler.hpp"
#include "cubic_scalar_simplified.hpp"
#include "cubic_scalar.hpp"

#include "boost/type_index.hpp"

// #ifndef DISABLE_CUDA
// #include <thrust/device_vector.h>
// #include "equations_cuda.cuh"
// #include "cufft.h"
// #include "fdm3d_cuda.cuh"
// #endif


void run_quadratic_eqn(void);
void run_full_eqn(void);

struct FixedPositionObserver {
  typedef Eigen::ArrayXd State;
  std::string dir;
  std::vector<long long int> positions;
  std::vector<double> t_list;
  std::vector<double> psi_list;

  FixedPositionObserver(const std::string &dir_, const std::vector<long long int> positions_) :
    dir(dir_), positions(positions_)
  {}
  
  void operator()(const State &x, const double t) {
    for(auto i : positions) {
      psi_list.push_back(x[i]);
    }
    t_list.push_back(t);
  }

  void save(void) const {
    write_to_file(psi_list, dir + "psi_list.dat");
    write_to_file(t_list, dir + "t_list.dat");    
  }
};


template<typename... Observers>
struct ObserverPack {
  typedef Eigen::ArrayXd State;
  std::tuple<Observers & ...> observers;
  
  ObserverPack(Observers & ... observers_) : observers(observers_...) {}

  void operator()(const State &x, const double t) {
    std::apply([&](auto &&... args) { ((args(x, t)), ...); }, observers);
  }
};


int main(int argc, char **argv) {
  // auto obj = Eigen::ArrayXd::Zero(10);
  // std::cout << boost::typeindex::type_id_runtime(obj) << '\n';


  // constexpr long long int lm_0 = std::get<0>(sqr_rsh_terms[81]);
  // constexpr long long int l_0 = RSH::idx_to_lm(lm_0).first;

  // std::cout << "lm_0 = " << lm_0 << std::endl
  // 	    << "l_0 = " << l_0 << std::endl;
  
  //const long long int lm_size = sqr_rsh_terms.size();
  //const long long int grid_size = 10;
  //Eigen::ArrayXd x = Eigen::ArrayXd::Ones(lm_size * grid_size);
  
  //auto expr1 = RSH::quadratic_rsh_expression<8, 40>(x, grid_size);
  //std::cout << boost::typeindex::type_id_runtime(expr1) << '\n';
  
  // std::vector<double> t;
  // std::vector<Eigen::ArrayXd> f;
  // auto obs1 = SaveAllObserver<Eigen::ArrayXd>(t, f);
  // auto obs2 = SaveAllObserver<Eigen::ArrayXd>(t, f);

  Eigen::ArrayXd test = Eigen::ArrayXd::Ones(10);
  auto obs1 = FixedPositionObserver("", {1});
  auto obs2 = FixedPositionObserver("", {2});
  
  auto obs = ObserverPack(obs1, obs2, obs1, obs2, obs1);
  std::cout << boost::typeindex::type_id_runtime(obs) << '\n';
  obs(test, 0);

  std::cout << obs1.t_list.size() << std::endl
	    << obs2.t_list.size() << std::endl;
  
  // run_quadratic_eqn();
  // run_full_eqn();
  
  return 0;
}


void run_quadratic_eqn(void) {
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using namespace std::numbers;
  using std::array;

  // const std::string dir = "/home/hypermania/Research/BHQuasinormalModes/output/outgoing_wavepacket_flat_2/";
  const std::string dir = "/home/hypermania/Research/BHQuasinormalModes/output/speed_test/";
  prepare_directory_for_output(dir);

  const double r0 = 1;
  const double lambda = 0.1;
  const double r_min = -400;
  const double r_max = 600;
  const long long int N = 1 << 15;
  
  const double r_source = 50;
  
  typedef NonlinearScalarPDE Equation;
  typedef Equation::Param Param;
  typedef Equation::State State;

  Param param;
  param.r0 = r0;
  param.lambda = lambda;
  param.r_min = r_min;
  param.r_max = r_max;
  param.N = N;
  param.t_start = 0;
  param.t_end = 20; //r_max - r_source;
  param.t_interval = 0.5;
  param.delta_t = 0.01;

  save_param_for_Mathematica(param, dir);
  
  // The equation object.
  Equation eqn(param);

  // Save constant intervals
  // ConstIntervalObserver observer(dir, param);
  auto observer = [](const State &x, double t){};

  // auto stepper = runge_kutta4_classic<State, double, State, double>();
  // auto stepper = adams_bashforth<5, State, double, State, double>();
  // auto stepper = adams_bashforth_moulton<5, State, double, State, double>();
  auto stepper = make_controlled(1e-15, 1e-15, runge_kutta_fehlberg78<State, double, State, double>());
  

  VectorXd state(4 * (N+1));
  state.array() = 0;
  const double sigma = 0.5;
  for(int i = 0; i < N+1; ++i) {
    const double r_ast = i_to_r_ast(r_min, r_max, N, i);
    // Out-going wavepacket
    state(i) = pow(2 * pi, -0.5) * (1 / sigma) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma));
    state(2*(N+1)+i) = pow(2 * pi, -0.5) * pow(sigma, -3) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma)) * (r_ast - r_source);
  }

  // Solve the equation.
  run_and_measure_time("Solving equation",
  		       [&](){
			 // int num_steps = integrate_const(stepper, eqn, state, t_start, t_end, delta_t, observer);
			 int num_steps = integrate_adaptive(stepper, std::ref(eqn), state, param.t_start, param.t_end, param.delta_t, std::ref(observer));
			 std::cout << "total number of steps = " << num_steps << '\n';
		       } );
  
  // write_to_file(observer.t_list, dir + "t_list.dat");
}


void run_linear_eqn(void){

  
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using boost::math::lambert_w0;
  using namespace std::numbers;

  // const std::string dir = "output/initial_test/";
  const std::string dir = "/home/hypermania/Research/BHQuasinormalModes/output/single_idx_test/";
  prepare_directory_for_output(dir);

  const long long int s = 2;
  const long long int l = 2;
  const double r0 = 1;
  
  const double r_min = -1000;
  const double r_max = 1000;
  const long long int N = 1 << 15;
  
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

  save_param_for_Mathematica(param, dir);
  
  // The equation object.
  Equation eqn(param);

  // Choose what to save in the course of simulation.
  // Here we save the field spectrum, density spectrum, and 2D density slices.
  std::vector<double> t_list;
  std::vector<State> x_list;
  
  // SaveAllObserver<State> observer(t_list, x_list);
  const long long int rIdx = r_ast_to_i(r_min, r_max, N, 50.0);
  std::cout << "using rIdx = " << rIdx << std::endl;
  std::vector<double> psi_list;
  std::vector<double> dt_psi_list;
  auto observer = [&](const State &x, double t)->void{
    psi_list.push_back(x[rIdx]);
    dt_psi_list.push_back(x[(N+1) + rIdx]);
    t_list.push_back(t);
  };

  
  // Choose the numerical integrator.
  // Here we use RK4, you can also use other methods.
  // See https://www.boost.org/doc/libs/1_85_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/overview.html .
  // auto stepper = runge_kutta4_classic<State, double, State, double>();
  auto stepper = make_controlled(1e-20, 1e-20, runge_kutta_fehlberg78<State, double, State, double>());

  const double r_source = 50.;
  VectorXd state(2 * (N+1));
  state.array() = 0;
  for(int i = 0; i < N+1; ++i) {
    const double r_ast = i_to_r_ast(r_min, r_max, N, i);
    state(i) = pow(pi, -0.5) * exp(-(r_ast - r_source)*(r_ast - r_source));
  }

  const double t_start = 0;
  const double t_end = 2 * (r_max - r_source); //1000; // t < 2 * (r_max-r_source)
  const double delta_t = 0.001;
  
  // Solve the equation.
  run_and_measure_time("Solving equation",
  		       [&](){
			 // int num_steps = integrate_const(stepper, eqn, state, t_start, t_end, delta_t, observer);
			 int num_steps = integrate_adaptive(stepper, eqn, state, t_start, t_end, delta_t, observer);
			 std::cout << "total number of steps = " << num_steps << '\n';
		       } );
  
  write_to_file(t_list, dir + "t_list.dat");
  // for(int i = 0; i < x_list.size(); ++i) {
  //   write_VectorXd_to_filename_template(x_list[i], dir + "state_%d.dat", i);
  // }
  write_to_file(psi_list, dir + "psi_list.dat");
  write_to_file(dt_psi_list, dir + "dt_psi_list.dat");
  
}  




void run_full_eqn(void) {
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using namespace std::numbers;
  using std::array;
    
  // const std::string dir = "/home/hypermania/Research/BHQuasinormalModes/output/animation_test_14/";
  const std::string dir = "/home/hypermania/Dropbox/Research/BHQuasinormalModes/QuasiNormalModes/output/quadratic_rsh_coupling_0001/";
  prepare_directory_for_output(dir);

  const double r0 = 1;
  const long long int l_max = 0;
  const double lambda = 0.001;
  const double r_min = -400;
  const double r_max = 600;
  const long long int N = 1 << 15;
  
  const double r_source = 50;
  
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
  param.t_end = 800; //r_max - r_source;
  param.t_interval = 0.5;
  param.delta_t = 0.001;

  save_param_for_Mathematica(param, dir);
  
  // The equation object.
  Equation eqn(param);
  
  // Choose the numerical integrator.
  auto stepper = make_controlled(1e-15, 1e-15, runge_kutta_fehlberg78<State, double, State, double>());

  const long long int rIdx = r_ast_to_i(r_min, r_max, N, 50.0);
  std::cout << "using rIdx = " << rIdx << std::endl;
  std::vector<long long int> positions;
  for(int i = 0; i < 2 * eqn.lm_size; ++i) {
    positions.push_back(eqn.grid_size * i + rIdx);
  }
  auto observer = FixedPositionObserver(dir, positions);

  /*
  std::vector<double> t_list;
  std::vector<double> psi_list;
  auto observer = [&](const State &x, double t)->void{
    for(int i = 0; i < eqn.lm_size * 2; ++i) {
      psi_list.push_back(x[eqn.grid_size * i + rIdx]);
    }
    t_list.push_back(t);
  };
  */

  
  // Initialize out-going wavepacket  
  State state(eqn.state_size);
  state = 0;
  
  ArrayXd r_ast(eqn.grid_size);
  for(int i = 0; i < eqn.grid_size; ++i) {
    r_ast[i] = i_to_r_ast(r_min, r_max, N, i);
  }

  const double sigma = 0.5;
  const long long int grid_begin = RSH::lm_to_idx(1, 1) * eqn.grid_size;
  state(seqN(grid_begin, eqn.grid_size)) = pow(2 * pi, -0.5) * (1 / sigma) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma));
  state(seqN(eqn.half_state_size + grid_begin, eqn.grid_size)) = pow(2 * pi, -0.5) * pow(sigma, -3) * exp(-(r_ast - r_source)*(r_ast - r_source) / (2 * sigma * sigma)) * (r_ast - r_source);
  
  // Solve the equation.
  run_and_measure_time("Solving equation",
  		       [&](){
			 // int num_steps = integrate_const(stepper, eqn, state, t_start, t_end, delta_t, observer);
			 int num_steps = integrate_adaptive(stepper, std::ref(eqn), state, param.t_start, param.t_end, param.delta_t, std::ref(observer));
			 std::cout << "total number of steps = " << num_steps << '\n';
		       } );
  write_to_file(state, dir + "final_state.dat");
  // write_to_file(psi_list, dir + "psi_list.dat");
  // write_to_file(t_list, dir + "t_list.dat");
  observer.save();

}

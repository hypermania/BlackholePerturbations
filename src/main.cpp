#include <vector>
#include <algorithm>
#include <typeinfo>
#include <thread>

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include "odeint_eigen/eigen_operations.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/lockfree/queue.hpp>

#include "utility.hpp"
#include "param.hpp"
#include "io.hpp"
#include "observer.hpp"
//#include "regge_wheeler.hpp"
#include "regge_wheeler_precise.hpp"
//#include "cubic_scalar_simplified.hpp"
#include "cubic_scalar.hpp"

#include "boost/type_index.hpp"
#include <boost/math/tools/roots.hpp>

void run_sourced_eqn(void);
void run_coupled_eqn(void);
void solve_rast(void);


int main(int argc, char **argv) {  
  // run_coupled_eqn();
  //run_sourced_eqn();
  solve_rast();
  return 0;
}


/*! 
  \brief Solve a bunch of sourced Regge-Wheeler equations in parallel.
*/
void run_sourced_eqn(void) {
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using boost::math::lambert_w0;
  using boost::multiprecision::cpp_bin_float_100;
  using namespace std::numbers;

  typedef QuasiNormalModePDEPrecise Equation;
  typedef QuasiNormalModePDEPreciseParam Param;
  typedef Equation::Scalar Scalar;
  typedef Equation::State State;
  typedef Equation::Vector Vector;

  // Solves the Regge-Wheeler equation for a particular set of (l,\beta),
  // where l and beta are the angular number and the fall off rate of the source.
  auto run_simulation = [](const long long int l, const long long int beta)->void {
    std::string format_string = "output/batched_precise/l_%d_beta_%d/";
      
    char dir_buffer[128];
    sprintf(dir_buffer, format_string.data(), l, beta);
    const std::string dir(dir_buffer);
    prepare_directory_for_output(dir);  
  
    const long long int s = 0;
    // const long long int l = 1;
    const Scalar r0 = 1;
  
    const Scalar r_min = -600; // -400
    const Scalar r_max = 1200; // -800
    const long long int N = static_cast<long long int>((r_max - r_min) / 0.03);

    const Scalar t_start = 0;
    const Scalar t_end = 1000;
    const Scalar delta_t = 0.01;

    Param param;
    param.r_min = r_min; //.convert_to<double>();
    param.r_max = r_max; //.convert_to<double>();
    param.N = N;
    param.s = s;
    param.l = l;
    param.r0 = r0; //.convert_to<double>();
    param.t_start = t_start; //.convert_to<double>();
    param.t_end = t_end; //.convert_to<double>();
    param.t_interval = 0.5;
    param.delta_t = delta_t; //.convert_to<double>();

    save_param_for_Mathematica(param, dir);
  
    Equation eqn(param);
  
    auto stepper = runge_kutta_fehlberg78<State, Scalar, State, Scalar>();
    
    // Extract the waveform at r_* = 50  
    const long long int rIdx = r_ast_to_i(param.r_min.convert_to<double>(), param.r_max.convert_to<double>(), N, 50.0);
    auto observer1 = FixedPositionObserver(dir, {rIdx, rIdx + (N+1)});
    auto observer2 = ApproximateTimeObserver(dir, {50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550.});
    auto observer = ObserverPack(observer1, observer2);

    Vector r_ast = eqn.compute_r_ast_vector(r_min, r_max, N);
    Vector r = eqn.compute_r_vector(r_min, r_max, N, r0);
      
    const Scalar sigma = Scalar(1) / Scalar(2);
    const Scalar prefactor = pow(Scalar(2 * pi), Scalar(-0.5)) * (1 / sigma);
    const Scalar exp_factor = Scalar(1) / (2 * sigma * sigma);

    Vector front_factor = r;
    front_factor = front_factor.pow(-beta);
    front_factor *= prefactor;
    front_factor.head(r_ast_to_i(param.r_min.convert_to<double>(), param.r_max.convert_to<double>(), N, 10.0)) = 0;

    // An outgoing Gaussian source
    eqn.Q = [&](const Scalar t)->Vector{
      // const Scalar sigma = Scalar(1) / Scalar(2);
      // return pow(Scalar(2 * pi), Scalar(-0.5)) * (1 / sigma) * exp(-(t - r_ast).abs2() / (2 * sigma * sigma)) * r_factor;
      //return prefactor * exp(-(t - r_ast).abs2() * exp_factor) * r_factor;
      return front_factor * exp(-(t - r_ast).abs2() * exp_factor);
    };
      
    Vector state = Vector::Zero(2 * (N+1));
      
    // Solve the equation.
    run_and_measure_time("Solving equation",
			 [&](){
			   int num_steps = integrate_adaptive(stepper, std::ref(eqn), state, t_start, t_end, delta_t, std::ref(observer));
			   std::cout << "total number of steps = " << num_steps << '\n';
			 } );
    observer.save();
  };


  // Solve the equations in parallel by calling run_simulation in different threads
  // l_beta_array is a queue of (l,\beta) parameters to solve
  std::vector<std::pair<long long int, long long int>> l_beta_array;
  // for(long long int l = 4; l <= 4; ++l) {
  //   for(long long int beta = 2; beta <= 6; ++beta) {
  //     l_beta_array.push_back(std::make_pair(l, beta));
  //   }
  // }
  l_beta_array.push_back(std::make_pair(4, 3));
	
  boost::lockfree::queue<int> q(10);
  for(size_t idx = 0; idx < l_beta_array.size(); ++idx){
    q.push(idx);
  }

  // Set the maximum number of threads to use
  size_t num_threads = std::thread::hardware_concurrency() / 2;
  std::cout << "num_threads = " << num_threads << '\n';
  auto threads = std::vector<std::thread>(0);
  for(size_t i = 0; i < num_threads; ++i){
    threads.push_back(std::thread([&](void){
      int idx;
      while(q.pop(idx)){
	auto [l, beta] = l_beta_array[idx];
	run_simulation(l, beta);
      }
    }));
  }

  for(size_t i = 0; i < num_threads; ++i){
    if(threads[i].joinable()){
      threads[i].join();
    }
  }

}

/*! 
  \brief Solve the blackhole perturbations for a scalar field with cubic self-interaction.
*/
void run_coupled_eqn(void) {
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using namespace std::numbers;
  using std::array;
    
  const std::string dir = "output/quadratic_rsh_coupling_0001_ingoing/";
  prepare_directory_for_output(dir);

  const double r0 = 1;
  const long long int l_max = 1;  // The cutoff angular number
  const double lambda = 0.001;
  const double r_min = -600; //-400;
  const double r_max = 1200; //600;
  //const long long int N = 1 << 15;
  const long long int N = static_cast<long long int>((r_max - r_min) / 0.03);
  
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
  param.t_end = 1200; //r_max - r_source;
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
  auto observer2 = ApproximateTimeObserver(dir, {50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550.});
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

void solve_rast(void) {
  using namespace Eigen;
  using boost::multiprecision::cpp_bin_float_100;
  using boost::math::lambert_w0;
  typedef boost::multiprecision::cpp_bin_float_100 HighPrecisionScalar;

  typedef double Scalar;
  typedef Eigen::Array<Scalar, -1, 1> Vector;
  typedef Vector State;
  
  const long long int N = 10000;
  const double rast_min = -10;
  const double rast_max = 10;
  const double M = 0.5;
  const double a = 0.1;
  HighPrecisionScalar rast_min_hp = static_cast<HighPrecisionScalar>(rast_min);
  HighPrecisionScalar rast_max_hp = static_cast<HighPrecisionScalar>(rast_max);
  HighPrecisionScalar M_hp = static_cast<HighPrecisionScalar>(M);
  HighPrecisionScalar a_hp = static_cast<HighPrecisionScalar>(a);
  HighPrecisionScalar r_minus_hp = M_hp - sqrt(pow(M_hp, 2) - pow(a_hp, 2));
  HighPrecisionScalar r_plus_hp = M_hp + sqrt(pow(M_hp, 2) - pow(a_hp, 2));
  HighPrecisionScalar h_hp = (rast_max_hp - rast_min_hp) / (N - 1);
  Vector r(N + 1);
  Vector r_ast(N + 1);

  for(long long int i = 0; i < N + 1; ++i) {
    const HighPrecisionScalar r_ast_hp = rast_min_hp + i * h_hp - h_hp / HighPrecisionScalar(2);

    auto r_to_rast =
      [&](HighPrecisionScalar r_hp) -> HighPrecisionScalar {
	return -r_ast_hp + r_hp
	  + (r_plus_hp * r_plus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_plus_hp)
	  - (r_minus_hp * r_minus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_minus_hp);
      };

    std::uintmax_t it = 1000;
    boost::math::tools::eps_tolerance<HighPrecisionScalar> tol(100);
    std::pair<HighPrecisionScalar,HighPrecisionScalar> result = boost::math::tools::bisect(r_to_rast, r_plus_hp,rast_max_hp, tol, it);
    
    HighPrecisionScalar r_hp = (result.first + result.second) / 2;
    
    r[i] = r_hp.convert_to<Scalar>();
    r_ast[i] = r_ast_hp.convert_to<Scalar>();
    
    std::cout << std::setprecision(100);
    std::cout << "r = " << r[i] << ", r_ast = " << r_ast[i] << std::endl;
  }

  // HighPrecisionScalar target = static_cast<HighPrecisionScalar>(-1000);

  // auto r_to_rast =
  //   [&](HighPrecisionScalar r_hp) -> HighPrecisionScalar {
  //     std::cout << std::setprecision(100);
  //     std::cout << "input = " << r_hp << std::endl;
  //     return -target + r_hp
  //     + (r_plus_hp * r_plus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_plus_hp)
  //     - (r_minus_hp * r_minus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_minus_hp);
  //   };
  
  // std::uintmax_t it = 10000;
  // int digits = std::numeric_limits<double>::digits;
  // boost::math::tools::eps_tolerance<HighPrecisionScalar> tol(200 - 3);
  // std::pair<HighPrecisionScalar,HighPrecisionScalar> result = boost::math::tools::bisect(r_to_rast, r_plus_hp, static_cast<HighPrecisionScalar>(100), tol, it);
  // std::cout << std::setprecision(100);
  // std::cout << "solution for r_ast = -100 is " << result.first << std::endl;

  // auto r_to_rast_newton =
  //   [&](HighPrecisionScalar r_hp) -> std::pair<HighPrecisionScalar, HighPrecisionScalar> {
  //     std::cout << std::setprecision(100);
  //     std::cout << "input = " << r_hp << std::endl;
  //     return std::make_pair(-target + r_hp
  // 			    + (r_plus_hp * r_plus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_plus_hp)
  // 			    - (r_minus_hp * r_minus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_minus_hp),
  // 			    (r_hp * r_hp + a_hp * a_hp) / (r_hp * r_hp - 2 * M_hp * r_hp + a_hp * a_hp) );
  //   };
  // HighPrecisionScalar result = boost::math::tools::newton_raphson_iterate(r_to_rast_newton, 3*M_hp, r_plus_hp, static_cast<HighPrecisionScalar>(100), 80);
  
  // std::cout << std::setprecision(100);
  // std::cout << "solution for r_ast = -100 is " << result << std::endl;
}

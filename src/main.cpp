#include <vector>
#include <algorithm>
#include <functional>
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
  
  // {
  //   typedef Eigen::Array<HighPrecisionScalar, -1, 1> Vector;
  //   const HighPrecisionScalar M = HighPrecisionScalar("1.0");
  //   const HighPrecisionScalar a = HighPrecisionScalar("0.5");
  //   // const std::complex<HighPrecisionScalar> ii(HighPrecisionScalar("0"), HighPrecisionScalar("1"));
  //   const long long int grid_size = 10;
  //   auto expr = atan((a)*((pow((pow(r(seqN((2),((-4)+(grid_size)))),4))+((pow(a,2))*((r(seqN((2),((-4)+(grid_size)))))*(((2)*(M))+(r(seqN((2),((-4)+(grid_size)))))))),-1))*((sqrt(r(seqN((2),((-4)+(grid_size))))))*((sqrt((pow(a,2))+((r(seqN((2),((-4)+(grid_size)))))*(((-2)*(M))+(r(seqN((2),((-4)+(grid_size)))))))))*(sqrt((pow(r(seqN((2),((-4)+(grid_size)))),3))+((pow(a,2))*(((2)*(M))+(r(seqN((2),((-4)+(grid_size)))))))))))));

  //   Vector expr_1 = sqrt((pow(a,2))+((r(seqN((0),(grid_size))))*(((-2)*(M))+(r(seqN((0),(grid_size)))))));
  //   Vector expr_2 = sqrt((pow(r(seqN((0),(grid_size))),3))+((pow(a,2))*(((2)*(M))+(r(seqN((0),(grid_size)))))));
  //   Vector expr_3 = sqrt(r(seqN((0),(grid_size))));
  //   Vector expr_4 = atan((a)*((expr_1)*((expr_2)*((expr_3)*(pow((pow(r(seqN((0),(grid_size))),4))+((pow(a,2))*((r(seqN((0),(grid_size))))*(((2)*(M))+(r(seqN((0),(grid_size))))))),-1))))));

  // }

}

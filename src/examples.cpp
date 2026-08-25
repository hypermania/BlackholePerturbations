#include "examples.hpp"

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
#include "regge_wheeler_precise.hpp"
#include "sds_precise.hpp"
#include "teukolsky_precise.hpp"
#include "cubic_scalar.hpp"

#include "boost/type_index.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <omp.h>


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


/*! 
  \brief Solve the Teukolsky equation in Schwarzschild using high precision CPU arithmetic.
  Runs production-quality simulations with OpenMP-parallel operator evaluation.
  Saves time series at r_* = 50 (psi and Pi) and periodic snapshots.
*/
void run_teukolsky_precise_eqn(void) {
  using namespace Eigen;
  using namespace boost::numeric::odeint;
  using boost::math::lambert_w0;
  using boost::multiprecision::cpp_bin_float_100;
  using namespace std::numbers;

  typedef TeukolskyPDEPrecise Equation;
  typedef TeukolskyPDEPreciseParam Param;
  typedef Equation::Scalar Scalar;
  typedef Equation::State State;
  typedef Equation::Vector Vector;

  auto run_simulation = [](const long long int l, const long long int s, const long long int beta100, const Scalar ko_epsilon, const Scalar t_end)->void {
    // std::string format_string = "output/teukolsky_precise_far/l_%d_s_%d_beta100_%d/";
    std::string format_string = "output/teukolsky_precise/l_%d_s_%d_beta100_%d/";


    char dir_buffer[128];
    sprintf(dir_buffer, format_string.data(), static_cast<int>(l), static_cast<int>(s), static_cast<int>(beta100));
    const std::string dir(dir_buffer);
    prepare_directory_for_output(dir);  

    const Scalar beta = beta100 / Scalar(100);
    const Scalar M = Scalar(1) / Scalar(2);

    const Scalar r_min = -500;
    const Scalar r_max =  1000;
    // const Scalar r_min = 0;
    // const Scalar r_max =  1500;
    const long long int N = static_cast<long long int>((r_max - r_min) / 0.03);

    const Scalar t_start = 0;
    const Scalar delta_t = 0.01;

    Param param;
    param.s = s;
    param.l = l;
    param.M = M;
    param.r_min = r_min;
    param.r_max = r_max;
    param.N = N;
    param.ko_epsilon = ko_epsilon;
    param.t_start = t_start;
    param.t_end = t_end;
    param.t_interval = Scalar(1) / Scalar(2);
    param.delta_t = delta_t;

    save_param_for_Mathematica(param, dir);
  
    Equation eqn(param);
  
    // auto stepper = runge_kutta_fehlberg78<State, Scalar, State, Scalar>();
    auto stepper = runge_kutta_dopri5<State, Scalar, State, Scalar>();
    //auto stepper = runge_kutta4_classic<State, Scalar, State, Scalar>();
    
    const long long int rIdx = r_ast_to_i(param.r_min.convert_to<double>(), param.r_max.convert_to<double>(), N, 50.0);
    auto observer1 = FixedPositionObserver(dir, {rIdx, rIdx + (N+1)});

    // Snapshots every 50 time units
    std::vector<double> snap_times;
    int n_snaps = static_cast<int>(static_cast<double>(t_end) / 50.0 + 1);
    for(int i = 0; i < n_snaps; ++i) snap_times.push_back(50.0 * i);
    auto observer2 = ApproximateTimeObserver(dir, snap_times);
    auto observer = ObserverPack(observer1, observer2);

    // Outgoing Gaussian source at r_* = 10
    Vector r_ast = eqn.compute_r_ast_vector(r_min, r_max, N);
      
    const Scalar r_source = Scalar(10);
    // const Scalar r_source = Scalar(500);
    const Scalar sigma = Scalar(1) / Scalar(2);
    const Scalar pf = pow(Scalar(2 * pi), Scalar(-0.5)) / sigma;
    const Scalar denom = Scalar(2) * sigma * sigma;
    
    Vector r = eqn.compute_r_vector(r_min, r_max, N, Scalar(2) * M);
    Vector r_beta = r;
    r_beta = r_beta.pow(-beta);
    r_beta.head(r_ast_to_i(param.r_min.convert_to<double>(), param.r_max.convert_to<double>(), N, r_source.convert_to<double>())) = 0;

    eqn.Q = [&](const Scalar &t, Vector &result)->void {
#pragma omp parallel for schedule(static)
      for(long long int i = 0; i <= N; ++i) {
        Scalar arg = t - r_ast[i] + r_source;
        result[i] = pf * boost::multiprecision::exp(-arg * arg / denom) * r_beta[i];
      }
    };
    
    Vector state = Vector::Zero(2 * (N+1));
      
    run_and_measure_time("Solving Teukolsky precise equation",
			 [&](){
			   int num_steps = integrate_const(stepper, std::ref(eqn), state, t_start, t_end, delta_t, std::ref(observer));
			   std::cout << "total number of steps = " << num_steps << '\n';
			 } );
    observer.save();
  };

  // run_simulation(1, -1, 190, Scalar("0.7"), Scalar(1000));
  // run_simulation(1, -1, 180, Scalar("0.7"), Scalar(1000));
  // run_simulation(1, -1, 170, Scalar("0.7"), Scalar(1000));
  // run_simulation(1, -1, 160, Scalar("0.7"), Scalar(1000));
  // run_simulation(1, -1, 150, Scalar("0.7"), Scalar(1000));

  // run_simulation(2, -1, 200, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 199, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 198, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 197, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 196, Scalar("0.7"), Scalar(1000));

  // run_simulation(2, -1, 300, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 280, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 260, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 240, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 230, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 220, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 210, Scalar("0.7"), Scalar(1000));

  // run_simulation(2, -1, 200, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 180, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 220, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 240, Scalar("0.7"), Scalar(1000));
  
  // run_simulation(2, -2, 195, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -2, 199, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -2, 220, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -2, 240, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -2, 260, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -2, 280, Scalar("0.7"), Scalar(1000));

  run_simulation(1, -1, 400, Scalar("0.7"), Scalar(1000));
  run_simulation(2, -2, 400, Scalar("0.7"), Scalar(1000));

  
  // run_simulation(2, -1, 195, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 190, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 180, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 170, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 160, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -1, 150, Scalar("0.7"), Scalar(1000));
  return;
  
  // Low-ℓ tail runs
  // run_simulation(1, -1, 2, Scalar("0.5"), Scalar(1000));
  // run_simulation(1, -1, 3, Scalar("0.5"), Scalar(1000));
  // run_simulation(1, -1, 0, Scalar("0.5"), Scalar(1000));
  // run_simulation(1, -1, 1, Scalar("0.5"), Scalar(1000));

  // run_simulation(2, -1, 0, Scalar("0.5"), Scalar(1000));
  // run_simulation(2, -1, 1, Scalar("0.5"), Scalar(1000));
  // run_simulation(2, -1, 2, Scalar("0.5"), Scalar(1000));
  // run_simulation(2, -1, 3, Scalar("0.5"), Scalar(1000));

  //run_simulation(2, -2, 0, Scalar("0.7"), Scalar(1000));
  //run_simulation(2, -2, 1, Scalar("0.7"), Scalar(1000));
  //run_simulation(2, -2, 2, Scalar("0.7"), Scalar(1000));
  // run_simulation(2, -2, 3, Scalar("0.7"), Scalar(1000));

  // run_simulation(3, -2, 2, Scalar("0.7"), Scalar(1000));
    
  
  // run_simulation(1, -1, Scalar("0.5"), Scalar(500));
  // run_simulation(3, -1, Scalar("0.5"), Scalar(500));
  // run_simulation(2, -2, Scalar("0.7"), Scalar(500));
  // run_simulation(3, -2, Scalar("0.7"), Scalar(500));


  auto run_dirac_simulation = [](const long long int l, const long long int s, const Scalar ko_epsilon, const Scalar t_end)->void {
    std::string format_string = "output/teukolsky_precise_dirac/l_%d_s_%d/";
      
    char dir_buffer[128];
    sprintf(dir_buffer, format_string.data(), static_cast<int>(l), static_cast<int>(s));
    const std::string dir(dir_buffer);
    prepare_directory_for_output(dir);  
  
    const Scalar M = Scalar(1) / Scalar(2);

    const Scalar r_min = -500;
    const Scalar r_max =  1000;
    const long long int N = static_cast<long long int>((r_max - r_min) / 0.03);

    const Scalar t_start = 0;
    const Scalar delta_t = 0.01;

    Param param;
    param.s = s;
    param.l = l;
    param.M = M;
    param.r_min = r_min;
    param.r_max = r_max;
    param.N = N;
    param.ko_epsilon = ko_epsilon;
    param.t_start = t_start;
    param.t_end = t_end;
    param.t_interval = Scalar(1) / Scalar(2);
    param.delta_t = delta_t;

    save_param_for_Mathematica(param, dir);
  
    Equation eqn(param);
  
    // auto stepper = runge_kutta_fehlberg78<State, Scalar, State, Scalar>();
    auto stepper = runge_kutta_dopri5<State, Scalar, State, Scalar>();
    //auto stepper = runge_kutta4_classic<State, Scalar, State, Scalar>();
    
    const long long int rIdx = r_ast_to_i(param.r_min.convert_to<double>(), param.r_max.convert_to<double>(), N, 50.0);
    auto observer1 = FixedPositionObserver(dir, {rIdx, rIdx + (N+1)});

    // Snapshots every 50 time units
    std::vector<double> snap_times;
    int n_snaps = static_cast<int>(static_cast<double>(t_end) / 50.0 + 1);
    for(int i = 0; i < n_snaps; ++i) snap_times.push_back(50.0 * i);
    auto observer2 = ApproximateTimeObserver(dir, snap_times);
    auto observer = ObserverPack(observer1, observer2);

    // Outgoing Gaussian source at r_* = 10
    Vector r_ast = eqn.compute_r_ast_vector(r_min, r_max, N);
      
    const Scalar r_source = Scalar(10);
    const Scalar sigma = Scalar(1) / Scalar(2);
    const Scalar pf = pow(Scalar(2 * pi), Scalar(-0.5)) / sigma;
    const Scalar denom = Scalar(2) * sigma * sigma;

    Vector spatial_source(N + 1);
#pragma omp parallel for schedule(static)
    for(long long int i = 0; i <= N; ++i) {
      const Scalar radial_offset = r_ast[i] - r_source;
      spatial_source[i] = pf * boost::multiprecision::exp(
          -(radial_offset * radial_offset) / denom);
    }

    eqn.set_separable_source(std::move(spatial_source),
                             [r_source, denom](const Scalar &t)->Scalar {
      const Scalar time_offset = t - r_source;
      return boost::multiprecision::exp(
          -(time_offset * time_offset) / denom);
    });
      
    Vector state = Vector::Zero(2 * (N+1));
      
    run_and_measure_time("Solving Teukolsky precise equation",
			 [&](){
			   int num_steps = integrate_const(stepper, std::ref(eqn), state, t_start, t_end, delta_t, std::ref(observer));
			   std::cout << "total number of steps = " << num_steps << '\n';
			 } );
    observer.save();
  };

  run_dirac_simulation(1, -1, Scalar("0.7"), Scalar(1000));
  run_dirac_simulation(2, -1, Scalar("0.7"), Scalar(1000));
  run_dirac_simulation(2, -2, Scalar("0.7"), Scalar(1000));
  run_dirac_simulation(3, -2, Scalar("0.7"), Scalar(1000));
  
  run_dirac_simulation(3, -1, Scalar("0.7"), Scalar(1000));
  // run_dirac_simulation(0, 0, Scalar("0.7"), Scalar(1000));


}


/*!
  \brief Solve one sourced Schwarzschild-de Sitter master equation.

  The configuration block below selects the spin, multipole, cosmological
  constant, source profile, and waveform. Evolution and source evaluation use
  binary128; observers intentionally save double arrays, matching the existing
  precise runners.
*/
void run_sds_precise_eqn(void) {
  using namespace boost::numeric::odeint;

  typedef SdSMasterPDEPrecise Equation;
  typedef SdSMasterPDEPreciseParam Param;
  typedef Equation::Scalar Scalar;
  typedef Equation::State State;

  // Production configuration.
  const long long int s = 0;
  const long long int l = 1;
  const Scalar M("0.5");
  const Scalar Lambda("1e-4");
  const Scalar beta(2);
  const SdSSourceProfile profile = SdSSourceProfile::TortoisePower;
  const SdSWaveform waveform = SdSWaveform::Gaussian;

  const Scalar r_min(-500);
  const Scalar r_max(1000);
  const long long int N = static_cast<long long int>(
      ((r_max - r_min) / Scalar("0.03")).convert_to<long long int>());
  const Scalar t_start(0);
  const Scalar t_end(1000);
  const Scalar delta_t("0.01");

  const std::string dir = std::string("output/sds_precise/")
      + sds_source_profile_name(profile) + "/";
  prepare_directory_for_output(dir);

  Param param;
  param.s = s;
  param.l = l;
  param.M = M;
  param.Lambda = Lambda;
  param.r_min = r_min;
  param.r_max = r_max;
  param.N = N;
  param.t_start = t_start;
  param.t_end = t_end;
  param.t_interval = Scalar("0.5");
  param.delta_t = delta_t;
  save_param_for_Mathematica(param, dir);

  SdSTranslatedSourceParam source;
  source.profile = profile;
  source.waveform = waveform;
  source.beta = beta;
  source.amplitude = 1;
  source.u_center = -10;
  source.sigma = Scalar("0.5");
  source.onset_time = t_start;
  source.cutoff_sigma = 12;
  source.L = 1;
  source.x0 = 100;
  source.X0 = 0;
  source.X1 = 20;
  Equation equation(param, SdSSource(source));

  const double observer_x = 50.0;
  const long long int observer_index = r_ast_to_i(
      param.r_min.convert_to<double>(), param.r_max.convert_to<double>(), N,
      observer_x);
  auto fixed_observer = FixedPositionObserver(
      dir, {observer_index, observer_index + N + 1});

  std::vector<double> snapshot_times;
  const int snapshot_count = static_cast<int>(t_end.convert_to<double>() / 50.0)
                             + 1;
  for(int i = 0; i < snapshot_count; ++i) {
    snapshot_times.push_back(50.0 * i);
  }
  auto snapshot_observer = ApproximateTimeObserver(dir, snapshot_times);
  auto observer = ObserverPack(fixed_observer, snapshot_observer);

  const Scalar effective_u_min = source.u_center
                                 - source.cutoff_sigma * source.sigma;
  const Scalar finite_domain_limit = Scalar(2) * r_max - Scalar(observer_x)
                                     + effective_u_min;
  {
    std::ofstream metadata(dir + "source_and_geometry.txt");
    metadata << std::setprecision(36)
             << "profile " << sds_source_profile_name(profile) << '\n'
             << "waveform " << sds_waveform_name(waveform) << '\n'
             << "beta " << beta << '\n'
             << "amplitude " << source.amplitude << '\n'
             << "u_center " << source.u_center << '\n'
             << "sigma " << source.sigma << '\n'
             << "onset_time " << source.onset_time << '\n'
             << "cutoff_sigma " << source.cutoff_sigma << '\n'
             << "L " << source.L << '\n'
             << "x0 " << source.x0 << '\n'
             << "X0 " << source.X0 << '\n'
             << "X1 " << source.X1 << '\n'
             << "r_black_hole " << equation.r_black_hole << '\n'
             << "r_cosmological " << equation.r_cosmological << '\n'
             << "r_negative " << equation.r_negative << '\n'
             << "kappa_black_hole " << equation.kappa_black_hole << '\n'
             << "kappa_cosmological " << equation.kappa_cosmological << '\n'
             << "tortoise_convention x(3M)=3M+2M*log(1/2)\n"
             << "observer_x " << observer_x << '\n'
             << "effective_u_min " << effective_u_min << '\n'
             << "finite_domain_fit_limit " << finite_domain_limit << '\n';
  }

  if(t_end >= finite_domain_limit) {
    std::cerr << "Warning: t_end reaches the conservative finite-domain "
                 "contamination bound t < " << finite_domain_limit << '\n';
  }

  State state = State::Zero(2 * (N + 1));
  auto stepper = runge_kutta_dopri5<State, Scalar, State, Scalar>();
  run_and_measure_time("Solving precise sourced SdS equation", [&]() {
    const int steps = integrate_const(
        stepper, std::ref(equation), state, t_start, t_end, delta_t,
        std::ref(observer));
    std::cout << "total number of steps = " << steps << '\n';
  });
  fixed_observer.save();
  snapshot_observer.save();
}


namespace {

using SdSScalar = SdSMasterPDEPrecise::Scalar;

long long int nearest_sds_grid_index(const SdSMasterPDEPrecise &equation,
                                     const double requested_x) {
  const SdSScalar continuous_index =
      (SdSScalar(requested_x) - equation.param.r_min) * equation.inv_h
      + SdSScalar("0.5");
  const long long int index = static_cast<long long int>(
      std::llround(continuous_index.convert_to<double>()));
  return std::clamp(index, 0LL, equation.grid_size - 1);
}


} // namespace


/*!
  \brief Run one member of the sourced SdS areal-radius scan.

  Here q = 9 Lambda M^2. The scan fixes M = 0.5 and uses the
  normalized nonzero-mean Gaussian source

    F(u) = A exp[-(u-u0)^2/(2 sigma^2)] / (sqrt(2 pi) sigma),  A = 1,

  multiplied by r^{-beta}. Every raw and derived artifact for one parameter
  set is stored in one output directory.
*/
void run_sds_areal_scan(const std::string &q_code, const long long int s,
                        const long long int l,
                        const std::string &beta_text) {
  using namespace boost::numeric::odeint;

  omp_set_dynamic(0);

  using Equation = SdSMasterPDEPrecise;
  using Param = SdSMasterPDEPreciseParam;
  using Scalar = Equation::Scalar;
  using State = Equation::State;

  const Scalar q(sds_q_code_to_decimal(q_code));
  const Scalar beta(beta_text);

  const Scalar M("0.5");
  const Scalar Lambda = q / (Scalar(9) * M * M);
  const Scalar x_min(-500);
  const Scalar x_max(1000);
  const Scalar nominal_dx("0.03");
  const long long int N = static_cast<long long int>(
      ((x_max - x_min) / nominal_dx).convert_to<long long int>());
  const Scalar t_start(0);
  const Scalar t_end(1000);
  const Scalar delta_t("0.01");

  Param param;
  param.s = s;
  param.l = l;
  param.M = M;
  param.Lambda = Lambda;
  param.r_min = x_min;
  param.r_max = x_max;
  param.N = N;
  param.t_start = t_start;
  param.t_end = t_end;
  param.t_interval = delta_t;
  param.delta_t = delta_t;
  SdSTranslatedSourceParam source;
  source.profile = SdSSourceProfile::ArealPower;
  source.waveform = SdSWaveform::Gaussian;
  source.beta = beta;
  source.amplitude = 1;
  source.u_center = -10;
  source.sigma = Scalar("0.5");
  source.onset_time = t_start;
  source.cutoff_sigma = 12;
  Equation equation(param, SdSSource(source));

  std::ostringstream beta_label;
  beta_label << std::setprecision(std::numeric_limits<Scalar>::max_digits10)
             << beta;
  std::ostringstream directory;
  directory << "output/sds_areal_scan/q_" << q_code
            << "_l_" << l << "_beta_" << beta_label.str() << "/";
  const std::string dir = directory.str();
  prepare_directory_for_output(dir);
  save_param_for_Mathematica(param, dir);

  constexpr double requested_observer = 50;
  const long long int observer_index = nearest_sds_grid_index(
      equation, requested_observer);
  const std::vector<long long int> observed_components = {
      observer_index, observer_index + equation.grid_size};

  auto fixed_observer = FixedPositionObserver(dir, observed_components);
  const std::vector<double> snapshot_times = {
      0, 10, 20, 30, 40, 50, 75, 100, 125, 150,
      200, 250, 300, 350, 400, 450, 500, 550, 600,
      650, 700, 750, 800, 850, 900, 950, 1000
  };
  auto snapshot_observer = ApproximateTimeObserver(dir, snapshot_times);
  auto observer = ObserverPack(fixed_observer, snapshot_observer);

  Eigen::ArrayXd x_grid(equation.grid_size);
  for(long long int i = 0; i < equation.grid_size; ++i) {
    x_grid[i] = equation.grid_coordinate(i).convert_to<double>();
  }
  const Eigen::ArrayXd r_grid = equation.r.cast<double>();
  const Eigen::ArrayXd f_grid = equation.f.cast<double>();
  const Eigen::ArrayXd potential_grid = equation.V.cast<double>();
  write_to_file(x_grid, dir + "x_grid.dat");
  write_to_file(r_grid, dir + "r_grid.dat");
  write_to_file(f_grid, dir + "f_grid.dat");
  write_to_file(potential_grid, dir + "potential_grid.dat");
  write_to_file(snapshot_times, dir + "snapshot_times_requested.dat");

  const Scalar effective_u_min = source.u_center
                                 - source.cutoff_sigma * source.sigma;
  const Scalar finite_domain_limit = Scalar(2) * x_max
                                     - Scalar(requested_observer)
                                     + effective_u_min;
  {
    std::ofstream metadata(dir + "metadata.txt");
    metadata << std::setprecision(36)
             << "run_name q_" << q_code << "_l_" << l
             << "_beta_" << beta_label.str() << '\n'
             << "q_9LambdaM2 " << q << '\n'
             << "s " << s << '\n'
             << "l " << l << '\n'
             << "beta " << beta << '\n'
             << "M " << M << '\n'
             << "Lambda " << Lambda << '\n'
             << "x_min " << x_min << '\n'
             << "x_max " << x_max << '\n'
             << "N " << N << '\n'
             << "grid_size " << equation.grid_size << '\n'
             << "dx " << equation.grid_space() << '\n'
             << "t_start " << t_start << '\n'
             << "t_end " << t_end << '\n'
             << "delta_t " << delta_t << '\n'
             << "openmp_threads " << omp_get_max_threads() << '\n'
             << "source_profile "
             << sds_source_profile_name(source.profile) << '\n'
             << "source_waveform " << sds_waveform_name(source.waveform)
             << '\n'
             << "source_amplitude_integral_A " << source.amplitude << '\n'
             << "source_u_center " << source.u_center << '\n'
             << "source_sigma " << source.sigma << '\n'
             << "source_onset_time " << source.onset_time << '\n'
             << "source_cutoff_sigma " << source.cutoff_sigma << '\n'
             << "r_black_hole " << equation.r_black_hole << '\n'
             << "r_cosmological " << equation.r_cosmological << '\n'
             << "r_negative " << equation.r_negative << '\n'
             << "kappa_black_hole " << equation.kappa_black_hole << '\n'
             << "kappa_cosmological " << equation.kappa_cosmological << '\n'
             << "tortoise_convention x(3M)=3M+2M*log(1/2)\n"
             << "finite_domain_fit_limit " << finite_domain_limit << '\n'
             << "time_series_layout row_major_[psi_x50,Pi_x50]\n"
             << "observer_requested_x " << requested_observer << '\n'
             << "observer_index " << observer_index << '\n'
             << "observer_actual_x "
             << equation.grid_coordinate(observer_index) << '\n';
  }

  State state = State::Zero(2 * equation.grid_size);
  auto stepper = runge_kutta_dopri5<State, Scalar, State, Scalar>();
  const auto wall_start = std::chrono::steady_clock::now();
  const int steps = integrate_const(
      stepper, std::ref(equation), state, t_start, t_end, delta_t,
      std::ref(observer));
  const double wall_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - wall_start).count();

  fixed_observer.save();
  snapshot_observer.save();
  const Eigen::ArrayXd final_state = state.cast<double>();
  write_to_file(final_state, dir + "final_state.dat");
  {
    std::ofstream summary(dir + "run_summary.txt");
    summary << std::setprecision(17)
            << "steps " << steps << '\n'
            << "time_samples " << fixed_observer.t_list.size() << '\n'
            << "snapshots_saved " << snapshot_observer.t_list.size() << '\n'
            << "wall_seconds " << wall_seconds << '\n';
  }
  {
    std::ofstream complete(dir + "COMPLETE");
    complete << "Simulation output complete\n";
  }

  std::cout << "SdS areal scan run complete: " << dir << '\n'
            << "steps = " << steps << ", wall time = " << wall_seconds
            << " s\n";
}

#include "teukolsky_scalar.hpp"
#include <map>
#include <iostream>
#include <thread>
#include <boost/lockfree/queue.hpp>

#include "teukolsky_scalar.gen"

// typedef std::function<void(TeukolskyScalarPDE::Vector &)> CoeffFunc;

// CoeffFunc f_11 = [](TeukolskyScalarPDE::Vector &r) {
//   r = 0;
//  };

// TeukolskyScalarPDE::Vector f_2(TeukolskyScalarPDE::Vector &r) {
//   r = 0;
//   return TeukolskyScalarPDE::Vector(10);
// }

// typedef decltype(&f_2) CoeffFunc2;


// std::unordered_map<long long int, std::unordered_map<long long int, CoeffFunc>> my_map;
// std::map<long long int, CoeffFunc2> my_map_2;

// void init_func_table(void) {
//   my_map[0];
//   my_map[0][0] = f_11;
//   my_map_2.emplace(0, f_2);
// }

/*!
  \brief Given info for all lm -> (lm1, coeff_idx), use cutoff l_max, keep only terms with lm and lm1 below cutoff.
 */
TeukolskyScalarPDE::CouplingInfo make_coupling_info_map(const std::vector<std::vector<std::pair<long long int, long long int>>> &info, const long long int l_max) {
  const long long int lm_size = (l_max + 1) * (l_max + 1);
  TeukolskyScalarPDE::CouplingInfo result;
  result.resize(lm_size);
  for(int idx = 0; idx < lm_size; ++idx){
    //result[idx];
    for(auto p : info[idx]){
      auto [idx1, coeff1] = p;
      if(idx1 < lm_size){
	result[idx][idx1] = coeff1;
      }
    }
  }
  return result;
}

TeukolskyScalarPDE::TeukolskyScalarPDE(Param param_) : param(param_) {
  using boost::multiprecision::cpp_bin_float_100;
  using namespace Eigen;
  // using boost::math::lambert_w0;
    
  const Scalar rast_min = param.rast_min;
  const Scalar rast_max = param.rast_max;
  const auto N = param.N;
  const Scalar M = param.M;
  const Scalar a = param.a;
  const auto s = param.s;
  const auto l_max = param.l_max;

  grid_size = N + 1;
  lm_size = (l_max + 1) * (l_max + 1);
    
  const Scalar h = (rast_max - rast_min) / (N - 1);

  //Vector r = compute_r_vector(rast_min, rast_max, N, M, a);
  
  psi_lm_map = make_coupling_info_map(psi_lm_coupling_info, l_max);
  dr_psi_lm_map = make_coupling_info_map(dr_psi_lm_coupling_info, l_max);
  drdr_psi_lm_map = make_coupling_info_map(drdr_psi_lm_coupling_info, l_max);
  dt_psi_lm_map = make_coupling_info_map(dt_psi_lm_coupling_info, l_max);

  auto a_hp = static_cast<HighPrecisionScalar>(a);
  auto M_hp = static_cast<HighPrecisionScalar>(M);
  auto rast_min_hp = static_cast<HighPrecisionScalar>(rast_min);
  auto rast_max_hp = static_cast<HighPrecisionScalar>(rast_max);

  auto t1 = std::chrono::system_clock::now();
  
  auto r_hp = compute_hp_r_vector(rast_min_hp, rast_max_hp, N, M_hp, a_hp);
  
  auto t2 = std::chrono::system_clock::now();
  
  coeffs = compute_coeffs(a_hp, M_hp, r_hp);

  auto t3 = std::chrono::system_clock::now();
  std::chrono::duration<double> time_diff_1 = t2 - t1;
  std::chrono::duration<double> time_diff_2 = t3 - t2;
  std::cout << std::setw(9) << "time spent 1 = " << time_diff_1.count() << " s" << '\n';
  std::cout << std::setw(9) << "time spent 2 = " << time_diff_2.count() << " s" << '\n';
  
  Q = [N](const Scalar t)->Vector{ return Vector::Zero(N+1); };
}


TeukolskyScalarPDE::Vector TeukolskyScalarPDE::compute_r_vector(const Scalar rast_min, const Scalar rast_max, const long long int N, const Scalar M, const Scalar a) {
  return compute_hp_r_vector(static_cast<HighPrecisionScalar>(rast_min),
			     static_cast<HighPrecisionScalar>(rast_max),
			     N,
			     static_cast<HighPrecisionScalar>(M),
			     static_cast<HighPrecisionScalar>(a)).cast<double>();
}


TeukolskyScalarPDE::HighPrecisionVector TeukolskyScalarPDE::compute_hp_r_vector(const TeukolskyScalarPDE::HighPrecisionScalar rast_min, const TeukolskyScalarPDE::HighPrecisionScalar rast_max, const long long int N, const TeukolskyScalarPDE::HighPrecisionScalar M, const TeukolskyScalarPDE::HighPrecisionScalar a) {
  using boost::multiprecision::cpp_bin_float_100;
  using boost::math::lambert_w0;
  typedef boost::multiprecision::cpp_bin_float_100 HighPrecisionScalar;
  HighPrecisionScalar rast_min_hp = static_cast<HighPrecisionScalar>(rast_min);
  HighPrecisionScalar rast_max_hp = static_cast<HighPrecisionScalar>(rast_max);
  HighPrecisionScalar M_hp = static_cast<HighPrecisionScalar>(M);
  HighPrecisionScalar a_hp = static_cast<HighPrecisionScalar>(a);
  HighPrecisionScalar r_minus_hp = M_hp - sqrt(pow(M_hp, 2) - pow(a_hp, 2));
  HighPrecisionScalar r_plus_hp = M_hp + sqrt(pow(M_hp, 2) - pow(a_hp, 2));
  HighPrecisionScalar h_hp = (rast_max_hp - rast_min_hp) / (N - 1);
  TeukolskyScalarPDE::HighPrecisionVector r(N + 1);
  for(long long int i = 0; i < N + 1; ++i) {
    const HighPrecisionScalar r_ast_hp = rast_min_hp + i * h_hp - h_hp / HighPrecisionScalar(2);

    auto r_to_rast =
      [&](HighPrecisionScalar r_hp) -> HighPrecisionScalar {
	return - r_ast_hp + r_hp
	  + (r_plus_hp * r_plus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_plus_hp)
	  - (r_minus_hp * r_minus_hp + a_hp * a_hp) / (r_plus_hp - r_minus_hp) * log(r_hp - r_minus_hp);
      };

    std::uintmax_t it = 1000;
    boost::math::tools::eps_tolerance<HighPrecisionScalar> tol(100);
    std::pair<HighPrecisionScalar,HighPrecisionScalar> result = boost::math::tools::bisect(r_to_rast, r_plus_hp, rast_max_hp, tol, it);
    
    HighPrecisionScalar r_hp = (result.first + result.second) / 2;
    
    r[i] = r_hp;
  }
  return r;
}


std::pair<long long int, long long int> TeukolskyScalarPDE::idx_to_lm(const long long int idx) {
  const long long int l = static_cast<long long int>(floor(sqrt(idx)));
  const long long int m = idx - l * l - l;
  return std::make_pair(l, m);
}

void TeukolskyScalarPDE::compute_derivatives(const State &x, ComplexVector &dr_psi_lm, ComplexVector &drdr_psi_lm) const {
  using namespace Eigen;
  const auto N = param.N;
  const Scalar r_min = param.rast_min;
  const Scalar r_max = param.rast_max;
  const Scalar h = (r_max - r_min) / (N - 1);
  for(long long int lm = 0; lm < lm_size; ++lm){
dr_psi_lm(lm*grid_size+0) = ((-2.0833333333333333333)*((pow(h,-1))*(x((lm)*grid_size+(0)))))+(((4.0000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+(1)))))+(((-3.0000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+(2)))))+(((1.3333333333333333333)*((pow(h,-1))*(x((lm)*grid_size+(3)))))+((-0.25000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+(4))))))));
dr_psi_lm(lm*grid_size+1) = ((-0.25000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+(0)))))+(((-0.83333333333333333333)*((pow(h,-1))*(x((lm)*grid_size+(1)))))+(((1.5000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+(2)))))+(((-0.50000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+(3)))))+((0.083333333333333333333)*((pow(h,-1))*(x((lm)*grid_size+(4))))))));
dr_psi_lm(seqN(lm*grid_size+2,grid_size-4)) = ((0.083333333333333333333)*((pow(h,-1))*(x(seqN((lm)*grid_size+(0),((-4)+(grid_size)))))))+(((-0.66666666666666666667)*((pow(h,-1))*(x(seqN((lm)*grid_size+(1),((-4)+(grid_size)))))))+(((0.66666666666666666667)*((pow(h,-1))*(x(seqN((lm)*grid_size+(3),((-4)+(grid_size)))))))+((-0.083333333333333333333)*((pow(h,-1))*(x(seqN((lm)*grid_size+(4),((-4)+(grid_size)))))))));
dr_psi_lm(lm*grid_size+grid_size-2) = ((-0.083333333333333333333)*((pow(h,-1))*(x((lm)*grid_size+((-5)+(grid_size))))))+(((0.50000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+((-4)+(grid_size))))))+(((-1.5000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+((-3)+(grid_size))))))+(((0.83333333333333333333)*((pow(h,-1))*(x((lm)*grid_size+((-2)+(grid_size))))))+((0.25000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+((-1)+(grid_size)))))))));
dr_psi_lm(lm*grid_size+grid_size-1) = ((0.25000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+((-5)+(grid_size))))))+(((-1.3333333333333333333)*((pow(h,-1))*(x((lm)*grid_size+((-4)+(grid_size))))))+(((3.0000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+((-3)+(grid_size))))))+(((-4.0000000000000000000)*((pow(h,-1))*(x((lm)*grid_size+((-2)+(grid_size))))))+((2.0833333333333333333)*((pow(h,-1))*(x((lm)*grid_size+((-1)+(grid_size)))))))));
drdr_psi_lm(lm*grid_size+0) = ((2.9166666666666666667)*((pow(h,-2))*(x((lm)*grid_size+(0)))))+(((-8.6666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+(1)))))+(((9.5000000000000000000)*((pow(h,-2))*(x((lm)*grid_size+(2)))))+(((-4.6666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+(3)))))+((0.91666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+(4))))))));
drdr_psi_lm(lm*grid_size+1) = ((0.91666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+(0)))))+(((-1.6666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+(1)))))+(((0.50000000000000000000)*((pow(h,-2))*(x((lm)*grid_size+(2)))))+(((0.33333333333333333333)*((pow(h,-2))*(x((lm)*grid_size+(3)))))+((-0.083333333333333333333)*((pow(h,-2))*(x((lm)*grid_size+(4))))))));
drdr_psi_lm(seqN(lm*grid_size+2,grid_size-4)) = ((-0.083333333333333333333)*((pow(h,-2))*(x(seqN((lm)*grid_size+(0),((-4)+(grid_size)))))))+(((1.3333333333333333333)*((pow(h,-2))*(x(seqN((lm)*grid_size+(1),((-4)+(grid_size)))))))+(((-2.5000000000000000000)*((pow(h,-2))*(x(seqN((lm)*grid_size+(2),((-4)+(grid_size)))))))+(((1.3333333333333333333)*((pow(h,-2))*(x(seqN((lm)*grid_size+(3),((-4)+(grid_size)))))))+((-0.083333333333333333333)*((pow(h,-2))*(x(seqN((lm)*grid_size+(4),((-4)+(grid_size))))))))));
drdr_psi_lm(lm*grid_size+grid_size-2) = ((-0.083333333333333333333)*((pow(h,-2))*(x((lm)*grid_size+((-5)+(grid_size))))))+(((0.33333333333333333333)*((pow(h,-2))*(x((lm)*grid_size+((-4)+(grid_size))))))+(((0.50000000000000000000)*((pow(h,-2))*(x((lm)*grid_size+((-3)+(grid_size))))))+(((-1.6666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+((-2)+(grid_size))))))+((0.91666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+((-1)+(grid_size)))))))));
drdr_psi_lm(lm*grid_size+grid_size-1) = ((0.91666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+((-5)+(grid_size))))))+(((-4.6666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+((-4)+(grid_size))))))+(((9.5000000000000000000)*((pow(h,-2))*(x((lm)*grid_size+((-3)+(grid_size))))))+(((-8.6666666666666666667)*((pow(h,-2))*(x((lm)*grid_size+((-2)+(grid_size))))))+((2.9166666666666666667)*((pow(h,-2))*(x((lm)*grid_size+((-1)+(grid_size)))))))));
    
  }
}
  

void TeukolskyScalarPDE::operator()(const State &x, State &dxdt, const Scalar t) {
  using namespace Eigen;
  const auto N = param.N;
  const long long int dt_grid_begin = lm_size * grid_size;

  dxdt(seqN(0, dt_grid_begin)) = x(seqN(dt_grid_begin, dt_grid_begin));
  
  ComplexVector dr_psi_lm(lm_size * grid_size);
  ComplexVector drdr_psi_lm(lm_size * grid_size);

  compute_derivatives(x, dr_psi_lm, drdr_psi_lm);

  // For each spherical harmonic, set the second order time derivative
  for(int lm = 0; lm < lm_size; ++lm){
    auto dtdt_psi_lm = dxdt(seqN(dt_grid_begin + lm * grid_size, grid_size));
    dtdt_psi_lm = 0;
    for(auto [lm1, idx1] : psi_lm_map[lm]){
      dtdt_psi_lm -= coeffs[idx1] * x(seqN(lm1 * grid_size, grid_size));
    }
    for(auto [lm1, idx1] : dt_psi_lm_map[lm]){
      dtdt_psi_lm -= coeffs[idx1] * x(seqN(dt_grid_begin + lm1 * grid_size, grid_size));
    }
    for(auto [lm1, idx1] : dr_psi_lm_map[lm]){
      dtdt_psi_lm -= coeffs[idx1] * dr_psi_lm(seqN(lm1 * grid_size, grid_size));
    }
    for(auto [lm1, idx1] : drdr_psi_lm_map[lm]){
      dtdt_psi_lm -= coeffs[idx1] * drdr_psi_lm(seqN(lm1 * grid_size, grid_size));
    }
  }
}

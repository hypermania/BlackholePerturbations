#include "teukolsky_scalar.hpp"
#include <map>
#include <iostream>

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

TeukolskyScalarPDE::CouplingInfo make_coupling_info_map(const std::vector<std::vector<std::pair<long long int, long long int>>> &info, const long long int l_max) {
  const long long int lm_size = (l_max + 1) * (l_max + 1);
  TeukolskyScalarPDE::CouplingInfo result;
  result.resize(lm_size);
  for(int idx = 0; idx < lm_size; ++idx){
    result[idx];
    for(auto p : info[idx]){
      auto [idx1, coeff1] = p;
      result[idx][idx1] = coeff1;
    }
  }
  return result;
}

TeukolskyScalarPDE::TeukolskyScalarPDE(Param param_) : param(param_) {
  using boost::multiprecision::cpp_bin_float_100;
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

  Vector r = compute_r_vector(rast_min, rast_max, N, M, a);
  
  psi_lm_map = make_coupling_info_map(psi_lm_coupling_info, l_max);
  dr_psi_lm_map = make_coupling_info_map(dr_psi_lm_coupling_info, l_max);
  drdr_psi_lm_map = make_coupling_info_map(drdr_psi_lm_coupling_info, l_max);
  dt_psi_lm_map = make_coupling_info_map(dt_psi_lm_coupling_info, l_max);

  std::cout << "point 0\n";
  
  auto hp_r = compute_hp_r_vector(static_cast<HighPrecisionScalar>(rast_min),
				  static_cast<HighPrecisionScalar>(rast_max),
				  N,
				  static_cast<HighPrecisionScalar>(M),
				  static_cast<HighPrecisionScalar>(a)
				  );

  std::cout << "point 1\n";
  coeffs = compute_coeffs(static_cast<HighPrecisionScalar>(a),
			  static_cast<HighPrecisionScalar>(M),
			  hp_r);
    std::cout << "point 2\n";
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

void TeukolskyScalarPDE::operator()(const State &x, State &dxdt, const Scalar t) {
  using namespace Eigen;
  const auto N = param.N;
  const long long int dt_grid_begin = lm_size * grid_size;
    
  // dxdt.head(grid_size) = x.tail(grid_size);
  // free_propagation_new(x, dxdt);
  // dxdt.tail(grid_size) += -x.head(grid_size) * V + Q(t);

  // psi_lm;
  // dt_psi_lm;
  Vector dr_psi_lm;
  Vector drdr_psi_lm;
  
  // For each spherical harmonic, set the second order time derivative
  for(int lm = 0; lm < lm_size; ++lm){
    auto dtdt_psi_lm = dxdt(seqN(dt_grid_begin + lm * grid_size, grid_size));
    dtdt_psi_lm = 0;
    // for(auto [lm1, idx1] : psi_lm_map[lm]){
    //   dtdt_psi_lm -= coeff[idx1] * x(seqN(lm1 * grid_size, grid_size));
    // }
    for(auto [lm1, idx1] : dt_psi_lm_map[lm]){
      dtdt_psi_lm -= coeffs[idx1] * x(seqN(dt_grid_begin + lm1 * grid_size, grid_size));
    }
  //   for(auto [lm1, idx1] : dr_psi_lm_map[lm]){
  //     dtdt_psi_lm -= coeff[idx1] * dr_psi_lm(seqN(lm1 * grid_size, grid_size));
  //   }
  //   for(auto [lm1, idx1] : drdr_psi_lm_map[lm]){
  //     dtdt_psi_lm -= coeff[idx1] * drdr_psi_lm(seqN(lm1 * grid_size, grid_size));
  //   }
  }

  
}

// void TeukolskyScalarPDE::free_propagation_new(const Vector &x, Vector &dxdt) const {
//     using namespace Eigen;
//     const auto N = param.N;
//     const Scalar r_min = param.rast_min;
//     const Scalar r_max = param.rast_max;
//     const Scalar h = (r_max - r_min) / (N - 1); // get_h(param.rast_min, param.rast_max, N);
    
//     const long long int grid_size = N + 1;
//     const long long int grid_begin = 0;
//     const long long int dt_grid_begin = grid_size;
    
//     dxdt(dt_grid_begin+0)
//       = (-25 * x(dt_grid_begin+0) + 48 * x(dt_grid_begin+1) - 36 * x(dt_grid_begin+2) + 16 * x(dt_grid_begin+3) - 3 * x(dt_grid_begin+4) ) / (12*h);
//     dxdt(dt_grid_begin+1)
//       = (11 * x(grid_begin+0) - 20 * x(grid_begin+1) + 6 * x(grid_begin+2) + 4 * x(grid_begin+3) - x(grid_begin+4) ) / (12*h*h);
//     dxdt(seqN(dt_grid_begin+2,grid_size-4))
//       = (-1 * x(seqN(grid_begin,grid_size-4)) + 16 * x(seqN(grid_begin+1,grid_size-4)) - 30 * x(seqN(grid_begin+2,grid_size-4)) + 16 * x(seqN(grid_begin+3,grid_size-4)) - 1 * x(seqN(grid_begin+4,grid_size-4)) ) / (12*h*h);
//     dxdt(dt_grid_begin+grid_size-2)
//       = (-1 * x(grid_begin+grid_size-1-4) + 4 * x(grid_begin+grid_size-1-3) + 6 * x(grid_begin+grid_size-1-2) - 20 * x(grid_begin+grid_size-1-1) + 11 * x(grid_begin+grid_size-1-0)  ) / (12*h*h);
//     dxdt(dt_grid_begin+grid_size-1)
//       = (-3 * x(dt_grid_begin+grid_size-1-4) + 16 * x(dt_grid_begin+grid_size-1-3) - 36 * x(dt_grid_begin+grid_size-1-2) + 48 * x(dt_grid_begin+grid_size-1-1) - 25 * x(dt_grid_begin+grid_size-1-0)  ) / (12*h);
//   }



// TeukolskyScalarPDE::Vector TeukolskyScalarPDE::compute_r_ast_vector(const Scalar r_min, const Scalar r_max, const long long int N) {
//   Scalar h = (r_max - r_min) / (N - 1);
//   Vector r_ast(N + 1);
//   for(long long int i = 0; i < N + 1; ++i) {
//     r_ast[i] = r_min + i * h - h / Scalar(2);
//   }
//   return r_ast;
// }

#include "teukolsky.hpp"

#include <thread>
#include <boost/lockfree/queue.hpp>

#include "teukolsky.gen"

namespace Teukolsky
{

  CouplingInfo make_coupling_info_map(const std::vector<std::vector<std::pair<long long int, long long int>>> &info, const long long int l_max) {
    const long long int lm_size = (l_max + 1) * (l_max + 1);
    CouplingInfo result;
    result.resize(lm_size);
    for(int idx = 0; idx < lm_size; ++idx){
      for(auto p : info[idx]){
	auto [idx1, coeff1] = p;
	if(idx1 < lm_size){
	  result[idx][idx1] = coeff1;
	}
      }
    }
    return result;
  }

  HighPrecisionVector compute_hp_r_vector(const HighPrecisionScalar rast_min, const HighPrecisionScalar rast_max, const long long int N, const HighPrecisionScalar M, const HighPrecisionScalar a) {
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
    HighPrecisionVector r(N + 1);
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

  Vector compute_r_vector(const Scalar rast_min, const Scalar rast_max, const long long int N, const Scalar M, const Scalar a) {
    return compute_hp_r_vector(static_cast<HighPrecisionScalar>(rast_min),
			       static_cast<HighPrecisionScalar>(rast_max),
			       N,
			       static_cast<HighPrecisionScalar>(M),
			       static_cast<HighPrecisionScalar>(a)).cast<double>();
  }

 
}

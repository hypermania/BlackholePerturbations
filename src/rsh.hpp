/*!
  \file rsh_expressions.hpp
  \brief Generation of Eigen expressions involving Real Spherical Harmonic (RSH) cross terms.
  
*/
#ifndef RSH_EXPRESSIONS_HPP
#define RSH_EXPRESSIONS_HPP

#include <Eigen/Dense>
#include "coeffs.hpp"

namespace RSH {

  constexpr long long int lm_to_idx(const long long int l, const long long int m) {
    return l*l + (m-(-l));
  }

  constexpr std::pair<long long int, long long int> idx_to_lm(const long long int idx) {
    const long long int l = static_cast<long long int>(sqrt(static_cast<double>(idx)));
    const long long int m = idx - l * l - l;
    return std::make_pair(l, m);
  }
  
  template<long long int l_max, long long int lm_idx, long long int term_idx>
  constexpr auto quadratic_rsh_expression_accumulator(const Eigen::ArrayXd &x, const long long int grid_size) {
    using namespace Eigen;
    if constexpr(term_idx >= sqr_rsh_terms.size()) {
      return ArrayXd::Zero(grid_size);
    } else {
      
      constexpr long long int lm_0 = std::get<0>(sqr_rsh_terms[term_idx]);
      // constexpr long long int l_0 = idx_to_lm(lm_0).first;

      // Assuming that sqr_rsh_terms is sorted in lm_0 (first tuple elem),
      // we terminate the recursion immediately when lm_0 > lm_idx
      if constexpr(lm_0 > lm_idx) {
	// static_assert(lm_0 == 0, "print lm_0");
	return ArrayXd::Zero(grid_size);
      } else {

	if constexpr(lm_0 != lm_idx) {
	  return quadratic_rsh_expression_accumulator<l_max, lm_idx, term_idx + 1>(x, grid_size);
	} else {
	  
	  const long long int lm_1 = std::get<1>(sqr_rsh_terms[term_idx]);
	  const long long int lm_2 = std::get<2>(sqr_rsh_terms[term_idx]);
	  const double coupling = std::get<3>(sqr_rsh_terms[term_idx]);
	  const long long int l_1 = idx_to_lm(lm_1).first;
	  const long long int l_2 = idx_to_lm(lm_2).first;
	  
	  if constexpr(l_1 <= l_max && l_2 <= l_max) {
	    return coupling * x(seqN(lm_1 * grid_size, grid_size)) * x(seqN(lm_2 * grid_size, grid_size))
	      + quadratic_rsh_expression_accumulator<l_max, lm_idx, term_idx + 1>(x, grid_size);
	  } else {
	    return quadratic_rsh_expression_accumulator<l_max, lm_idx, term_idx + 1>(x, grid_size);
	  }
	  
	}
      }
    }
  }

  template<long long int l_max, long long int lm_idx>
  constexpr auto quadratic_rsh_expression(const Eigen::ArrayXd &x, const long long int grid_size) {
    static_assert(l_max <= l_max_for_terms, "l_max is higher than the currently supported value.");
    return quadratic_rsh_expression_accumulator<l_max, lm_idx, 0>(x, grid_size);
  }

}

#endif

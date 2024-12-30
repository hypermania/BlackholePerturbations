/*!
  \file discretization.hpp
  \brief Utilities for spatial discretization.
  
*/
#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

constexpr double get_h(const double r_min, const double r_max, const long long int N) {
  const double h = (r_max - r_min) / (N - 1);
  return h;
}

constexpr double i_to_r_ast(const double r_min, const double r_max, const long long int N, const long long int i) {
  const double h = (r_max - r_min) / (N - 1);
  return r_min + i * h - h / 2.0;
}

constexpr long long int r_ast_to_i(const double r_min, const double r_max, const long long int N, const double r) {
  const double h = (r_max - r_min) / (N - 1);
  return (long long int)(0.5 + (r - r_min) / h);
}

#endif

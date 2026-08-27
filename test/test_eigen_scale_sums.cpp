#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

#include <omp.h>

#include <Eigen/Dense>
#include <boost/multiprecision/float128.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

#include "odeint_eigen/eigen_operations.hpp"

namespace {

int failures = 0;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

template <class State>
bool values_equal(const State &actual, const State &expected) {
  if (actual.size() != expected.size()) return false;
  for (Eigen::Index i = 0; i < actual.size(); ++i) {
#ifdef __FAST_MATH__
    using std::abs;
    using Scalar = typename State::Scalar;
    const Scalar scale = Scalar(1) + abs(expected[i]);
    if (abs(actual[i] - expected[i]) >
        Scalar(64) * std::numeric_limits<Scalar>::epsilon() * scale) {
      return false;
    }
#else
    if (actual[i] != expected[i]) return false;
#endif
  }
  return true;
}

template <class State>
bool values_exactly_equal(const State &actual, const State &expected) {
  if (actual.size() != expected.size()) return false;
  for (Eigen::Index i = 0; i < actual.size(); ++i) {
    if (actual[i] != expected[i]) return false;
  }
  return true;
}

template <class State>
void check_scale_sums(Eigen::Index size, int threads, const std::string &label) {
  using Scalar = typename State::Scalar;
  using Operations = boost::numeric::odeint::eigen_operations<State>;
  using LegacyOperations = boost::numeric::odeint::eigen_operations<State, false>;
  omp_set_num_threads(threads);

  std::array<State, 7> input;
  for (State &vector : input) vector.resize(size);
  for (Eigen::Index i = 0; i < size; ++i) {
    for (int vector = 0; vector < 7; ++vector) {
      const long long pattern = ((i * (17 + vector * 2) + vector * 11) % 97) - 48;
      input[vector][i] = Scalar(pattern) / Scalar(37 + vector);
    }
  }
  const std::array<Scalar, 7> alpha{
      Scalar(2) / Scalar(7), Scalar(-3) / Scalar(11),
      Scalar(5) / Scalar(13), Scalar(-7) / Scalar(17),
      Scalar(11) / Scalar(19), Scalar(-13) / Scalar(23),
      Scalar(17) / Scalar(29)};

  State actual(size);
  State expected(size);

  {
    typename Operations::template scale_sum1<> operation(alpha[0]);
    operation(actual, input[0]);
    for (Eigen::Index i = 0; i < size; ++i) expected[i] = alpha[0] * input[0][i];
    require(values_equal(actual, expected), label + " scale_sum1");
  }
  {
    typename Operations::template scale_sum2<> operation(alpha[0], alpha[1]);
    operation(actual, input[0], input[1]);
    for (Eigen::Index i = 0; i < size; ++i)
      expected[i] = alpha[0] * input[0][i] + alpha[1] * input[1][i];
    require(values_equal(actual, expected), label + " scale_sum2");

    State aliased = input[0];
    operation(aliased, aliased, input[1]);
    require(values_equal(aliased, expected), label + " scale_sum2 aliased output");
  }
  {
    typename Operations::template scale_sum3<> operation(alpha[0], alpha[1], alpha[2]);
    operation(actual, input[0], input[1], input[2]);
    for (Eigen::Index i = 0; i < size; ++i)
      expected[i] = alpha[0] * input[0][i] + alpha[1] * input[1][i]
                    + alpha[2] * input[2][i];
    require(values_equal(actual, expected), label + " scale_sum3");
  }
  {
    typename Operations::template scale_sum4<> operation(
        alpha[0], alpha[1], alpha[2], alpha[3]);
    operation(actual, input[0], input[1], input[2], input[3]);
    for (Eigen::Index i = 0; i < size; ++i)
      expected[i] = alpha[0] * input[0][i] + alpha[1] * input[1][i]
                    + alpha[2] * input[2][i] + alpha[3] * input[3][i];
    require(values_equal(actual, expected), label + " scale_sum4");
  }
  {
    typename Operations::template scale_sum5<> operation(
        alpha[0], alpha[1], alpha[2], alpha[3], alpha[4]);
    operation(actual, input[0], input[1], input[2], input[3], input[4]);
    for (Eigen::Index i = 0; i < size; ++i)
      expected[i] = alpha[0] * input[0][i] + alpha[1] * input[1][i]
                    + alpha[2] * input[2][i] + alpha[3] * input[3][i]
                    + alpha[4] * input[4][i];
    require(values_equal(actual, expected), label + " scale_sum5");
  }
  {
    typename Operations::template scale_sum6<> operation(
        alpha[0], alpha[1], alpha[2], alpha[3], alpha[4], alpha[5]);
    operation(actual, input[0], input[1], input[2], input[3], input[4], input[5]);
    for (Eigen::Index i = 0; i < size; ++i)
      expected[i] = alpha[0] * input[0][i] + alpha[1] * input[1][i]
                    + alpha[2] * input[2][i] + alpha[3] * input[3][i]
                    + alpha[4] * input[4][i] + alpha[5] * input[5][i];
    require(values_equal(actual, expected), label + " scale_sum6");
  }
  {
    typename Operations::template scale_sum7<> operation(
        alpha[0], alpha[1], alpha[2], alpha[3], alpha[4], alpha[5], alpha[6]);
    operation(actual, input[0], input[1], input[2], input[3], input[4], input[5],
              input[6]);
    for (Eigen::Index i = 0; i < size; ++i)
      expected[i] = alpha[0] * input[0][i] + alpha[1] * input[1][i]
                    + alpha[2] * input[2][i] + alpha[3] * input[3][i]
                    + alpha[4] * input[4][i] + alpha[5] * input[5][i]
                    + alpha[6] * input[6][i];
    require(values_equal(actual, expected), label + " scale_sum7");

    State aliased = input[0];
    operation(aliased, aliased, input[1], input[2], input[3], input[4], input[5],
              input[6]);
    require(values_equal(aliased, expected), label + " scale_sum7 aliased output");
  }


  State legacy(size);
  const Scalar one = 1;
  {
    typename Operations::template scale_sum1<> optimized(one);
    typename LegacyOperations::template scale_sum1<> reference(one);
    optimized(actual, input[0]);
    reference(legacy, input[0]);
    require(values_exactly_equal(actual, legacy), label + " unit scale_sum1 exact");
  }
  {
    typename Operations::template scale_sum2<> optimized(one, alpha[1]);
    typename LegacyOperations::template scale_sum2<> reference(one, alpha[1]);
    optimized(actual, input[0], input[1]);
    reference(legacy, input[0], input[1]);
    require(values_exactly_equal(actual, legacy), label + " unit scale_sum2 exact");
  }
  {
    typename Operations::template scale_sum3<> optimized(one, alpha[1], alpha[2]);
    typename LegacyOperations::template scale_sum3<> reference(one, alpha[1], alpha[2]);
    optimized(actual, input[0], input[1], input[2]);
    reference(legacy, input[0], input[1], input[2]);
    require(values_exactly_equal(actual, legacy), label + " unit scale_sum3 exact");
  }
  {
    typename Operations::template scale_sum4<> optimized(
        one, alpha[1], alpha[2], alpha[3]);
    typename LegacyOperations::template scale_sum4<> reference(
        one, alpha[1], alpha[2], alpha[3]);
    optimized(actual, input[0], input[1], input[2], input[3]);
    reference(legacy, input[0], input[1], input[2], input[3]);
    require(values_exactly_equal(actual, legacy), label + " unit scale_sum4 exact");
  }
  {
    typename Operations::template scale_sum5<> optimized(
        one, alpha[1], alpha[2], alpha[3], alpha[4]);
    typename LegacyOperations::template scale_sum5<> reference(
        one, alpha[1], alpha[2], alpha[3], alpha[4]);
    optimized(actual, input[0], input[1], input[2], input[3], input[4]);
    reference(legacy, input[0], input[1], input[2], input[3], input[4]);
    require(values_exactly_equal(actual, legacy), label + " unit scale_sum5 exact");
  }
  {
    typename Operations::template scale_sum6<> optimized(
        one, alpha[1], alpha[2], alpha[3], alpha[4], alpha[5]);
    typename LegacyOperations::template scale_sum6<> reference(
        one, alpha[1], alpha[2], alpha[3], alpha[4], alpha[5]);
    optimized(actual, input[0], input[1], input[2], input[3], input[4], input[5]);
    reference(legacy, input[0], input[1], input[2], input[3], input[4], input[5]);
#ifdef __FAST_MATH__
    // GCC 11 and GCC 15 choose different permitted reassociation trees for
    // the compiler-generated legacy expression.  The optimized path fixes a
    // stable tree explicitly, so under fast-math the portable requirement is
    // binary128 roundoff equivalence rather than compiler-version-dependent
    // bitwise identity.
    require(values_equal(actual, legacy),
            label + " unit scale_sum6 roundoff equivalence");
#else
    require(values_exactly_equal(actual, legacy), label + " unit scale_sum6 exact");
#endif
  }
  {
    typename Operations::template scale_sum7<> optimized(
        one, alpha[1], alpha[2], alpha[3], alpha[4], alpha[5], alpha[6]);
    typename LegacyOperations::template scale_sum7<> reference(
        one, alpha[1], alpha[2], alpha[3], alpha[4], alpha[5], alpha[6]);
    optimized(actual, input[0], input[1], input[2], input[3], input[4], input[5],
              input[6]);
    reference(legacy, input[0], input[1], input[2], input[3], input[4], input[5],
              input[6]);
    require(values_exactly_equal(actual, legacy), label + " unit scale_sum7 exact");
  }
}

}  // namespace

int main() {
  using DoubleArray = Eigen::Array<double, -1, 1>;
  using QuadArray = Eigen::Array<boost::multiprecision::float128, -1, 1>;
  using DoubleMatrix = Eigen::Matrix<double, -1, 1>;
  using QuadMatrix = Eigen::Matrix<boost::multiprecision::float128, -1, 1>;
  omp_set_dynamic(0);

  check_scale_sums<DoubleArray>(31, 6, "double array serial threshold");
  check_scale_sums<DoubleArray>(4095, 6, "double array just below threshold");
  check_scale_sums<DoubleArray>(4096, 2, "double array threshold");
  check_scale_sums<DoubleArray>(5003, 6, "double array uneven parallel partition");
  check_scale_sums<QuadArray>(31, 6, "binary128 array serial threshold");
  check_scale_sums<QuadArray>(4096, 2, "binary128 array threshold");
  check_scale_sums<QuadArray>(5003, 6, "binary128 array uneven parallel partition");
  check_scale_sums<DoubleMatrix>(31, 6, "double matrix serial threshold");
  check_scale_sums<DoubleMatrix>(4096, 6, "double matrix threshold");
  check_scale_sums<QuadMatrix>(4096, 6, "binary128 matrix threshold");

  if (failures != 0) {
    std::cerr << failures << " Eigen scale-sum checks failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: scale_sum1 through scale_sum7 for Eigen arrays and matrices, "
               "double and binary128\n";
  return EXIT_SUCCESS;
}

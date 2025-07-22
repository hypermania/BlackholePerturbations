#include <vector>
#include <algorithm>
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
//#include "cubic_scalar.hpp"
#include "boost/type_index.hpp"

#include "teukolsky_scalar.hpp"

#include "examples.hpp"


int main(int argc, char **argv) {
  TeukolskyScalarPDEParam param;
  param.s = 0; //.convert_to<double>();
  param.l_max = 1; //.convert_to<double>();
  param.M = 1.0;
  param.a = 0.1;

  param.rast_min = -10;
  param.rast_max = 10;
  param.N = 1000;
  
  param.t_start = 0; //.convert_to<double>();
  param.t_end = 10; //.convert_to<double>();
  param.t_interval = 0.5;
  param.delta_t = 0.1; //.convert_to<double>();

  TeukolskyScalarPDE eqn(param);
  typedef TeukolskyScalarPDE::HighPrecisionScalar HP;

  auto c = std::complex<HP>(HP(1.0), HP(1.0));
  
  //run_coupled_eqn();
  //run_sourced_eqn();

  return 0;
}


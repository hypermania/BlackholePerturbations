/*! 
  \file teukolsky_scalar_cuda.cuh
  \author Siyang Ling
  \brief Header for Teukolsky equations that runs on the GPU.

*/
#ifndef TEUKOLSKY_SCALAR_CUDA_CUH
#define TEUKOLSKY_SCALAR_CUDA_CUH

#include <Eigen/Dense>
// #include <boost/math/special_functions/lambert_w.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
// #include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/eigen.hpp>
//#include "discretization.hpp"

#include "teukolsky_scalar.hpp"
#include "teukolsky.hpp"
#include <thrust/device_vector.h>
#include <thrust/complex.h>
#include <cuda_runtime.h>


/*! 
  \brief Teukolsky equation with \f$ s=0 \f$ that runs with CUDA.
*/
struct CudaTeukolskyScalarPDE {
  typedef TeukolskyScalarPDEParam Param;
  typedef double Scalar;
  typedef thrust::complex<Scalar> ComplexScalar;
  typedef thrust::device_vector<double> Vector;
  typedef thrust::device_vector<thrust::complex<double>> ComplexVector;

  // typedef boost::multiprecision::cpp_bin_float_100 HighPrecisionScalar;
  // typedef std::complex<HighPrecisionScalar> HighPrecisionComplex;
  // typedef Eigen::Array<HighPrecisionScalar, -1, 1> HighPrecisionVector;
  // typedef Eigen::Array<HighPrecisionComplex, -1, 1> HighPrecisionComplexVector;

  typedef Teukolsky::HighPrecisionScalar HighPrecisionScalar;
  typedef Teukolsky::HighPrecisionComplex HighPrecisionComplex;
  typedef Teukolsky::HighPrecisionVector HighPrecisionVector;
  typedef Teukolsky::HighPrecisionComplexVector HighPrecisionComplexVector;

  typedef ComplexVector State;
  
  Param param;
  // Vector V;
  // std::function<Vector(const Scalar)> Q;
  long long int grid_size;
  long long int lm_size;

  typedef Teukolsky::CouplingInfo CouplingInfo;
  CouplingInfo psi_lm_map;
  CouplingInfo dr_psi_lm_map;
  CouplingInfo drdr_psi_lm_map;
  CouplingInfo dt_psi_lm_map;

  std::vector<ComplexVector> coeffs;
  ComplexVector drdr_psi_lm;
  ComplexVector dr_psi_lm;

  // CUDA graph for evaluating the system ( for operator() )
  cudaGraphExec_t system_graph_exec;
  
  CudaTeukolskyScalarPDE(Param param_);
  // ~CudaTeukolskyScalarPDE();

  /*!
    \brief The function called by odeint library.
    \param[in] x The current state of the system.
    \param[out] dxdt The time derivative, dxdt of the system.
    \param t The current time parameter.
  */
  void operator()(const State &x, State &dxdt, const Scalar t);
  
  //void compute_derivatives(const State &x, ComplexVector &dr_psi_lm, ComplexVector &drdr_psi_lm) const;

  cudaGraph_t prepare_cuda_graph(const State &x, State &dxdt);

};


#endif

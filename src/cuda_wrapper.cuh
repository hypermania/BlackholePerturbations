/*! 
  \file cuda_wrapper.cuh
  \author Siyang Ling
  \brief Wrapper for CUDA Toolkit.
*/
#ifndef CUDA_WRAPPER_CUH
#define CUDA_WRAPPER_CUH

#include <iostream>

#include <Eigen/Dense>

#include <thrust/device_vector.h>
#include <thrust/complex.h>
// #include <thrust/host_vector.h>
// #include <thrust/execution_policy.h>
// #include <thrust/reduce.h>
// #include <thrust/functional.h>
// #include <thrust/fill.h>
// #include <thrust/transform.h>

#include "cufft.h"
#include "cufftXt.h"
#include <cuda_runtime.h>



typedef decltype(Eigen::VectorXd().begin()) eigen_iterator;
typedef decltype(thrust::device_vector<double>().begin()) thrust_iterator;
typedef thrust::detail::normal_iterator<thrust::device_ptr<const double>> thrust_const_iterator;
typedef Eigen::internal::pointer_based_stl_iterator<Eigen::Matrix<double, -1, 1>> eigen_iterator_2;


/*
  Explicit template instantiation declarations for the thrust library.
  They are declared here so that they are instantiatiated in cuda_wrapper.cu (and compiled with nvcc),
  and don't get instantiated in other translation units.
  This is necessary since we want to call thrust functions in translation units compiled by other compilers (g++ / icpx).
*/
extern template class thrust::device_vector<double>;
extern template class thrust::device_ptr<double>;
extern template thrust::device_ptr<double> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<double>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<double>>);
extern template eigen_iterator thrust::copy(const thrust::detail::execution_policy_base<thrust::cuda_cub::cross_system<thrust::cuda_cub::tag, thrust::system::cpp::detail::tag>> &, thrust_const_iterator, thrust_const_iterator, eigen_iterator);

extern template thrust_iterator thrust::copy(eigen_iterator, eigen_iterator, thrust_iterator);
extern template eigen_iterator thrust::copy(thrust_iterator, thrust_iterator, eigen_iterator);

template<typename Derived, typename ThrustScalar>
void copy_vector(Eigen::PlainObjectBase<Derived> &out, const thrust::device_vector<ThrustScalar> &in)
{
  assert(out.size() * sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar) >= in.size() * sizeof(ThrustScalar));
  cudaMemcpy((void *)out.data(), (const void *)thrust::raw_pointer_cast(in.data()), in.size() * sizeof(ThrustScalar), cudaMemcpyDeviceToHost);
}

template<typename ThrustScalar, typename Derived>
void copy_vector(thrust::device_vector<ThrustScalar> &out, const Eigen::PlainObjectBase<Derived> &in)
{
  assert(out.size() * sizeof(ThrustScalar) >= in.size() * sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar));
  cudaMemcpy((void *)thrust::raw_pointer_cast(out.data()), (const void *)in.data(), in.size() * sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar), cudaMemcpyHostToDevice);
}

extern template void copy_vector(thrust::device_vector<thrust::complex<double>> &out, const Eigen::PlainObjectBase<Eigen::ArrayXcd> &in);
extern template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXcd> &out, const thrust::device_vector<thrust::complex<double>> &in);
extern template void copy_vector(thrust::device_vector<double> &out, const Eigen::PlainObjectBase<Eigen::ArrayXd> &in);
extern template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXd> &out, const thrust::device_vector<double> &in);


void show_gpu_memory_usage(void);


/*
  new declarations
 */
extern template class thrust::device_vector<thrust::complex<double>>;
// extern template class thrust::device_ptr<thrust::complex<double>>;
// extern template class thrust::device_reference<thrust::complex<double>>;
// extern template class thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>;
extern template thrust::device_ptr<thrust::complex<double>> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<thrust::complex<double>>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<thrust::complex<double>>>);
extern template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::strip_const_get_value<thrust::cuda_cub::tag>(thrust::cuda_cub::tag const &) const;
extern template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::convert_to_value_type<thrust::cuda_cub::tag>(thrust::cuda_cub::tag *) const;
extern template thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::operator thrust::complex<double>() const;

#endif

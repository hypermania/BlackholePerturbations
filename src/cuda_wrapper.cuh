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
#include <cuda_runtime.h>

#include <cufile.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>



void show_gpu_memory_usage(void);
void show_cuda_info(void);


/*
  Functions for copying between CUDA device and CPU memory.
 */
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

template<typename Scalar, typename ThrustScalar>
void copy_vector(std::vector<Scalar> &out, const thrust::device_vector<ThrustScalar> &in)
{
  assert(out.size() * sizeof(Scalar) >= in.size() * sizeof(ThrustScalar));
  cudaMemcpy((void *)out.data(), (const void *)thrust::raw_pointer_cast(in.data()), in.size() * sizeof(ThrustScalar), cudaMemcpyDeviceToHost);
}

template<typename ThrustScalar, typename Scalar>
void copy_vector(thrust::device_vector<ThrustScalar> &out, const std::vector<Scalar> &in)
{
  assert(out.size() * sizeof(ThrustScalar) >= in.size() * sizeof(Scalar));
  cudaMemcpy((void *)thrust::raw_pointer_cast(out.data()), (const void *)in.data(), in.size() * sizeof(Scalar), cudaMemcpyHostToDevice);
}

extern template void copy_vector(thrust::device_vector<thrust::complex<double>> &out, const Eigen::PlainObjectBase<Eigen::ArrayXcd> &in);
extern template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXcd> &out, const thrust::device_vector<thrust::complex<double>> &in);
extern template void copy_vector(thrust::device_vector<double> &out, const Eigen::PlainObjectBase<Eigen::ArrayXd> &in);
extern template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXd> &out, const thrust::device_vector<double> &in);
extern template void copy_vector(thrust::device_vector<thrust::complex<double>> &out, const std::vector<std::complex<double>> &in);
extern template void copy_vector(std::vector<std::complex<double>> &out, const thrust::device_vector<thrust::complex<double>> &in);
extern template void copy_vector(thrust::device_vector<double> &out, const std::vector<double> &in);
extern template void copy_vector(std::vector<double> &out, const thrust::device_vector<double> &in);


/*
  Wrapper for low level CUDA functions.
*/
void cudaMemcpy2DWrapper(void* dst, size_t dpitch, const void* src, size_t spitch, size_t width, size_t height);

/*
  Functions for copying between CUDA device and storage.
 */
template<typename ThrustScalar>
void gpudirect_write(const thrust::device_vector<ThrustScalar> &vec, const std::string &filename)
{
  CUfileHandle_t cfHandle;
  CUfileDescr_t cfDescr = {};
  int fd = open(filename.data(), O_CREAT | O_RDWR, 0664);
  // assert(fd >= 0, "file open failed");
  
  // Set up GDS descriptor
  cfDescr.handle.fd = fd;
  cfDescr.type = CU_FILE_HANDLE_TYPE_OPAQUE_FD;
  CUfileError_t status = cuFileHandleRegister(&cfHandle, &cfDescr);
  // assert(status.err == CU_FILE_SUCCESS, "cuFileHandleRegister failed");
  
  // Perform the write
  ssize_t bytes = cuFileWrite(cfHandle,
			      reinterpret_cast<const void *>(thrust::raw_pointer_cast(vec.data())),
			      sizeof(ThrustScalar) * vec.size(),
			      0, 0);
  // assert(bytes >= 0, "GPUDirect file write failed");
  
  // Clean up
  cuFileHandleDeregister(cfHandle);
  close(fd);
}

template<typename ThrustScalar>
void gpudirect_read(thrust::device_vector<ThrustScalar> &vec, const std::string &filename)
{
  CUfileHandle_t cfHandle;
  CUfileDescr_t cfDescr = {};
  int fd = open(filename.data(), O_RDONLY, 0);
  // assert(fd >= 0, "file open failed");

  struct stat sb;
  fstat(fd, &sb);
  const size_t size = sb.st_size;
  
  // Set up GDS descriptor
  cfDescr.handle.fd = fd;
  cfDescr.type = CU_FILE_HANDLE_TYPE_OPAQUE_FD;
  CUfileError_t status = cuFileHandleRegister(&cfHandle, &cfDescr);
  // assert(status.err == CU_FILE_SUCCESS, "cuFileHandleRegister failed");
  
  // Perform the write
  vec.resize(size / sizeof(ThrustScalar));
  ssize_t bytes = cuFileRead(cfHandle,
			     reinterpret_cast<void *>(thrust::raw_pointer_cast(vec.data())),
			     size,
			     0, 0);
  // assert(bytes >= 0, "GPUDirect file read failed");
  
  // Clean up
  cuFileHandleDeregister(cfHandle);
  close(fd);
}

extern template void gpudirect_write(const thrust::device_vector<double> &vec, const std::string &filename);
extern template void gpudirect_write(const thrust::device_vector<thrust::complex<double>> &vec, const std::string &filename);
extern template void gpudirect_read(thrust::device_vector<double> &vec, const std::string &filename);
extern template void gpudirect_read(thrust::device_vector<thrust::complex<double>> &vec, const std::string &filename);


/*
  Explicit template instantiation declarations for the thrust library.
  They are declared here so that they are instantiatiated in cuda_wrapper.cu (and compiled with nvcc),
  and don't get instantiated in other translation units.
  This is necessary since we want to call thrust functions in translation units compiled by other compilers (g++ / icpx).
*/
typedef decltype(Eigen::VectorXd().begin()) eigen_iterator;
typedef decltype(thrust::device_vector<double>().begin()) thrust_iterator;
typedef thrust::detail::normal_iterator<thrust::device_ptr<const double>> thrust_const_iterator;
typedef Eigen::internal::pointer_based_stl_iterator<Eigen::Matrix<double, -1, 1>> eigen_iterator_2;


extern template class thrust::device_vector<double>;
extern template class thrust::device_ptr<double>;
extern template thrust::device_ptr<double> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<double>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<double>>);
extern template eigen_iterator thrust::copy(const thrust::detail::execution_policy_base<thrust::cuda_cub::cross_system<thrust::cuda_cub::tag, thrust::system::cpp::detail::tag>> &, thrust_const_iterator, thrust_const_iterator, eigen_iterator);


extern template thrust_iterator thrust::copy(eigen_iterator, eigen_iterator, thrust_iterator);
extern template eigen_iterator thrust::copy(thrust_iterator, thrust_iterator, eigen_iterator);


extern template class thrust::device_vector<thrust::complex<double>>;
// extern template class thrust::device_ptr<thrust::complex<double>>;
// extern template class thrust::device_reference<thrust::complex<double>>;
// extern template class thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>;
extern template thrust::device_ptr<thrust::complex<double>> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<thrust::complex<double>>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<thrust::complex<double>>>);
extern template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::strip_const_get_value<thrust::cuda_cub::tag>(thrust::cuda_cub::tag const &) const;
extern template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::convert_to_value_type<thrust::cuda_cub::tag>(thrust::cuda_cub::tag *) const;
extern template thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::operator thrust::complex<double>() const;

#endif

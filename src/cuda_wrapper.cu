#include "cuda_wrapper.cuh"

template void copy_vector(thrust::device_vector<thrust::complex<double>> &out, const Eigen::PlainObjectBase<Eigen::ArrayXcd> &in);
template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXcd> &out, const thrust::device_vector<thrust::complex<double>> &in);
template void copy_vector(thrust::device_vector<double> &out, const Eigen::PlainObjectBase<Eigen::ArrayXd> &in);
template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXd> &out, const thrust::device_vector<double> &in);
template void copy_vector(thrust::device_vector<thrust::complex<double>> &out, const std::vector<std::complex<double>> &in);
template void copy_vector(std::vector<std::complex<double>> &out, const thrust::device_vector<thrust::complex<double>> &in);
template void copy_vector(thrust::device_vector<double> &out, const std::vector<double> &in);
template void copy_vector(std::vector<double> &out, const thrust::device_vector<double> &in);

/*
  Wrapper for low level CUDA functions.
*/
void cudaMemcpy2DWrapper(void* dst, size_t dpitch, const void* src, size_t spitch, size_t width, size_t height)
{
  cudaMemcpy2D(dst, dpitch, src, spitch, width, height, cudaMemcpyDefault);
}



template void gpudirect_write(const thrust::device_vector<double> &vec, const std::string &filename);
template void gpudirect_write(const thrust::device_vector<thrust::complex<double>> &vec, const std::string &filename);
template void gpudirect_read(thrust::device_vector<double> &vec, const std::string &filename);
template void gpudirect_read(thrust::device_vector<thrust::complex<double>> &vec, const std::string &filename);

/*
  Thrust library template instantiation definitions
 */
template class thrust::device_vector<double>;
template class thrust::device_ptr<double>;
template thrust::device_ptr<double> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<double>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<double>>);
template eigen_iterator thrust::copy(const thrust::detail::execution_policy_base<thrust::cuda_cub::cross_system<thrust::cuda_cub::tag, thrust::system::cpp::detail::tag>> &, thrust_const_iterator, thrust_const_iterator, eigen_iterator);

template thrust_iterator thrust::copy(eigen_iterator, eigen_iterator, thrust_iterator);
template eigen_iterator thrust::copy(thrust_iterator, thrust_iterator, eigen_iterator);

template class thrust::device_vector<thrust::complex<double>>;
// template class thrust::device_ptr<thrust::complex<double>>;
// template class thrust::device_reference<thrust::complex<double>>;
// template class thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>;
template thrust::device_ptr<thrust::complex<double>> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<thrust::complex<double>>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<thrust::complex<double>>>);
template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::strip_const_get_value<thrust::cuda_cub::tag>(thrust::cuda_cub::tag const &) const;
template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::convert_to_value_type<thrust::cuda_cub::tag>(thrust::cuda_cub::tag *) const;
template thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::operator thrust::complex<double>() const;

void show_gpu_memory_usage(void)
{
  size_t free, total;
  cudaMemGetInfo(&free, &total);
  std::cout << "free / total = "
	    << free << " B / " << total << " B ("
	    << free / (1024*1024) << " MB / " << total / (1024*1024) << " MB)\n";
}

void show_cuda_info(void)
{
  int cudaDevAttrMaxThreadsPerBlock_val;
  int cudaDevAttrMaxSharedMemoryPerBlock_val;
  int cudaDevAttrTotalConstantMemory_val;
  int cudaDevAttrClockRate_val;
  
  int cudaDevAttrMemoryClockRate_val;
  int cudaDevAttrGlobalMemoryBusWidth_val;
  int cudaDevAttrL2CacheSize_val;
  
  int cudaDevAttrMaxSharedMemoryPerMultiprocessor_val;
  int cudaDevAttrMaxRegistersPerMultiprocessor_val;
  int cudaDevAttrMaxBlocksPerMultiprocessor_val;
  
  cudaDeviceGetAttribute(&cudaDevAttrMaxThreadsPerBlock_val, cudaDevAttrMaxThreadsPerBlock, 0);
  cudaDeviceGetAttribute(&cudaDevAttrMaxSharedMemoryPerBlock_val, cudaDevAttrMaxSharedMemoryPerBlock, 0);
  cudaDeviceGetAttribute(&cudaDevAttrTotalConstantMemory_val, cudaDevAttrTotalConstantMemory, 0);
  cudaDeviceGetAttribute(&cudaDevAttrClockRate_val, cudaDevAttrClockRate, 0);

  cudaDeviceGetAttribute(&cudaDevAttrMemoryClockRate_val, cudaDevAttrMemoryClockRate, 0);
  cudaDeviceGetAttribute(&cudaDevAttrGlobalMemoryBusWidth_val, cudaDevAttrGlobalMemoryBusWidth, 0);
  cudaDeviceGetAttribute(&cudaDevAttrL2CacheSize_val, cudaDevAttrL2CacheSize, 0);
  
  cudaDeviceGetAttribute(&cudaDevAttrMaxSharedMemoryPerMultiprocessor_val, cudaDevAttrMaxSharedMemoryPerMultiprocessor, 0);
  cudaDeviceGetAttribute(&cudaDevAttrMaxRegistersPerMultiprocessor_val, cudaDevAttrMaxRegistersPerMultiprocessor, 0);
  cudaDeviceGetAttribute(&cudaDevAttrMaxBlocksPerMultiprocessor_val, cudaDevAttrMaxBlocksPerMultiprocessor, 0);

  std::cout << "cudaDevAttrMaxThreadsPerBlock_val = " << cudaDevAttrMaxThreadsPerBlock << std::endl;
  std::cout << "cudaDevAttrMaxSharedMemoryPerBlock_val = " << cudaDevAttrMaxSharedMemoryPerBlock << std::endl;
  std::cout << "cudaDevAttrTotalConstantMemory_val = " << cudaDevAttrTotalConstantMemory << std::endl;
  std::cout << "cudaDevAttrClockRate_val = " << cudaDevAttrClockRate << std::endl;

  
  std::cout << "cudaDevAttrMemoryClockRate_val = " << cudaDevAttrMemoryClockRate_val << std::endl;
  std::cout << "cudaDevAttrGlobalMemoryBusWidth_val = " << cudaDevAttrGlobalMemoryBusWidth_val << std::endl;
  std::cout << "cudaDevAttrL2CacheSize_val = " << cudaDevAttrL2CacheSize_val << std::endl;

  std::cout << "cudaDevAttrMaxSharedMemoryPerMultiprocessor_val = " << cudaDevAttrMaxSharedMemoryPerMultiprocessor << std::endl;
  std::cout << "cudaDevAttrMaxRegistersPerMultiprocessor_val = " << cudaDevAttrMaxRegistersPerMultiprocessor << std::endl;
  std::cout << "cudaDevAttrMaxBlocksPerMultiprocessor_val = " << cudaDevAttrMaxBlocksPerMultiprocessor << std::endl;
  
}

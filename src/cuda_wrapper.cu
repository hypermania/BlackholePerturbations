#include "cuda_wrapper.cuh"

template void copy_vector(thrust::device_vector<thrust::complex<double>> &out, const Eigen::PlainObjectBase<Eigen::ArrayXcd> &in);
template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXcd> &out, const thrust::device_vector<thrust::complex<double>> &in);
template void copy_vector(thrust::device_vector<double> &out, const Eigen::PlainObjectBase<Eigen::ArrayXd> &in);
template void copy_vector(Eigen::PlainObjectBase<Eigen::ArrayXd> &out, const thrust::device_vector<double> &in);

template class thrust::device_vector<double>;
template class thrust::device_ptr<double>;
template thrust::device_ptr<double> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<double>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<double>>);
template eigen_iterator thrust::copy(const thrust::detail::execution_policy_base<thrust::cuda_cub::cross_system<thrust::cuda_cub::tag, thrust::system::cpp::detail::tag>> &, thrust_const_iterator, thrust_const_iterator, eigen_iterator);

template thrust_iterator thrust::copy(eigen_iterator, eigen_iterator, thrust_iterator);
template eigen_iterator thrust::copy(thrust_iterator, thrust_iterator, eigen_iterator);
void show_gpu_memory_usage(void)
{
  size_t free, total;
  cudaMemGetInfo(&free, &total);
  std::cout << "free / total = "
	    << free << " B / " << total << " B ("
	    << free / (1024*1024) << " MB / " << total / (1024*1024) << " MB)\n";
}

/*
  New definitions
 */
template class thrust::device_vector<thrust::complex<double>>;
// template class thrust::device_ptr<thrust::complex<double>>;
// template class thrust::device_reference<thrust::complex<double>>;
// template class thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>;
template thrust::device_ptr<thrust::complex<double>> thrust::for_each_n(const thrust::detail::execution_policy_base<thrust::cuda_cub::tag> &, thrust::device_ptr<thrust::complex<double>>, unsigned long, thrust::detail::device_generate_functor<thrust::detail::fill_functor<thrust::complex<double>>>);
template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::strip_const_get_value<thrust::cuda_cub::tag>(thrust::cuda_cub::tag const &) const;
template thrust::complex<double> thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::convert_to_value_type<thrust::cuda_cub::tag>(thrust::cuda_cub::tag *) const;
template thrust::reference<thrust::complex<double>, thrust::device_ptr<thrust::complex<double>>, thrust::device_reference<thrust::complex<double>>>::operator thrust::complex<double>() const;

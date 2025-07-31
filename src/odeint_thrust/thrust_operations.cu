
/*! 
  \file thrust_operations.cu
  \author Siyang Ling
  \brief Auto-generated translation unit for compatibility between thrust library and odeint.
*/
#ifndef THRUST_OPERATIONS_CU
#define THRUST_OPERATIONS_CU
#include "thrust_operations.cuh"

namespace boost {
namespace numeric {
namespace odeint {


__global__
void scale_sum1_kernel(double *v0, const double *v1, const double alpha1, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i];
  }
}


__global__
void scale_sum2_kernel(double *v0, const double *v1, const double *v2, const double alpha1, const double alpha2, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i];
  }
}


__global__
void scale_sum3_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double alpha1, const double alpha2, const double alpha3, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i];
  }
}


__global__
void scale_sum4_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i];
  }
}


__global__
void scale_sum5_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i];
  }
}


__global__
void scale_sum6_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i];
  }
}


__global__
void scale_sum7_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i];
  }
}


__global__
void scale_sum8_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double *v8, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const double alpha8, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i] + alpha8 * v8[i];
  }
}


__global__
void scale_sum9_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double *v8, const double *v9, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const double alpha8, const double alpha9, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i] + alpha8 * v8[i] + alpha9 * v9[i];
  }
}


__global__
void scale_sum10_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double *v8, const double *v9, const double *v10, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const double alpha8, const double alpha9, const double alpha10, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i] + alpha8 * v8[i] + alpha9 * v9[i] + alpha10 * v10[i];
  }
}


__global__
void scale_sum11_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double *v8, const double *v9, const double *v10, const double *v11, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const double alpha8, const double alpha9, const double alpha10, const double alpha11, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i] + alpha8 * v8[i] + alpha9 * v9[i] + alpha10 * v10[i] + alpha11 * v11[i];
  }
}


__global__
void scale_sum12_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double *v8, const double *v9, const double *v10, const double *v11, const double *v12, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const double alpha8, const double alpha9, const double alpha10, const double alpha11, const double alpha12, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i] + alpha8 * v8[i] + alpha9 * v9[i] + alpha10 * v10[i] + alpha11 * v11[i] + alpha12 * v12[i];
  }
}


__global__
void scale_sum13_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double *v8, const double *v9, const double *v10, const double *v11, const double *v12, const double *v13, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const double alpha8, const double alpha9, const double alpha10, const double alpha11, const double alpha12, const double alpha13, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i] + alpha8 * v8[i] + alpha9 * v9[i] + alpha10 * v10[i] + alpha11 * v11[i] + alpha12 * v12[i] + alpha13 * v13[i];
  }
}


__global__
void scale_sum14_kernel(double *v0, const double *v1, const double *v2, const double *v3, const double *v4, const double *v5, const double *v6, const double *v7, const double *v8, const double *v9, const double *v10, const double *v11, const double *v12, const double *v13, const double *v14, const double alpha1, const double alpha2, const double alpha3, const double alpha4, const double alpha5, const double alpha6, const double alpha7, const double alpha8, const double alpha9, const double alpha10, const double alpha11, const double alpha12, const double alpha13, const double alpha14, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v1[i] + alpha2 * v2[i] + alpha3 * v3[i] + alpha4 * v4[i] + alpha5 * v5[i] + alpha6 * v6[i] + alpha7 * v7[i] + alpha8 * v8[i] + alpha9 * v9[i] + alpha10 * v10[i] + alpha11 * v11[i] + alpha12 * v12[i] + alpha13 * v13[i] + alpha14 * v14[i];
  }
}


__global__
void scale_sum_swap_kernel(double *v0, const double *v2, const double alpha1, const double alpha2, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    v0[i] = alpha1 * v0[i] + alpha2 * v2[i];
  }
}

__global__
void rel_error_kernel(double *x_err, const double *x_old, const double *dxdt_old, const double eps_abs, const double eps_rel, const double a_x, const double a_dxdt, const int size)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < size) {
    x_err[i] = fabs(x_err[i]) / (eps_abs + eps_rel * (a_x * fabs(x_old[i]) + a_dxdt * fabs(dxdt_old[i])));
  }
}



  template<typename Fac1>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum1<Fac1>::scale_sum1(const Fac1 alpha1)
    : m_alpha1(alpha1) {}
  
  template<typename Fac1>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum1<Fac1>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum1_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), m_alpha1, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum2<Fac1, Fac2>::scale_sum2(const Fac1 alpha1, const Fac2 alpha2)
    : m_alpha1(alpha1), m_alpha2(alpha2) {}
  
  template<typename Fac1, typename Fac2>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum2<Fac1, Fac2>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum2_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), m_alpha1, m_alpha2, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum3<Fac1, Fac2, Fac3>::scale_sum3(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3) {}
  
  template<typename Fac1, typename Fac2, typename Fac3>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum3<Fac1, Fac2, Fac3>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum3_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), m_alpha1, m_alpha2, m_alpha3, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum4<Fac1, Fac2, Fac3, Fac4>::scale_sum4(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum4<Fac1, Fac2, Fac3, Fac4>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum4_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum5<Fac1, Fac2, Fac3, Fac4, Fac5>::scale_sum5(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum5<Fac1, Fac2, Fac3, Fac4, Fac5>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum5_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum6<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6>::scale_sum6(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum6<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum6_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum7<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7>::scale_sum7(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum7<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum7_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum8<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8>::scale_sum8(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum8<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum8_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v8.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, m_alpha8, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum9<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9>::scale_sum9(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum9<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum9_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v8.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v9.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, m_alpha8, m_alpha9, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum10<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10>::scale_sum10(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum10<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum10_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v8.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v9.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v10.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, m_alpha8, m_alpha9, m_alpha10, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum11<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11>::scale_sum11(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum11<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum11_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v8.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v9.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v10.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v11.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, m_alpha8, m_alpha9, m_alpha10, m_alpha11, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11, typename Fac12>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum12<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11, Fac12>::scale_sum12(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11, const Fac12 alpha12)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11), m_alpha12(alpha12) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11, typename Fac12>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum12<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11, Fac12>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11, const thrust::device_vector<thrust::complex<double>> &v12) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum12_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v8.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v9.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v10.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v11.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v12.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, m_alpha8, m_alpha9, m_alpha10, m_alpha11, m_alpha12, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11, typename Fac12, typename Fac13>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum13<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11, Fac12, Fac13>::scale_sum13(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11, const Fac12 alpha12, const Fac13 alpha13)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11), m_alpha12(alpha12), m_alpha13(alpha13) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11, typename Fac12, typename Fac13>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum13<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11, Fac12, Fac13>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11, const thrust::device_vector<thrust::complex<double>> &v12, const thrust::device_vector<thrust::complex<double>> &v13) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum13_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v8.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v9.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v10.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v11.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v12.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v13.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, m_alpha8, m_alpha9, m_alpha10, m_alpha11, m_alpha12, m_alpha13, 2 * v0.size());
  }


  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11, typename Fac12, typename Fac13, typename Fac14>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum14<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11, Fac12, Fac13, Fac14>::scale_sum14(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11, const Fac12 alpha12, const Fac13 alpha13, const Fac14 alpha14)
    : m_alpha1(alpha1), m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11), m_alpha12(alpha12), m_alpha13(alpha13), m_alpha14(alpha14) {}
  
  template<typename Fac1, typename Fac2, typename Fac3, typename Fac4, typename Fac5, typename Fac6, typename Fac7, typename Fac8, typename Fac9, typename Fac10, typename Fac11, typename Fac12, typename Fac13, typename Fac14>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum14<Fac1, Fac2, Fac3, Fac4, Fac5, Fac6, Fac7, Fac8, Fac9, Fac10, Fac11, Fac12, Fac13, Fac14>::operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11, const thrust::device_vector<thrust::complex<double>> &v12, const thrust::device_vector<thrust::complex<double>> &v13, const thrust::device_vector<thrust::complex<double>> &v14) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    scale_sum14_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v1.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v3.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v4.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v5.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v6.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v7.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v8.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v9.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v10.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v11.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v12.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v13.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v14.data())), m_alpha1, m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6, m_alpha7, m_alpha8, m_alpha9, m_alpha10, m_alpha11, m_alpha12, m_alpha13, m_alpha14, 2 * v0.size());
  }


  template<class Fac1, class Fac2>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum_swap2<Fac1, Fac2>::scale_sum_swap2(const Fac1 alpha1, const Fac2 alpha2)
    : m_alpha1(alpha1), m_alpha2(alpha2) {}

  template<class Fac1, class Fac2>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum_swap2<Fac1, Fac2>::operator()(thrust::device_vector<thrust::complex<double>> &v0, thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (v0.size() + threadsPerBlock - 1) / threadsPerBlock;
    
    v0.swap(v1);
    scale_sum_swap_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(v0.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(v2.data())), m_alpha1, m_alpha2, v0.size());
  }


  template<class Fac>
  thrust_operations<thrust::device_vector<thrust::complex<double>>>::rel_error<Fac>::rel_error(Fac eps_abs, Fac eps_rel, Fac a_x, Fac a_dxdt)
    : m_eps_abs(eps_abs), m_eps_rel(eps_rel), m_a_x(a_x), m_a_dxdt(a_dxdt) {}

  template<class Fac>
  void thrust_operations<thrust::device_vector<thrust::complex<double>>>::rel_error<Fac>::operator()(thrust::device_vector<thrust::complex<double>> &x_err, const thrust::device_vector<thrust::complex<double>> &x_old, const thrust::device_vector<thrust::complex<double>> &dxdt_old) const
  {
    const int threadsPerBlock = 512;
    const int numBlocks = (2 * x_err.size() + threadsPerBlock - 1) / threadsPerBlock;

    rel_error_kernel<<<numBlocks, threadsPerBlock>>>(reinterpret_cast<double *>(thrust::raw_pointer_cast(x_err.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(x_old.data())), reinterpret_cast<const double *>(thrust::raw_pointer_cast(dxdt_old.data())), m_eps_abs, m_eps_rel, m_a_x, m_a_dxdt, 2 * x_err.size());
  }


template class thrust_operations<thrust::device_vector<thrust::complex<double>>>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum1<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum2<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum3<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum4<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum5<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum6<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum7<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum8<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum9<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum10<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum11<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum12<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum13<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum14<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::scale_sum_swap2<double>;
template class thrust_operations<thrust::device_vector<thrust::complex<double>>>::rel_error<double>;

}
}
}

#endif

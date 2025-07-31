
/*! 
  \file thrust_operations.cuh
  \author Siyang Ling
  \brief Auto-generated header for compatibility between thrust library and odeint. Adapated from boost/numeric/odeint/external/thrust/thrust_operations.hpp
*/
#ifndef THRUST_OPERATIONS_CUH
#define THRUST_OPERATIONS_CUH
#include <thrust/device_vector.h>
#include <thrust/complex.h>

namespace boost {
namespace numeric {
namespace odeint {

template<typename Vector>
struct thrust_operations {};


template<>
struct thrust_operations<thrust::device_vector<thrust::complex<double>>> {

  template<typename Fac1 = double>
  struct scale_sum1
  {
    const Fac1 m_alpha1;

    scale_sum1(const Fac1 alpha1);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1) const;
  };


  template<typename Fac1 = double, typename Fac2 = double>
  struct scale_sum2
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2;

    scale_sum2(const Fac1 alpha1, const Fac2 alpha2);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double>
  struct scale_sum3
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3;

    scale_sum3(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double>
  struct scale_sum4
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4;

    scale_sum4(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double>
  struct scale_sum5
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5;

    scale_sum5(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double>
  struct scale_sum6
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6;

    scale_sum6(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double>
  struct scale_sum7
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7;

    scale_sum7(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double, typename Fac8 = double>
  struct scale_sum8
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7; const Fac8 m_alpha8;

    scale_sum8(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double, typename Fac8 = double, typename Fac9 = double>
  struct scale_sum9
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7; const Fac8 m_alpha8; const Fac9 m_alpha9;

    scale_sum9(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double, typename Fac8 = double, typename Fac9 = double, typename Fac10 = double>
  struct scale_sum10
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7; const Fac8 m_alpha8; const Fac9 m_alpha9; const Fac10 m_alpha10;

    scale_sum10(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double, typename Fac8 = double, typename Fac9 = double, typename Fac10 = double, typename Fac11 = double>
  struct scale_sum11
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7; const Fac8 m_alpha8; const Fac9 m_alpha9; const Fac10 m_alpha10; const Fac11 m_alpha11;

    scale_sum11(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double, typename Fac8 = double, typename Fac9 = double, typename Fac10 = double, typename Fac11 = double, typename Fac12 = double>
  struct scale_sum12
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7; const Fac8 m_alpha8; const Fac9 m_alpha9; const Fac10 m_alpha10; const Fac11 m_alpha11; const Fac12 m_alpha12;

    scale_sum12(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11, const Fac12 alpha12);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11, const thrust::device_vector<thrust::complex<double>> &v12) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double, typename Fac8 = double, typename Fac9 = double, typename Fac10 = double, typename Fac11 = double, typename Fac12 = double, typename Fac13 = double>
  struct scale_sum13
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7; const Fac8 m_alpha8; const Fac9 m_alpha9; const Fac10 m_alpha10; const Fac11 m_alpha11; const Fac12 m_alpha12; const Fac13 m_alpha13;

    scale_sum13(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11, const Fac12 alpha12, const Fac13 alpha13);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11, const thrust::device_vector<thrust::complex<double>> &v12, const thrust::device_vector<thrust::complex<double>> &v13) const;
  };


  template<typename Fac1 = double, typename Fac2 = double, typename Fac3 = double, typename Fac4 = double, typename Fac5 = double, typename Fac6 = double, typename Fac7 = double, typename Fac8 = double, typename Fac9 = double, typename Fac10 = double, typename Fac11 = double, typename Fac12 = double, typename Fac13 = double, typename Fac14 = double>
  struct scale_sum14
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2; const Fac3 m_alpha3; const Fac4 m_alpha4; const Fac5 m_alpha5; const Fac6 m_alpha6; const Fac7 m_alpha7; const Fac8 m_alpha8; const Fac9 m_alpha9; const Fac10 m_alpha10; const Fac11 m_alpha11; const Fac12 m_alpha12; const Fac13 m_alpha13; const Fac14 m_alpha14;

    scale_sum14(const Fac1 alpha1, const Fac2 alpha2, const Fac3 alpha3, const Fac4 alpha4, const Fac5 alpha5, const Fac6 alpha6, const Fac7 alpha7, const Fac8 alpha8, const Fac9 alpha9, const Fac10 alpha10, const Fac11 alpha11, const Fac12 alpha12, const Fac13 alpha13, const Fac14 alpha14);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, const thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2, const thrust::device_vector<thrust::complex<double>> &v3, const thrust::device_vector<thrust::complex<double>> &v4, const thrust::device_vector<thrust::complex<double>> &v5, const thrust::device_vector<thrust::complex<double>> &v6, const thrust::device_vector<thrust::complex<double>> &v7, const thrust::device_vector<thrust::complex<double>> &v8, const thrust::device_vector<thrust::complex<double>> &v9, const thrust::device_vector<thrust::complex<double>> &v10, const thrust::device_vector<thrust::complex<double>> &v11, const thrust::device_vector<thrust::complex<double>> &v12, const thrust::device_vector<thrust::complex<double>> &v13, const thrust::device_vector<thrust::complex<double>> &v14) const;
  };


  template<typename Fac1 = double, typename Fac2 = double>
  struct scale_sum_swap2
  {
    const Fac1 m_alpha1; const Fac2 m_alpha2;

    scale_sum_swap2(const Fac1 alpha1, const Fac2 alpha2);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &v0, thrust::device_vector<thrust::complex<double>> &v1, const thrust::device_vector<thrust::complex<double>> &v2) const;
  };


  template<class Fac = double>
  struct rel_error
  {
    const Fac m_eps_abs, m_eps_rel, m_a_x, m_a_dxdt;

    rel_error(Fac eps_abs, Fac eps_rel, Fac a_x, Fac a_dxdt);
    
    void operator()(thrust::device_vector<thrust::complex<double>> &x_err, const thrust::device_vector<thrust::complex<double>> &x_old, const thrust::device_vector<thrust::complex<double>> &dxdt_old) const;
  };

};


}
}
}

#endif

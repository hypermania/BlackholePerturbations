/*
  Modifications by Siyang Ling.
  Adapted so that Eigen can work with runge_kutta_dopri5 and runge_kutta_fehlberg78.
  The generic Runge Kutta algorithms implemented by odeint uses Boost.fusion, which introduces extra memory allocation/free. The "scale_sum*" operations are reimplemented here using native Eigen to avoid this overhead. Exact unit coefficients are handled without a redundant per-element multiplication.
*/

#ifndef EIGEN_OPERATIONS_HPP
#define EIGEN_OPERATIONS_HPP

#include <type_traits>

#include <Eigen/Dense>
#include <boost/multiprecision/float128.hpp>

#if defined(__GNUC__)
// GCC lowers software binary128 arithmetic to these libgcc routines.  The
// explicit calls below preserve the former -ffast-math addition tree while
// removing only multiplication by the exact unit coefficient.
extern "C" float128_type __addtf3(float128_type, float128_type);
extern "C" float128_type __multf3(float128_type, float128_type);
#endif

namespace boost {
  namespace numeric {
    namespace odeint {

      namespace detail {

#if defined(__GNUC__)
      __attribute__((always_inline))
#endif
      inline ::boost::multiprecision::float128 unit_scale_sum6(
          const ::boost::multiprecision::float128 &v1,
          const ::boost::multiprecision::float128 &v2,
          const ::boost::multiprecision::float128 &v3,
          const ::boost::multiprecision::float128 &v4,
          const ::boost::multiprecision::float128 &v5,
          const ::boost::multiprecision::float128 &v6,
          const ::boost::multiprecision::float128 &alpha2,
          const ::boost::multiprecision::float128 &alpha3,
          const ::boost::multiprecision::float128 &alpha4,
          const ::boost::multiprecision::float128 &alpha5,
          const ::boost::multiprecision::float128 &alpha6) {
#if defined(__FAST_MATH__) && defined(__GNUC__)
        using Raw = ::float128_type;
        const Raw sum34 = ::__addtf3(
            ::__multf3(alpha3.backend().value(), v3.backend().value()),
            ::__multf3(alpha4.backend().value(), v4.backend().value()));
        const Raw left = ::__addtf3(v1.backend().value(), sum34);
        const Raw sum56 = ::__addtf3(
            ::__multf3(alpha5.backend().value(), v5.backend().value()),
            ::__multf3(alpha6.backend().value(), v6.backend().value()));
        const Raw right = ::__addtf3(
            ::__multf3(alpha2.backend().value(), v2.backend().value()), sum56);
        return ::boost::multiprecision::float128(::__addtf3(left, right));
#else
        return v1 + alpha2 * v2 + alpha3 * v3 + alpha4 * v4
               + alpha5 * v5 + alpha6 * v6;
#endif
      }

      }  // namespace detail
      
      template<typename D, bool OptimizeUnitCoefficient = true>
      struct eigen_operations {
	// typedef Eigen::internal::traits<D>::Scalar Scalar;
	// typedef Eigen::internal::traits<D>::XprKind XprKind;
	// Eigen::internal::traits<D>::RowsAtCompileTime;
	// Eigen::internal::traits<D>::ColsAtCompileTime;
      };
      
      template<template<typename,int,int> typename EigenObj, typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime, bool OptimizeUnitCoefficient>
      struct eigen_operations<EigenObj<Scalar, RowsAtCompileTime, ColsAtCompileTime>, OptimizeUnitCoefficient> {
	typedef EigenObj<Scalar, RowsAtCompileTime, ColsAtCompileTime> State;
	static constexpr bool EnableUnitOptimization =
	    OptimizeUnitCoefficient
	    && std::is_same_v<Scalar, ::boost::multiprecision::float128>;
	
	template<class Fac = Scalar>
	struct scale_sum1
	{
	  const Fac m_alpha1;

	  scale_sum1(const Fac alpha1)
	    : m_alpha1(alpha1) {}

	  void operator()(State &v0, const State &v1) const
	  {
	    if constexpr(EnableUnitOptimization) if(m_alpha1 == Fac(1)) {
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	      for(Eigen::Index i = 0; i < v0.size(); ++i) v0[i] = v1[i];
	      return;
	    }
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	    for(Eigen::Index i = 0; i < v0.size(); ++i) v0[i] = m_alpha1 * v1[i];
	  }
	};

	template<class Fac = Scalar, class Time1 = Fac>
	struct scale_sum2
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;

	  scale_sum2(const Fac alpha1, const Fac alpha2)
	    : m_alpha1(alpha1) , m_alpha2(alpha2) {}

	  void operator()(State &v0, const State &v1, const State &v2) const
	  {
	    if constexpr(EnableUnitOptimization) if(m_alpha1 == Fac(1)) {
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	      for(Eigen::Index i = 0; i < v0.size(); ++i)
	        v0[i] = v1[i] + m_alpha2 * v2[i];
	      return;
	    }
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	    for(Eigen::Index i = 0; i < v0.size(); ++i)
	      v0[i] = m_alpha1 * v1[i] + m_alpha2 * v2[i];
	  }
	};

	// The extra Time? template parameters are required for dopri5 to work.
	template<class Fac = Scalar, class Time1 = Fac, class Time2 = Fac>
	struct scale_sum3
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;

	  scale_sum3(const Fac alpha1, const Fac alpha2, const Fac alpha3)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3) const
	  {
	    if constexpr(EnableUnitOptimization) if(m_alpha1 == Fac(1)) {
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	      for(Eigen::Index i = 0; i < v0.size(); ++i)
	        v0[i] = v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i];
	      return;
	    }
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	    for(Eigen::Index i = 0; i < v0.size(); ++i)
	      v0[i] = m_alpha1 * v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i];
	  }
	};

	
	//template<class Fac = Scalar, class Time = Scalar>
	template<class Fac = Scalar, class Time1 = Fac, class Time2 = Fac, class Time3 = Fac>
	struct scale_sum4
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;

	  scale_sum4(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4) const
	  {
	    if constexpr(EnableUnitOptimization) if(m_alpha1 == Fac(1)) {
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	      for(Eigen::Index i = 0; i < v0.size(); ++i)
	        v0[i] = v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i]
	                + m_alpha4 * v4[i];
	      return;
	    }
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	    for(Eigen::Index i = 0; i < v0.size(); ++i)
	      v0[i] = m_alpha1 * v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i] + m_alpha4 * v4[i];
	  }
	};


	//template<class Fac = Scalar, class Time = Scalar>
	template<class Fac = Scalar, class Time1 = Fac, class Time2 = Fac, class Time3 = Fac, class Time4 = Fac>
	struct scale_sum5
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;

	  scale_sum5(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5) const
	  {
	    if constexpr(EnableUnitOptimization) if(m_alpha1 == Fac(1)) {
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	      for(Eigen::Index i = 0; i < v0.size(); ++i)
	        v0[i] = v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i]
	                + m_alpha4 * v4[i] + m_alpha5 * v5[i];
	      return;
	    }
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	    for(Eigen::Index i = 0; i < v0.size(); ++i)
	      v0[i] = m_alpha1 * v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i]
	              + m_alpha4 * v4[i] + m_alpha5 * v5[i];
	  }
	};

	
	//template<class Fac = Scalar, class Time = Scalar>
	template<class Fac = Scalar, class Time1 = Fac, class Time2 = Fac, class Time3 = Fac, class Time4 = Fac, class Time5 = Fac>
	struct scale_sum6
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;

	  scale_sum6(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6) const
	  {
	    if constexpr(EnableUnitOptimization) if(m_alpha1 == Fac(1)) {
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	      for(Eigen::Index i = 0; i < v0.size(); ++i)
	        v0[i] = detail::unit_scale_sum6(
	            v1[i], v2[i], v3[i], v4[i], v5[i], v6[i],
	            m_alpha2, m_alpha3, m_alpha4, m_alpha5, m_alpha6);
	      return;
	    }
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	    for(Eigen::Index i = 0; i < v0.size(); ++i)
	      v0[i] = m_alpha1 * v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i]
	              + m_alpha4 * v4[i] + m_alpha5 * v5[i] + m_alpha6 * v6[i];
	  }
	};


	//template<class Fac = Scalar, class Time = Scalar>
	template<class Fac = Scalar, class Time1 = Fac, class Time2 = Fac, class Time3 = Fac, class Time4 = Fac, class Time5 = Fac, class Time6 = Fac>
	struct scale_sum7
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;

	  scale_sum7(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7) const
	  {
#pragma omp parallel for schedule(static) if(v0.size() >= 4096)
	    for(Eigen::Index i = 0; i < v0.size(); ++i)
	      v0[i] = m_alpha1 * v1[i] + m_alpha2 * v2[i] + m_alpha3 * v3[i]
	              + m_alpha4 * v4[i] + m_alpha5 * v5[i] + m_alpha6 * v6[i]
	              + m_alpha7 * v7[i];
	  }
	};


	template<class Fac = Scalar, class Time = Scalar>
	struct scale_sum8
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;
	  const Fac m_alpha8;

	  scale_sum8(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7, const Fac alpha8)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7, const State &v8) const
	  {
	    v0 = m_alpha1 * v1 + m_alpha2 * v2 + m_alpha3 * v3 + m_alpha4 * v4 + m_alpha5 * v5 + m_alpha6 * v6 + m_alpha7 * v7 + m_alpha8 * v8;
	  }
	};


	template<class Fac = Scalar, class Time = Scalar>
	struct scale_sum9
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;
	  const Fac m_alpha8;
	  const Fac m_alpha9;

	  scale_sum9(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7, const Fac alpha8, const Fac alpha9)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7, const State &v8, const State &v9) const
	  {
	    v0 = m_alpha1 * v1 + m_alpha2 * v2 + m_alpha3 * v3 + m_alpha4 * v4 + m_alpha5 * v5 + m_alpha6 * v6 + m_alpha7 * v7 + m_alpha8 * v8 + m_alpha9 * v9;
	  }
	};


	template<class Fac = Scalar, class Time = Scalar>
	struct scale_sum10
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;
	  const Fac m_alpha8;
	  const Fac m_alpha9;
	  const Fac m_alpha10;

	  scale_sum10(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7, const Fac alpha8, const Fac alpha9, const Fac alpha10)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7, const State &v8, const State &v9, const State &v10) const
	  {
	    v0 = m_alpha1 * v1 + m_alpha2 * v2 + m_alpha3 * v3 + m_alpha4 * v4 + m_alpha5 * v5 + m_alpha6 * v6 + m_alpha7 * v7 + m_alpha8 * v8 + m_alpha9 * v9 + m_alpha10 * v10;
	  }
	};


	template<class Fac = Scalar, class Time = Scalar>
	struct scale_sum11
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;
	  const Fac m_alpha8;
	  const Fac m_alpha9;
	  const Fac m_alpha10;
	  const Fac m_alpha11;

	  scale_sum11(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7, const Fac alpha8, const Fac alpha9, const Fac alpha10, const Fac alpha11)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7, const State &v8, const State &v9, const State &v10, const State &v11) const
	  {
	    v0 = m_alpha1 * v1 + m_alpha2 * v2 + m_alpha3 * v3 + m_alpha4 * v4 + m_alpha5 * v5 + m_alpha6 * v6 + m_alpha7 * v7 + m_alpha8 * v8 + m_alpha9 * v9 + m_alpha10 * v10 + m_alpha11 * v11;
	  }
	};


	template<class Fac = Scalar, class Time = Scalar>
	struct scale_sum12
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;
	  const Fac m_alpha8;
	  const Fac m_alpha9;
	  const Fac m_alpha10;
	  const Fac m_alpha11;
	  const Fac m_alpha12;

	  scale_sum12(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7, const Fac alpha8, const Fac alpha9, const Fac alpha10, const Fac alpha11, const Fac alpha12)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11), m_alpha12(alpha12) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7, const State &v8, const State &v9, const State &v10, const State &v11, const State &v12) const
	  {
	    v0 = m_alpha1 * v1 + m_alpha2 * v2 + m_alpha3 * v3 + m_alpha4 * v4 + m_alpha5 * v5 + m_alpha6 * v6 + m_alpha7 * v7 + m_alpha8 * v8 + m_alpha9 * v9 + m_alpha10 * v10 + m_alpha11 * v11 + m_alpha12 * v12;
	  }
	};


	template<class Fac = Scalar, class Time = Scalar>
	struct scale_sum13
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;
	  const Fac m_alpha8;
	  const Fac m_alpha9;
	  const Fac m_alpha10;
	  const Fac m_alpha11;
	  const Fac m_alpha12;
	  const Fac m_alpha13;

	  scale_sum13(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7, const Fac alpha8, const Fac alpha9, const Fac alpha10, const Fac alpha11, const Fac alpha12, const Fac alpha13)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11), m_alpha12(alpha12), m_alpha13(alpha13) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7, const State &v8, const State &v9, const State &v10, const State &v11, const State &v12, const State &v13) const
	  {
	    v0 = m_alpha1 * v1 + m_alpha2 * v2 + m_alpha3 * v3 + m_alpha4 * v4 + m_alpha5 * v5 + m_alpha6 * v6 + m_alpha7 * v7 + m_alpha8 * v8 + m_alpha9 * v9 + m_alpha10 * v10 + m_alpha11 * v11 + m_alpha12 * v12 + m_alpha13 * v13;
	  }
	};


	template<class Fac = Scalar, class Time = Scalar>
	struct scale_sum14
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;
	  const Fac m_alpha3;
	  const Fac m_alpha4;
	  const Fac m_alpha5;
	  const Fac m_alpha6;
	  const Fac m_alpha7;
	  const Fac m_alpha8;
	  const Fac m_alpha9;
	  const Fac m_alpha10;
	  const Fac m_alpha11;
	  const Fac m_alpha12;
	  const Fac m_alpha13;
	  const Fac m_alpha14;

	  scale_sum14(const Fac alpha1, const Fac alpha2, const Fac alpha3, const Fac alpha4, const Fac alpha5, const Fac alpha6, const Fac alpha7, const Fac alpha8, const Fac alpha9, const Fac alpha10, const Fac alpha11, const Fac alpha12, const Fac alpha13, const Fac alpha14)
	    : m_alpha1(alpha1) , m_alpha2(alpha2), m_alpha3(alpha3), m_alpha4(alpha4), m_alpha5(alpha5), m_alpha6(alpha6), m_alpha7(alpha7), m_alpha8(alpha8), m_alpha9(alpha9), m_alpha10(alpha10), m_alpha11(alpha11), m_alpha12(alpha12), m_alpha13(alpha13), m_alpha14(alpha14) {}

	  void operator()(State &v0, const State &v1, const State &v2, const State &v3, const State &v4, const State &v5, const State &v6, const State &v7, const State &v8, const State &v9, const State &v10, const State &v11, const State &v12, const State &v13, const State &v14) const
	  {
	    v0 = m_alpha1 * v1 + m_alpha2 * v2 + m_alpha3 * v3 + m_alpha4 * v4 + m_alpha5 * v5 + m_alpha6 * v6 + m_alpha7 * v7 + m_alpha8 * v8 + m_alpha9 * v9 + m_alpha10 * v10 + m_alpha11 * v11 + m_alpha12 * v12 + m_alpha13 * v13 + m_alpha14 * v14;
	  }
	};


	template<class Fac = Scalar, class Time1 = Fac>
	struct scale_sum_swap2
	{
	  const Fac m_alpha1;
	  const Fac m_alpha2;

	  scale_sum_swap2(const Fac alpha1, const Fac alpha2)
	    : m_alpha1(alpha1) , m_alpha2(alpha2) {}

	  void operator()(State &t1, State &t2, const State &t3) const
	  {
	    t1.swap(t2);
	    t1 = m_alpha1 * t1 + m_alpha2 * t3;
	  }
	};

	
	template<class Fac = Scalar>
	struct rel_error
	{
	  const Fac m_eps_abs, m_eps_rel, m_a_x, m_a_dxdt;

	  rel_error(Fac eps_abs, Fac eps_rel, Fac a_x, Fac a_dxdt)
	    : m_eps_abs(eps_abs), m_eps_rel(eps_rel), m_a_x(a_x), m_a_dxdt(a_dxdt) {}

	  void operator()(State &x_err, const State &x_old, const State &dxdt_old) const
	  {
	    typedef typename Eigen::internal::traits<State>::XprKind XprKind;
	    constexpr bool is_matrix = std::is_same<XprKind, Eigen::MatrixXpr>::value;
	    // Dispatch depending on whether State is Matrix or Array
	    if constexpr(is_matrix) {
	      x_err.array() = x_err.array().abs() / (m_eps_abs + m_eps_rel * (m_a_x * x_old.array().abs() + m_a_dxdt * dxdt_old.array().abs()));
	    } else {
	      x_err = x_err.abs() / (m_eps_abs + m_eps_rel * (m_a_x * x_old.abs() + m_a_dxdt * dxdt_old.abs()));
	    }
	  }
	};
	
      };

      template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
      struct operations_dispatcher<typename Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>>
      {
	typedef eigen_operations<Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>> operations_type;
      };

      template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
      struct operations_dispatcher<typename Eigen::Array<Scalar, RowsAtCompileTime, ColsAtCompileTime>>
      {
	typedef eigen_operations<Eigen::Array<Scalar, RowsAtCompileTime, ColsAtCompileTime>> operations_type;
      };

      
      // Adapted from odeint/external/eigen/eigen_algebra.hpp
      template<typename B,int S1,int S2,int O, int M1, int M2>
      struct vector_space_norm_inf< Eigen::Array<B,S1,S2,O,M1,M2> >
      {
	typedef B result_type;
	result_type operator()( const Eigen::Array<B,S1,S2,O,M1,M2> &m ) const
	{
	  return m.matrix().template lpNorm<Eigen::Infinity>();
	}
      };

      
    } // namespace odeint
  } // namespace numeric
} // namespace boost


namespace Eigen {
  // TODO
  // inline const Eigen::VectorXd
  // operator/(const int &i, const Eigen::VectorXd &m) {
  //   return ((double)i) * m.cwiseInverse();
  // }

  /*
  template<typename D>
  inline const 
  auto
  pow(const Eigen::MatrixBase<D> &m, const double exp) {
    return m.array().pow(exp).matrix();
  }
  */
  
  /*
  template< typename D >
  inline const 
  typename Eigen::CwiseUnaryOp<
    typename Eigen::internal::scalar_abs_op<
      typename Eigen::internal::traits< D >::Scalar > ,
    const D >
  pow(const Eigen::MatrixBase<D> &m, const double exp) {
    return m.array().pow(exp).matrix();
  }
  */
  
  /*
  template<typename D1,typename D2>
  inline const
  typename Eigen::CwiseBinaryOp<
    typename Eigen::internal::scalar_quotient_op<
      typename Eigen::internal::traits<D1>::Scalar>,
    const D1, const D2>
  operator=(Eigen::MatrixBase<D1> &x1, const Eigen::MatrixBase<D2> &x2) {
    return x1.cwiseQuotient(x2);
  }
  */
    


} // end Eigen namespace

#endif

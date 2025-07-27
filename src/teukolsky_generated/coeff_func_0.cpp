
#include "../teukolsky_scalar.hpp"
#include <map>

typedef TeukolskyScalarPDE::HighPrecisionScalar HighPrecisionScalar;


void compute_coeffs_0(const TeukolskyScalarPDE::HighPrecisionScalar a, const TeukolskyScalarPDE::HighPrecisionScalar M, const TeukolskyScalarPDE::HighPrecisionVector &r, const std::vector<TeukolskyScalarPDE::HighPrecisionVector> &expr, std::vector<TeukolskyScalarPDE::ComplexVector> &coeffs) {
coeffs[0] = std::complex<double>(1.0,0.0)*(((-2)*((pow(a,-1))*((pow(r,-3))*(((pow(a,2))+((-1)*((M)*(r))))*((pow((pow(r,3))+((pow(a,2))*(((2)*(M))+(r))),-1))*((expr[0])*((expr[1])*((expr[2])*(expr[3]))))))))).cast<double>());
}

void compute_coeffs_1(const TeukolskyScalarPDE::HighPrecisionScalar a, const TeukolskyScalarPDE::HighPrecisionScalar M, const TeukolskyScalarPDE::HighPrecisionVector &r, const std::vector<TeukolskyScalarPDE::HighPrecisionVector> &expr, std::vector<TeukolskyScalarPDE::ComplexVector> &coeffs) {
coeffs[1] = std::complex<double>(1.0,0.0)*(((2)*((a)*((pow((pow(r,5))+((pow(a,2))*((pow(r,2))*(((2)*(M))+(r)))),-1))*((expr[0])*((expr[1])*((expr[2])*(expr[3]))))))).cast<double>());
}

void compute_coeffs_2(const TeukolskyScalarPDE::HighPrecisionScalar a, const TeukolskyScalarPDE::HighPrecisionScalar M, const TeukolskyScalarPDE::HighPrecisionVector &r, const std::vector<TeukolskyScalarPDE::HighPrecisionVector> &expr, std::vector<TeukolskyScalarPDE::ComplexVector> &coeffs) {
coeffs[2] = std::complex<double>(1.0,0.0)*(((-1)*((pow(a,-1))*((pow(r,-1))*((pow((pow(a,2))+(pow(r,2)),2))*((pow(((pow(r,4))*(((-2)*(M))+(r)))+(((pow(a,4))*(((2)*(M))+(r)))+((pow(a,2))*(((-4)*((pow(M,2))*(r)))+((2)*(pow(r,3)))))),-1))*((expr[0])*((expr[1])*((expr[2])*(expr[3]))))))))).cast<double>());
}

void compute_coeffs_3(const TeukolskyScalarPDE::HighPrecisionScalar a, const TeukolskyScalarPDE::HighPrecisionScalar M, const TeukolskyScalarPDE::HighPrecisionVector &r, const std::vector<TeukolskyScalarPDE::HighPrecisionVector> &expr, std::vector<TeukolskyScalarPDE::ComplexVector> &coeffs) {
coeffs[3] = std::complex<double>(1.0,0.0)*(((HighPrecisionScalar("-0.5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"))*((HighPrecisionScalar("2.236067977499789696409173668731276235440618359611525724270897245410520925637804899414414408378782275"))*((pow(a,-3))*((pow(r,-3))*(((pow(a,2))+((-1)*((M)*(r))))*((pow(((pow(r,4))*(((-2)*(M))+(r)))+(((pow(a,4))*(((2)*(M))+(r)))+((pow(a,2))*(((-4)*((pow(M,2))*(r)))+((2)*(pow(r,3)))))),-1))*((expr[2])*(((3)*((a)*((((pow(r,4))*(((-2)*(M))+(r)))+(((pow(a,4))*(((2)*(M))+(r)))+(((pow(a,2))*(((-4)*((pow(M,2))*(r)))+((2)*(pow(r,3)))))+((pow(expr[0],2))*(pow(expr[1],2))))))*(expr[2]))))+((-2)*(((pow(a,4))+(((3)*(pow(r,4)))+((4)*((pow(a,2))*((r)*((M)+(r)))))))*((expr[0])*((expr[1])*(expr[3]))))))))))))).cast<double>());
}

void compute_coeffs_4(const TeukolskyScalarPDE::HighPrecisionScalar a, const TeukolskyScalarPDE::HighPrecisionScalar M, const TeukolskyScalarPDE::HighPrecisionVector &r, const std::vector<TeukolskyScalarPDE::HighPrecisionVector> &expr, std::vector<TeukolskyScalarPDE::ComplexVector> &coeffs) {
coeffs[4] = std::complex<double>(1.0,0.0)*(((HighPrecisionScalar("2.236067977499789696409173668731276235440618359611525724270897245410520925637804899414414408378782275"))*((pow(r,-2))*((pow(((a)*((pow(r,4))*(((-2)*(M))+(r))))+(((pow(a,5))*(((2)*(M))+(r)))+((pow(a,3))*(((-4)*((pow(M,2))*(r)))+((2)*(pow(r,3)))))),-1))*((expr[0])*((expr[1])*((expr[2])*(((3)*((a)*((expr[0])*((expr[1])*(expr[2])))))+((-1)*(((pow(a,4))+(((3)*(pow(r,4)))+((4)*((pow(a,2))*((r)*((M)+(r)))))))*(expr[3])))))))))).cast<double>());
}

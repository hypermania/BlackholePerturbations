#ifndef IO_HPP
#define IO_HPP
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>



template<typename Scalar>
void write_to_file(const std::vector<Scalar> &vector, std::string filename){
  //char *memblock = (char *)&vector[0];
  const char *memblock = reinterpret_cast<const char *>(vector.data());
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write(memblock, vector.size() * sizeof(Scalar));
  }
}

#ifndef NOT_USING_EIGEN

#include <Eigen/Dense>
template<typename Derived>
void write_to_file(const Eigen::PlainObjectBase<Derived> &obj, std::string filename){
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    // file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::DenseBase<Derived>::Scalar));
    file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar));
  }
}

template<typename Derived>
void write_to_filename_template(const Eigen::PlainObjectBase<Derived> &obj, const std::string format_string, const int idx)
{
  char filename[128];
  sprintf(filename, format_string.data(), idx);
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    // file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::DenseBase<Derived>::Scalar));
    file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar));
  }
}

#endif

// #ifndef DISABLE_CUDA
// #include <thrust/device_vector.h>
// template<typename ThrustScalar>
// void write_to_file(const thrust::device_vector<ThrustScalar> &obj, std::string filename){
//   std::vector<ThrustScalar> buffer(obj.size());
//   copy_vector(buffer, obj);
//   std::ofstream file(filename, std::ios::binary);
//   if(file.is_open()){
//     file.write((char *)buffer.data(), obj.size() * sizeof(ThrustScalar));
//   }
// }
// #endif

#endif

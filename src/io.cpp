#include "io.hpp"

void write_data_to_file(const char *buf, ssize_t size, std::string filename){
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write(buf, size);
  }
}

void write_vector_to_file(std::vector<double> &vector, std::string filename){
  char *memblock = (char *)&vector[0];

  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write(memblock, vector.size() * sizeof(double));
  }
}

std::vector<double> load_vector_from_file(std::string filename){
  std::streampos size = 0;
  char *memblock = NULL;

  std::ifstream file(filename, std::ios::in | std::ios::binary | std::ios::ate);
  if(file.is_open()){
    size = file.tellg();
    memblock = new char[size];
    file.seekg(0, std::ios::beg);
    file.read(memblock, size);
  }
  
  file.close();
  
  //std::cout << "Loading " << filename << ". Size of file is " << size << " bytes.\n";

  double *double_values = (double *)memblock;
  unsigned long long int N = (unsigned long long int)size / sizeof(double);
  
  std::vector<double> v(N);
  for(unsigned long long int i = 0; i < N; i++){
    v[i] = double_values[i];
  }

  return v;
}


boost::math::interpolators::quintic_hermite<std::vector<double>> load_quintic_interpolant_from_files(std::string file_x, std::string file_y, std::string file_dydx, std::string file_d2yd2x){
  std::vector<double> vec_x = load_vector_from_file(file_x);
  std::vector<double> vec_y = load_vector_from_file(file_y);
  std::vector<double> vec_dydx = load_vector_from_file(file_dydx);
  std::vector<double> vec_d2yd2x = load_vector_from_file(file_d2yd2x);
  auto spline = boost::math::interpolators::quintic_hermite(std::move(vec_x), std::move(vec_y), std::move(vec_dydx), std::move(vec_d2yd2x));
  //auto spline = boost::math::interpolators::quintic_hermite(vec_x, vec_y, vec_dydx, vec_d2yd2x);
  return spline;
}


boost::math::interpolators::cubic_hermite<std::vector<double>> load_cubic_interpolant_from_files(std::string file_x, std::string file_y, std::string file_dydx){
  std::vector<double> vec_x = load_vector_from_file(file_x);
  std::vector<double> vec_y = load_vector_from_file(file_y);
  std::vector<double> vec_dydx = load_vector_from_file(file_dydx);
  auto spline = boost::math::interpolators::cubic_hermite(std::move(vec_x), std::move(vec_y), std::move(vec_dydx));
  //auto spline = boost::math::interpolators::cubic_hermite(vec_x, vec_y, vec_dydx);
  return spline;
}


std::vector<double> load_mList_from_file(std::string filename){
  std::vector<double> mList(0);
  std::ifstream mList_file(filename, std::ios::in); 

  if(mList_file.is_open()){
    std::string line;
    while(std::getline(mList_file, line)){
      mList.push_back(std::stod(line));
    }
  }
  mList_file.close();

  return mList;
}


std::vector<std::vector<double>> load_kListList_from_file(std::string filename){
  std::vector<std::vector<double>> kListList(0);
  std::ifstream kListList_file(filename, std::ios::in);

  if(kListList_file.is_open()){
    std::string line;
    while(std::getline(kListList_file, line)){
      std::stringstream lineStream(line);
      std::string cell;
    
      kListList.push_back(std::vector<double>(0));
      while(std::getline(lineStream, cell, ',')){
	kListList.back().push_back(std::stod(cell));
      }
    }
  }
  kListList_file.close();

  return kListList;
}


std::vector<std::pair<double, double>> load_kmList_from_dir(std::string dir){
  std::vector<std::pair<double, double>> kmList(0);
  
  auto mList = load_mList_from_file(dir + "/mList.csv");
  auto kListList = load_kListList_from_file(dir + "/kListList.csv");
  
  for(unsigned int i = 0; i < mList.size(); i++){
    for(unsigned int j = 0; j < kListList[i].size(); j++){
      kmList.push_back(std::make_pair(kListList[i][j], mList[i]));
    }
  }

  return kmList;
}


void write_VectorXd_to_file(const Eigen::VectorXd &vector, std::string filename){
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write((char *)vector.data(), vector.size() * sizeof(double));
  }
}

void write_VectorXd_to_filename_template(const Eigen::VectorXd &vector, const std::string format_string, const int idx)
{
  char filename[128];
  sprintf(filename, format_string.data(), idx);
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write((char *)vector.data(), vector.size() * sizeof(double));
  }
}

Eigen::VectorXd load_VectorXd_from_file(const std::string &filename){
  std::streampos size = 0;

  std::ifstream file(filename, std::ios::in | std::ios::binary | std::ios::ate);
  if(file.is_open()){
    size = file.tellg();
    file.seekg(0, std::ios::beg);
  }

  unsigned long long int N = (unsigned long long int)size / sizeof(double);
  Eigen::VectorXd v(N);
  
  file.read((char *)v.data(), size);
  
  return v;
}


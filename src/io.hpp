#ifndef IO_HPP
#define IO_HPP
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include <boost/math/interpolators/quintic_hermite.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>

#include <Eigen/Dense>

void write_data_to_file(const char *buf, ssize_t size, std::string filename);
void write_vector_to_file(std::vector<double> &vector, std::string filename);
std::vector<double> load_vector_from_file(std::string filename);

void write_VectorXd_to_file(const Eigen::VectorXd &vector, std::string filename);
void write_VectorXd_to_filename_template(const Eigen::VectorXd &vector, const std::string format_string, const int idx);
Eigen::VectorXd load_VectorXd_from_file(const std::string &filename);

boost::math::interpolators::quintic_hermite<std::vector<double>> load_quintic_interpolant_from_files(std::string file_x, std::string file_y, std::string file_dydx, std::string file_d2yd2x);
boost::math::interpolators::cubic_hermite<std::vector<double>> load_cubic_interpolant_from_files(std::string file_x, std::string file_y, std::string file_dydx);

// std::vector<double> load_mList_from_file(std::string filename);
// std::vector<std::vector<double>> load_kListList_from_file(std::string filename);
// std::vector<std::pair<double, double>> load_kmList_from_dir(std::string dir);

#endif

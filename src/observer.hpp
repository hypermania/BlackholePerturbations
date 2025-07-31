/*! 
  \file observer.hpp
  \brief Implements "observers", which controls what gets saved during simulations.
  
*/

#ifndef OBSERVER_HPP
#define OBSERVER_HPP

#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>

#include "Eigen/Dense"


template<typename State>
struct SaveAllObserver {

  std::vector<double> &t;
  std::vector<State> &f;
  std::function<void(const State &, double)> callback;
  
  SaveAllObserver(std::vector<double> &t_, std::vector<State> &f_, std::function<void(const State &, double)> callback_) :
    t(t_), f(f_), callback(callback_)
  {}

  SaveAllObserver(std::vector<double> &t_, std::vector<State> &f_) :
    t(t_), f(f_)
  {}
  
  void operator()(const State &x, double t_current)
  {
    t.push_back(t_current);
    f.push_back(x);
    if(callback){
      callback(x, t_current);
    }
  }
};


struct ConstIntervalObserver {
  
  int idx;
  std::vector<double> t_list;
  std::string dir;
  double t_start; /*!< Start time of simulation. */
  double t_end; /*!< End time of simulation. */
  double t_interval; /*!< Save to file in `dir` every `t_interval`. */
  double t_last;


  template<typename Param>
  ConstIntervalObserver(const std::string &dir_, const Param &param) :
    idx(0), dir(dir_),
    t_start(param.t_start), t_end(param.t_end), t_interval(param.t_interval), t_last(param.t_start) {}

  ConstIntervalObserver(const ConstIntervalObserver &) = default;

  void operator()(const Eigen::VectorXd &x, double t)
  {
    if(t >= t_last + t_interval || t == t_end || t == t_start) {
      write_to_filename_template(x, dir + "state_%d.dat", idx);
      t_list.push_back(t);

      t_last = t;
      ++idx;
    }
  }
};


struct FixedPositionObserver {
  // typedef Eigen::ArrayXd State;
  std::string dir;
  std::vector<long long int> positions;
  std::vector<double> t_list;
  std::vector<double> psi_list;

  FixedPositionObserver(const std::string &dir_, const std::vector<long long int> positions_) :
    dir(dir_), positions(positions_)
  {}
  
  // void operator()(const State &x, const double t) {
  //   for(auto i : positions) {
  //     psi_list.push_back(x[i]);
  //   }
  //   t_list.push_back(t);
  // }

  template<typename State, typename Scalar>
  void operator()(const State &x, const Scalar t) {
    for(auto i : positions) {
      psi_list.push_back(static_cast<double>(x[i]));
    }
    t_list.push_back(static_cast<double>(t));
  }
  
  void save(void) const {
    write_to_file(psi_list, dir + "psi_list.dat");
    write_to_file(t_list, dir + "t_list.dat");    
  }
};


template<typename State = Eigen::ArrayXd, typename TimeScalar = double>
struct GenericFixedPositionObserver {
  typedef Eigen::DenseBase<State>::Scalar StateScalar;
  std::string dir;
  std::vector<long long int> positions;
  std::vector<TimeScalar> t_list;
  std::vector<StateScalar> psi_list;

  GenericFixedPositionObserver(const std::string &dir_, const std::vector<long long int> &positions_) :
    dir(dir_), positions(positions_)
  {}
  
  void operator()(const State &x, const TimeScalar t) {
    for(auto i : positions) {
      psi_list.push_back(x[i]);
    }
    t_list.push_back(t);
  }
  
  void save(void) const {
    write_to_file(psi_list, dir + "psi_list.dat");
    write_to_file(t_list, dir + "t_list.dat");    
  }
};


struct ApproximateTimeObserver {
  // typedef Eigen::ArrayXd State;
  std::string dir;
  std::vector<double> times;
  std::vector<double> t_list;
  size_t current_idx;

  ApproximateTimeObserver(const std::string &dir_, const std::vector<double> times_) :
    dir(dir_), times(times_), current_idx(0)
  {}

  // void operator()(const State &x, const double t) {
  //   if(current_idx < times.size() && t >= times[current_idx]) {
  //     write_to_filename_template(x, dir + "state_%d.dat", current_idx);
  //     ++current_idx;
  //     t_list.push_back(t);
  //   }
  // }

  template<typename State, typename Scalar>
  inline void operator()(const State &x, const Scalar t) {
    if(current_idx < times.size() && t >= times[current_idx]) {
      Eigen::ArrayXd x_temp = x.template cast<double>();
      write_to_filename_template(x_temp, dir + "state_%d.dat", current_idx);
      ++current_idx;
      t_list.push_back(static_cast<double>(t));
    }
  }
  
  void save(void) const {
    write_to_file(t_list, dir + "t_list_snapshots.dat");
  }
};

template<>
inline void ApproximateTimeObserver::operator()<Eigen::ArrayXcd, double>(const Eigen::ArrayXcd &x, const double t) {
  if(current_idx < times.size() && t >= times[current_idx]) {
    //if constexpr(std::is_same_v<Eigen::internal::traits<State>::Scalar, std::complex<double>>){
    Eigen::ArrayXcd x_temp = x;
    write_to_filename_template(x_temp, dir + "state_%d.dat", current_idx);
    ++current_idx;
    t_list.push_back(static_cast<double>(t));
  }
}


template<typename... Observers>
struct ObserverPack {
  //typedef Eigen::ArrayXd State;
  std::tuple<Observers & ...> observers;
  
  ObserverPack(Observers & ... observers_) : observers(observers_...) {}

  template<typename State, typename Scalar>
  //void operator()(const State &x, const double t) {
  void operator()(const State &x, const Scalar t) {
    std::apply([&](auto &&... args) { ((args(x, t)), ...); }, observers);
  }

  void save(void) const {
    std::apply([&](auto &&... args) { ((args.save()), ...); }, observers);
  }
};


#endif

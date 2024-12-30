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



#endif

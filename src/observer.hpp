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

// #include <boost/numeric/odeint.hpp>
// #include <boost/numeric/odeint/external/eigen/eigen.hpp>

template<typename Equation>
struct SaveFDTFObserver {
  typedef typename Equation::State State;
  
  std::vector<double> &t;
  std::vector<double> &f;
  std::vector<double> &dt_f;
  std::function<void(const State &, double)> callback;
  
  SaveFDTFObserver(std::vector<double> &t_, std::vector<double> &f_, std::vector<double> &dt_f_, std::function<void(const State &, double)> callback_) :
    t(t_), f(f_), dt_f(dt_f_), callback(callback_)
  {}

  SaveFDTFObserver(std::vector<double> &t_, std::vector<double> &f_, std::vector<double> &dt_f_) :
    t(t_), f(f_), dt_f(dt_f_)
  {}
  
  void operator()(const State &x, double t_current)
  {
    t.push_back(t_current);
    f.push_back(x[0]);
    dt_f.push_back(x[1]);
    if(callback) {
      callback(x, t_current);
    }
  }
};

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



#endif

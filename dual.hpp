#ifndef INCLUDED_ADOOPP_DUAL
#define INCLUDED_ADOOPP_DUAL

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>

#ifdef SPARSE_DAE
template <typename T>
using DualMap = std::map<int, T>;
#endif

namespace adoopp {

// template <typename T>
class Dual {
 public:
  Dual() : real_(0), dual_(0), index_(0), grad_(1) {}
  Dual(double r) : real_(r), dual_(0), index_(0), grad_(1) {}
  Dual(double r, double d) : real_(r), dual_(d), index_(0), grad_(1) {}
  Dual(const Dual& t)
      : real_(t.real_), dual_(t.dual_), index_(t.index_), grad_(t.grad_) {}
#ifdef SPARSE_DAE
  Dual(const Dual&& t) : real_(t.real_), dual(t.dual_) {
    dual_map_ = std::move(t.dual_map_);
  }
#else
  Dual(const Dual&& t)
      : real_(t.real_), dual_(t.dual_), index_(t.index_), grad_(t.grad_) {}

  Dual(const int& index, std::vector<int> grad)
      : index_(index), grad_(grad), real_(index), dual_(grad[index]) {}
#endif

  Dual& operator=(Dual&& t) = default;
#ifdef SPARSE_DAE
  template <typename T>
  T operator[](unsigned int i) const {
    if (dual_map_.find(i) != dual_map_.end())
      return dual_map_.at(i);
    else
      return 0
  }
#endif

  // Binary Operators
  friend Dual operator+(const Dual& t1, const Dual& t2);
  friend Dual operator-(const Dual& t1, const Dual& t2);
  friend Dual operator*(const Dual& t1, const Dual& t2);
  friend Dual operator/(const Dual& t1, const Dual& t2);

  Dual& operator+=(const Dual& t1);
  Dual& operator-=(const Dual& t1);
  Dual& operator*=(const Dual& t1);
  Dual& operator/=(const Dual& t1);

  // sinh,cosh,tanh

  // There is a clever way to overload these by only using one func and
  // strapping it to each ref I'll have to think about it, it's pretty low
  // priority right now.
  /*
          bool operator==(const Dual& t1, const Dual& t2);
          bool operator!=(const Dual& t1, const Dual& t2);
          bool operator<(const Dual& t1, const Dual& t2);
          bool operator<=(const Dual& t1, const Dual& t2);
          bool operator>(const Dual& t1, const Dual& t2);
          bool operator>=(const Dual& t1, const Dual& t2);
          bool operator==(const Dual& t1, const T& val);
          bool operator!=(const Dual& t1, const T& val);
          bool operator<(const Dual& t1, const T& val);
          bool operator<=(const Dual& t1, const T& val);
          bool operator>(const Dual& t1, const T& val);
          bool operator>=(const Dual& t1, const T& val);
          bool operator==(const T& val, const Dual& t1);
          bool operator!=(const T& val, const Dual& t1);
          bool operator<(const T& val, const Dual& t1);
          bool operator<=(const T& val, const Dual& t1);
          bool operator>(const T& val, const Dual& t1);
          bool operator>=(const T& val, const Dual& t1);
  */
  friend Dual operator+(const Dual& t) { return t; }
  friend Dual operator-(const Dual& t) { return Dual(-t.real_, -t.dual_); }

  friend Dual pow(const Dual& t, double d);
  friend Dual pow(const Dual& t1, const Dual& t2);
  friend Dual sin(const Dual& t);
  friend Dual cos(const Dual& t);
  friend Dual tan(const Dual& t);
  friend Dual sqrt(const Dual& t);
  friend Dual exp(const Dual& t);
  friend Dual log(const Dual& t);
  friend Dual asin(const Dual& t);
  friend Dual acos(const Dual& t);
  friend Dual atan(const Dual& t);
  friend Dual sqr(const Dual& t);
  //   friend Dual Diff(const Dual& t, int d);

  const double& real() const { return real_; }
  const double& dual() const { return dual_; }
  const std::vector<int>& grad() const { return grad_; }
  const int& index() const { return index_; }

  void setReal(double val) { real_ = val; }
  void setDual(double val) { dual_ = val; }

  void setIndex(int i) { index_ = i; }

  //   const bool& size() const { return }

 private:
#ifdef SPARSE_DAE
  bool is_const() const { return dual_map_.size() == 0; }
#else
  bool is_const() const { return dual_ == 0; }
#endif
  bool is_zero() const { return real_ == 0; }
  double real_;
  double dual_;
  int index_;
  std::vector<int> grad_;
#ifdef SPASE_DAE
  DualMap<T> dual_map_;
#endif
};
}  // namespace adoopp

#endif
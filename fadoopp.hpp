#ifndef INCLUDED_FADOOPP_H
#define INCLUDED_FADOOPP_H
#include <vector>
#include "dual.hpp"

typedef void (*ADFcn)(adoopp::Dual* dual);

template <typename T>
using VecT = std::vector<T>;      // template vector
using VecD = VecT<double>;        // double vector
using VecV = VecT<adoopp::Dual>;  // dual vector

namespace adoopp {
class Fadoopp {
  // Constructor: make new vector element and gradient element
 public:
  Fadoopp(const VecD& points, const VecV& duals)
      : points_(points), duals_(duals) {
    VecD jacElem_ = VecD();
  }
  ~Fadoopp() {}

  double grad() {
    grad_ = 0;
    if (jacElems_.size() == 0) {
      computeJacobian();
    }
    for (int i = 0; i < jacElems_.size(); i++) {
      grad_ += jacElems_[i];
    }
    return grad_;
  }
  VecD points() { return points_; }
  VecV duals() { return duals_; }
  void computeJacobian() {
    for (int wrt = 0; wrt < points_.size(); wrt++) {
      for (int i = 0; i < points_.size(); i++)
        if (i == wrt) {
          duals_[i].setDual(1);
        } else {
          duals_[i].setDual(0);
        }
      Dual temp_f;
      temp_f = pow(duals_[0], 2) + duals_[1] * duals_[2];
      jacElems_.push_back(temp_f.dual());
    }
  }
  double jacElem(int wrt) {
    if (jacElems_.size() == 0)
      computeJacobian();
    return jacElems_[wrt];
  }

  void setPoints(const VecD& points) {
    points_.clear();
    points_ = points;
  }

  int size() { return jacElems_.size(); }

 private:
  double grad_;
  VecD points_;
  VecV duals_;
  VecD jacElems_;
};
}  // namespace adoopp
#endif
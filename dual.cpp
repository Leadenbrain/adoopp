#include "dual.hpp"
#include <cassert>
namespace adoopp {
#define AS_OP(OP)                                    \
  Dual operator OP(const Dual& t1, const Dual& t2) { \
    Dual temp(t1);                                   \
    temp.real_ OP## = t2.real_;                      \
    temp.dual_ OP## = t2.dual_;                      \
    return temp;                                     \
  }

AS_OP(+)
AS_OP(-)
/*operator**/
Dual operator*(const Dual& t1, const Dual& t2) {
  Dual temp(t1);
  temp.real_ *= t2.real_;
  temp.dual_ = t1.real_ * t2.dual_ + t1.dual_ * t2.real_;
  return temp;
}
/*operator**/
Dual operator/(const Dual& t1, const Dual& t2) {
  Dual temp(t1);
  temp.real_ /= t2.real_;
  temp.dual_ =
      -1 * (t1.real_ * t2.dual_ - t1.dual_ * t2.real_) / (t2.real_ * t2.real_);
  return temp;
}

Dual& Dual::operator+=(const Dual& t) {
  // both constant
  if (this->is_const() && t.is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ == 0);
    this->real_ += t.real_;
    return *this;
  }
  // one is zero
  if (this->is_zero() || t.is_zero()) {
    assert(this->real_ == 0 || t.real_ == 0);
    this->real_ += t.real_;
    this->dual_ += t.dual_;
    return *this;
  }
  // t non constant
  if (this->is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ != 0);
    this->real_ += t.real_;
    this->dual_ = t.dual_;
    return *this;
  }
  // this non constant
  if (t.is_const()) {
    assert(this->dual_ != 0);
    assert(t.dual_ == 0);
    this->real_ += t.real_;
    return *this;
  }
  this->real_ += t.real_;
  this->dual_ += t.dual_;
  return *this;
}

Dual& Dual::operator-=(const Dual& t) {
  if (this->is_const() && t.is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ == 0);
    this->real_ -= t.real_;
    return *this;
  }
  // one is zero
  if (this->is_zero() || t.is_zero()) {
    assert(this->real_ == 0 || t.real_ == 0);
    this->real_ -= t.real_;
    this->dual_ -= t.dual_;
    return *this;
  }
  // t non constant
  if (this->is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ != 0);
    this->real_ -= t.real_;
    this->dual_ = t.dual_;
    return *this;
  }
  // this non constant
  if (t.is_const()) {
    assert(this->dual_ != 0);
    assert(t.dual_ == 0);
    this->real_ -= t.real_;
    return *this;
  }
  this->real_ -= t.real_;
  this->dual_ -= t.dual_;
  return *this;
}
Dual& Dual::operator*=(const Dual& t) {
  if (this->is_const() && t.is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ == 0);
    this->real_ *= t.real_;
    return *this;
  }
  // one is zero
  if (this->is_zero() || t.is_zero()) {
    assert(this->real_ == 0 || t.real_ == 0);
    this->real_ = 0;
    this->dual_ = this->real_ * t.dual_ + this->dual_ * t.real_;
    return *this;
  }
  // t non constant
  if (this->is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ != 0);
    this->real_ *= t.real_;
    this->dual_ = this->real_ * t.dual_;
    return *this;
  }
  // this non constant
  if (t.is_const()) {
    assert(this->dual_ != 0);
    assert(t.dual_ == 0);
    this->real_ *= t.real_;
    this->dual_ *= t.real_;
    return *this;
  }
  this->real_ *= t.real_;
  this->dual_ = this->real_ * t.dual_ + this->dual_ * t.real_;
  return *this;
}
Dual& Dual::operator/=(const Dual& t) {
  if (this->is_const() && t.is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ == 0);
    this->real_ /= t.real_;
    return *this;
  }
  // one is zero
  if (this->is_zero() || t.is_zero()) {
    assert(this->real_ == 0 || t.real_ == 0);
    this->real_ = 0;
    this->dual_ = -1.0 * (this->real_ * t.dual_ - this->dual_ * t.real_) /
                  (t.real_ * t.real_);
    return *this;
  }
  // t non constant
  if (this->is_const()) {
    assert(this->dual_ == 0);
    assert(t.dual_ != 0);
    this->real_ /= t.real_;
    this->dual_ = t.dual_;
    return *this;
  }
  // this non constant
  if (t.is_const()) {
    assert(this->dual_ != 0);
    assert(t.dual_ == 0);
    this->real_ /= t.real_;
    this->dual_ /= -1.0 * this->dual_ / t.real_;
    return *this;
  }
  this->real_ /= t.real_;
  this->dual_ = -1.0 * (this->real_ * t.dual_ - this->dual_ * t.real_) /
                (t.real_ * t.real_);
  return *this;
}

Dual pow(const Dual& t, double d) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::pow(t.real_, d), 0);
  }
  if (t.is_zero()) {
    assert(t.real_ == 0);
    return Dual(0, d * std::pow(t.real_, d - 1) * t.dual_);
  }
  double real_out = std::pow(t.real_, d);
  double dual_out = d * std::pow(t.real_, d - 1) * t.dual_;
  return Dual(real_out, dual_out);
}

Dual pow(const Dual& t1, const Dual& t2) {
  if (t1.is_const() && t2.is_const()) {
    assert(t1.dual_ == 0);
    assert(t2.dual_ == 0);
    return Dual(std::pow(t1.real_, t2.real_), 0);
  }
  if (t2.is_zero()) {
    // assert(t1.real_ == 1);
    assert(t2.dual_ == 0);
    return Dual(1, 0);
  }
  if (t1.is_zero()) {
    assert(t1.dual_ == 0);
    return Dual(0, 0);
  }
  double real_out = std::pow(t1.real_, t2.real_);
  double dual_out = t2.real_ * std::pow(t1.real_, t2.real_ - 1) * t1.dual_;
  return Dual(real_out, dual_out);
}

/*sin*/
Dual sin(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::sin(t.real_), 0);
  }
  double real_out = std::sin(t.real_);
  double dual_out = std::cos(t.real_) * t.dual_;
  return Dual(real_out, dual_out);
}
/*sin*/
Dual cos(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::cos(t.real_), 0);
  }
  double real_out = std::cos(t.real_);
  double dual_out = -1 * std::sin(t.real_) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual tan(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::tan(t.real_), 0);
  }
  double real_out = std::tan(t.real_);
  double dual_out = (1 + std::tan(t.real_) * std::tan(t.real_)) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual sqrt(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::sqrt(t.real_), 0);
  }
  double real_out = std::sqrt(t.real_);
  double dual_out = (0.5 / std::sqrt(t.real_)) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual exp(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::exp(t.real_), 0);
  }
  double real_out = std::exp(t.real_);
  double dual_out = std::exp(t.real_) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual log(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::log(t.real_), 0);
  }
  double real_out = std::log(t.real_);
  double dual_out = (1.0 / t.real_) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual asin(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::asin(t.real_), 0);
  }
  double real_out = std::asin(t.real_);
  double dual_out = (1.0 / (std::sqrt(1 - t.real_ * t.real_))) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual acos(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::acos(t.real_), 0);
  }
  double real_out = std::acos(t.real_);
  double dual_out = (-1.0 / (std::sqrt(1 - t.real_ * t.real_))) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual atan(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::atan(t.real_), 0);
  }
  double real_out = std::atan(t.real_);
  double dual_out = (1 / (1 + t.real_ * t.real_)) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual sinh(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::sinh(t.real_), 0);
  }
  double real_out = std::sinh(t.real_);
  double dual_out = std::cosh(t.real_) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual cosh(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::cosh(t.real_), 0);
  }
  double real_out = std::cosh(t.real_);
  double dual_out = -1 * std::sinh(t.real_) * t.dual_;
  return Dual(real_out, dual_out);
}
Dual tanh(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(std::tanh(t.real_), 0);
  }
  double real_out = std::tanh(t.real_);
  double dual_out = (1 - (std::tanh(t.real_) * std::tanh(t.real_)) * t.dual_);
  return Dual(real_out, dual_out);
}
Dual sqr(const Dual& t) {
  if (t.is_const()) {
    assert(t.dual_ == 0);
    return Dual(t.real_ * t.real_, 0);
  }
  double real_out = t.real_ * t.real_;
  double dual_out = 2 * t.real_ * t.dual_;
  return Dual(real_out, dual_out);
}
// Dual Diff(const Dual& t, int d) {}
}  // namespace adoopp
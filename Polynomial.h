#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

// Polynomial template class. Requires template arguments satisfying
// std::is_floating_point<T> condition.
// A constructor from the vector of  coefficients is provided, as well as
// comparison operators (equal and not equal), sum, difference and product
// operators. << operator is overloaded for printing polynomials.
// Methods for differentiation, searching for some of the roots are
// implemented. Their descriptions (with time asymptotics estimate) are
// located below.

#include "Utils/Utils.h"

namespace polynomial {

using matrix::matrix_utils::Equal;

template<class T>
concept PolynomialCoefficient =
std::is_floating_point_v<T> || std::is_integral_v<T>;

template<class T> requires PolynomialCoefficient<T>
class Polynomial {
 public:
  // Coefficients vector [a_0, a_1, ..., a_n] represents polynomial
  // a_0 x^{n} + a_1 x^{n - 1} + ... + a_{n - 1} x + a_n.
  explicit Polynomial(std::vector<T> coefficients,
                      T epsilon = std::numeric_limits<T>::epsilon());
  ~Polynomial() = default;

// ---------------------------------------------------------------------------
// Comparison operators. Use predefined epsilons to compare floating point
// values properly.

  bool operator==(const Polynomial& other) const;
  bool operator!=(const Polynomial& other) const;

// ---------------------------------------------------------------------------
// Arithmetic operators.

  Polynomial operator+(const Polynomial& other) const;
  Polynomial operator-(const Polynomial& other) const;
  Polynomial operator*(const Polynomial& other) const;
  Polynomial operator*(T number) const;
  template<class U>
  friend Polynomial<U> operator*(U number, const Polynomial<U>& polynomial);

// ---------------------------------------------------------------------------
// More difficult functions.

  Polynomial Differentiate() const;

// ---------------------------------------------------------------------------
// Printing the polynomial.

  std::string ToString() const;

 private:
  std::vector<T> coefficients_;
  std::vector<double> real_roots_ = {};
  T epsilon_;
};

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T>::Polynomial(std::vector<T> coefficients, T epsilon)
    : coefficients_(std::move(coefficients)), epsilon_(epsilon) {
}

// ---------------------------------------------------------------------------
// Comparison operators.

template<class T>
requires PolynomialCoefficient<T>
bool Polynomial<T>::operator==(const Polynomial& other) const {
  return (coefficients_ == other.coefficients_);
}

template<class T>
requires PolynomialCoefficient<T>
bool Polynomial<T>::operator!=(const Polynomial& other) const {
  return !operator==(other);
}

// ---------------------------------------------------------------------------
// Arithmetic operators.

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T> Polynomial<T>::operator+(const Polynomial& other) const {
  if (other.coefficients_.size() > coefficients_.size()) {
    return other + (*this);
  }

  auto result = *this;
  result.epsilon_ = std::max(epsilon_, other.epsilon_);

  int other_degree = other.coefficients_.size();
  for (int i = 0; i < other_degree; ++i) {
    result.coefficients_[i] += other.coefficients_[i];
  }

  while (result.coefficients_.size() > 1
      && Equal(result.coefficients_.back(), T(), epsilon_)) {
    result.coefficients_.pop_back();
  }

  return result;
}

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T> Polynomial<T>::operator-(const Polynomial& other) const {
  if (other.coefficients_.size() > coefficients_.size()) {
    return other + (*this);
  }

  auto result = *this;
  result.epsilon_ = std::max(epsilon_, other.epsilon_);

  int other_degree = other.coefficients_.size();
  for (int i = 0; i < other_degree; ++i) {
    result.coefficients_[i] -= other.coefficients_[i];
  }

  while (result.coefficients_.size() > 1
      && Equal(result.coefficients_.back(), T(), epsilon_)) {
    result.coefficients_.pop_back();
  }

  return result;
}

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T> Polynomial<T>::operator*(const Polynomial& other) const {
  std::vector<T> result(coefficients_.size() * other.coefficients_.size(), 0);
  int this_degree = coefficients_.size();
  int other_degree = other.coefficients_.size();
  for (int i = 0; i < this_degree; ++i) {
    for (int j = 0; j < other_degree; ++j) {
      result[i + j] += coefficients_[i] * other.coefficients_[j];
    }
  }

  while (result.size() > 1 && Equal(result.back(), T(), epsilon_)) {
    result.pop_back();
  }

  return Polynomial(result, std::max(epsilon_, other.epsilon_));
}

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T> Polynomial<T>::operator*(T number) const {
  return operator*(Polynomial({number}));
}

template<class U>
Polynomial<U> operator*(U number, const Polynomial<U>& polynomial) {
  return polynomial * number;
}

// ---------------------------------------------------------------------------
// More difficult functions.

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T> Polynomial<T>::Differentiate() const {
  if (coefficients_.size() == 1) {
    return Polynomial({0}, epsilon_);
  }
  std::vector<T> result(coefficients_.size() - 1, 0);
  int degree = coefficients_.size();
  for (int i = 1; i < degree; ++i) {
    result[i - 1] = i * coefficients_[i];
  }
  return Polynomial(result, epsilon_);
}

// ---------------------------------------------------------------------------
// Printing the polynomial.

template<class T>
requires PolynomialCoefficient<T>
std::string Polynomial<T>::ToString() const {
  if (coefficients_.size() == 1) {
    return std::to_string(coefficients_[0]);
  }

  std::string result = {};
  int degree = coefficients_.size();
  for (int i = degree - 1; i >= 0; --i) {
    if (Equal(coefficients_[i], T(), epsilon_)) continue;
    if (i != degree - 1) {
      result += " ";
      result += (coefficients_[i] > 0) ? "+" : "-";
      result += " ";
    }
    result += std::to_string(std::abs(coefficients_[i]));
    if (i == 0) continue;
    result += "x";
    if (i == 1) continue;
    result += "^{";
    result += std::to_string(i);
    result += "}";
  }
  return result;
}

template<class U>
std::ostream& operator<<(std::ostream& out, const Polynomial<U>& polynomial) {
  return out << polynomial.ToString();
}

}  // namespace polynomial

#endif  // POLYNOMIAL_H_

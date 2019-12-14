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

#include <algorithm>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Utils/Utils.h"

namespace polynomial {

using matrix::matrix_utils::Equal;

template<class T>
concept PolynomialCoefficient = std::is_floating_point_v<T>;

template<class T> requires PolynomialCoefficient<T>
class Polynomial {
 public:
  // Coefficients vector [a_0, a_1, ..., a_n] represents polynomial
  // a_0 x^{n} + a_1 x^{n - 1} + ... + a_{n - 1} x + a_n.
  explicit Polynomial(std::vector<T> coefficients,
                      T epsilon = std::numeric_limits<T>::epsilon());
  // Default constructor initializes polynomial with 0 value and minimum
  // possible epsilon.
  Polynomial();
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

  // This function returns f(x).
  T GetValue(T x) const;

  // This function returns derivative.
  Polynomial Differentiate() const;

// ---------------------------------------------------------------------------
// Searching for the roots. Bisection and Newton methods are used. Only real
// roots are counted.

 private:
  // This private function runs the bisection method for the passed interval
  // and returns the new interval (which is smaller and fits in epsilon).
  std::tuple<T, T> RunBisectionAlgorithm(T left_border, T right_border,
                                         T epsilon) const;

  // This private function runs the Newton method for the passed interval
  // and returns the root, which is lying in the interval.
  // Sometimes (for bad functions) the results can be unpredictable.
  // If the degree of some root is greater than 1, it can appear im the results
  // any number of times (1, 2, ...).
  std::optional<T> RunNewtonAlgorithm(
      T left_border, T right_border, T epsilon,
      int max_number_of_iterations = 500000) const;

 public:
  // This function counts the roots using firstly the bisection method, and
  // then the Newton method. Brackets are chosen with help of the roots of
  // derivative (so actually to find the roots of n-degree polynomial we need
  // to make recursive calls for (n-1)-degree polynomial, ..., 1-degree
  // polynomial).
  // If some roots have degree greater than 1, they will be counted only once.
  void FindRoots();

  // This function returns the internal vector of roots.
  std::vector<T> GetRoots() const;

// ---------------------------------------------------------------------------
// Printing the polynomial.

  std::string ToString() const;

 private:
  std::vector<T> coefficients_;
  std::vector<T> real_roots_;
  T epsilon_;
};

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T>::Polynomial(std::vector<T> coefficients, T epsilon)
    : coefficients_(std::move(coefficients)), epsilon_(epsilon) {
}

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T>::Polynomial() : Polynomial({0}) {
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
T Polynomial<T>::GetValue(T x) const {
  T value = T();
  T current_degree = 1;
  int degree = coefficients_.size();
  for (int i = 0; i < degree; ++i) {
    value += coefficients_[i] * current_degree;
    current_degree *= x;
  }
  return value;
}

template<class T>
requires PolynomialCoefficient<T>
Polynomial<T> Polynomial<T>::Differentiate() const {
  if (coefficients_.size() == 1) {
    return Polynomial({T()}, epsilon_);
  }
  std::vector<T> result(coefficients_.size() - 1, T());
  int degree = coefficients_.size();
  for (int i = 1; i < degree; ++i) {
    result[i - 1] = i * coefficients_[i];
  }
  return Polynomial(result, epsilon_);
}

// ---------------------------------------------------------------------------
// Searching for the roots.

template<class T>
requires PolynomialCoefficient<T>
std::tuple<T, T> Polynomial<T>::RunBisectionAlgorithm(T left_border,
                                                      T right_border,
                                                      T epsilon) const {
  while (right_border - left_border > epsilon) {
    T central_x = (left_border + right_border) / 2;
    auto function_value_central = GetValue(central_x);
    auto function_value_left = GetValue(left_border);
    auto function_value_right = GetValue(right_border);
    if (function_value_left * function_value_central > 0) {
      left_border = central_x;
    } else {
      right_border = central_x;
    }
  }
  return {left_border, right_border};
}

template<class T>
requires PolynomialCoefficient<T>
std::optional<T> Polynomial<T>::RunNewtonAlgorithm(
    T left_border, T right_border, T epsilon,
    int max_number_of_iterations) const {
  if (GetValue(left_border) * GetValue(right_border) > 0) {
    return std::nullopt;
  }

  auto derivative = Differentiate();
  auto second_derivative = derivative.Differentiate();

  if (Equal(GetValue(left_border), T(), epsilon)) {
    return left_border;
  } else if (Equal(GetValue(right_border), T(), epsilon)) {
    return right_border;
  }

  T previous_point;
  if (GetValue(left_border) * second_derivative.GetValue(left_border) > T()) {
    previous_point = left_border;
  } else if (GetValue(right_border) * second_derivative.GetValue(right_border)
      > T()) {
    previous_point = right_border;
  } else {
    return std::nullopt;
  }

  int number_of_iterations = 0;
  T current_point = previous_point;
  do {
    ++number_of_iterations;
    previous_point = current_point;
    current_point -=
        GetValue(current_point) / derivative.GetValue(current_point);
  } while (number_of_iterations < max_number_of_iterations &&
      std::abs(std::abs(current_point) - std::abs(previous_point)) > epsilon);

  return current_point;
}

template<class T>
requires PolynomialCoefficient<T>
void Polynomial<T>::FindRoots() {
  if (coefficients_.size() == 1) {
    return;
  }
  if (coefficients_.size() == 2) {
    real_roots_.push_back(-coefficients_[0] / coefficients_[1]);
    return;
  }
  if (coefficients_.size() == 3) {
    T D = coefficients_[1] * coefficients_[1]
        - 4 * coefficients_[0] * coefficients_[2];
    if (D < 0) return;
    if (Equal(D, T(), epsilon_)) {
      real_roots_.push_back((-coefficients_[1]) / (2 * coefficients_[2]));
      return;
    }
    real_roots_.push_back(
        (-coefficients_[1] + std::sqrt(D)) / (2 * coefficients_[2]));
    real_roots_.push_back(
        (-coefficients_[1] - std::sqrt(D)) / (2 * coefficients_[2]));
    return;
  }

  auto derivative = Differentiate();
  derivative.FindRoots();
  auto derivative_roots = derivative.GetRoots();

  std::vector<T> brackets = {-1000000000};
  std::copy(derivative_roots.begin(), derivative_roots.end(),
            std::back_inserter(brackets));
  brackets.push_back(1000000000);
  std::sort(brackets.begin(), brackets.end());

  std::vector<std::tuple<T, T>> intervals;
  int number_of_intervals = brackets.size() - 1;
  for (int i = 0; i < number_of_intervals; ++i) {
    intervals.emplace_back(brackets[i], brackets[i + 1]);
  }

  T bisection_epsilon = 0.01;
  for (auto& interval : intervals) {
    interval = RunBisectionAlgorithm(std::get<0>(interval),
                                     std::get<1>(interval),
                                     bisection_epsilon);
  }

  for (const auto& interval : intervals) {
    auto root = RunNewtonAlgorithm(std::get<0>(interval),
                                   std::get<1>(interval),
                                   epsilon_);
    if (!root.has_value()) continue;
    if (!real_roots_.empty()
        && Equal(root.value(), real_roots_.back(), epsilon_)) {
      continue;
    }
    real_roots_.push_back(root.value());
  }
}

template<class T>
requires PolynomialCoefficient<T>
std::vector<T> Polynomial<T>::GetRoots() const {
  return real_roots_;
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
      result += std::to_string(std::abs(coefficients_[i]));
    } else {
      result += std::to_string(coefficients_[i]);
    }
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

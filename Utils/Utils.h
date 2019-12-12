#ifndef UTILS_UTILS_H_
#define UTILS_UTILS_H_

#include <fstream>
#include <iomanip>
#include <random>

namespace matrix::matrix_utils {

extern std::mt19937_64 random_generator;
void SetRandomSeed(unsigned long long seed);

template<class T>
requires std::is_floating_point_v<T> || std::is_integral_v<T>
T Random(int modulo = 5000, bool force_positive = false) {
  T value = static_cast<T>(1) + random_generator() % modulo +
      static_cast<T>(1) / random_generator();
  if (!force_positive && random_generator() % 2 == 1) {
    value = -value;
  }
  return value;
}

template<class T>
requires std::is_floating_point_v<T> || std::is_integral_v<T>
std::vector<std::vector<T>> RandomSquareMatrix(int size) {
  SetRandomSeed(42);
  std::vector<std::vector<T>> result(size, std::vector<T>(size));
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      result[i][j] = Random<T>();
    }
  }
  return result;
}

template<class T>
requires std::is_floating_point_v<T> || std::is_integral_v<T>
inline bool Equal(T a, T b, T epsilon) {
  return (std::abs(a - b) <= epsilon);
}

}  // namespace matrix::matrix_utils

#endif  // UTILS_UTILS_H_

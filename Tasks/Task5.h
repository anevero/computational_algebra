#ifndef TASKS_TASK5_H_
#define TASKS_TASK5_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "../Matrix.h"

namespace matrix::matrix_tasks {

template<class T>
class Task5 {
 public:
  struct Task5Data {
    Matrix<T> matrix;
    int number_of_iterations;
    T first_x;
    T median_x;
    T last_x;
  };

  [[nodiscard]] static Task5Data SolveSystem(int n, T omega = 1.0,
                                             T epsilon = 1e-10);
  [[nodiscard]] static Matrix<T> GetAMatrix(int n);
  [[nodiscard]] static std::vector<std::vector<int>> GetTime(T min_omega,
                                                             T max_omega);
};

template<class T>
typename Task5<T>::Task5Data Task5<T>::SolveSystem(int n, T omega, T epsilon) {
  T old_first_x = 0;
  T old_median_x = 0;
  T old_last_x = 0;
  T new_first_x, new_median_x, new_last_x;
  T first_y, median_y, last_y;

  int max_number_of_iterations =
    std::log(epsilon * n) / std::log(static_cast<double>(n - 1) / n);
  int real_number_of_iterations = 0;

  for (int i = 0; i <= max_number_of_iterations; ++i) {
    new_first_x = (1 - old_last_x - (n - 2) * old_median_x) / n;
    new_median_x =
      (1 - omega * new_first_x - (1 - omega) * old_first_x - old_last_x) / n;
    new_last_x =
      (1 - (n - 2) * (omega * new_median_x + (1 - omega) * old_median_x)
        - omega * new_first_x - (1 - omega) * old_first_x) / n;
    old_first_x = new_first_x;
    old_median_x = new_median_x;
    old_last_x = new_last_x;

    first_y = n * old_first_x + (n - 2) * old_median_x + old_last_x;
    median_y = old_first_x + n * old_median_x + old_last_x;
    last_y = old_first_x + (n - 2) * old_median_x + n * old_last_x;

    if (std::abs(first_y - 1) >= epsilon) continue;
    if (std::abs(median_y - 1) >= epsilon) continue;
    if (std::abs(last_y - 1) >= epsilon) continue;

    real_number_of_iterations = i + 1;
    break;
  }

  std::vector<std::vector<T>> result(n, std::vector<T>(1));
  result[0][0] = old_first_x;
  result[n - 1][0] = old_last_x;
  for (int i = 1; i < n - 1; ++i) {
    result[i][0] = old_median_x;
  }

  return {Matrix(result, epsilon), real_number_of_iterations,
          old_first_x, old_median_x, old_last_x};
}

template<class T>
Matrix<T> Task5<T>::GetAMatrix(int n) {
  std::vector<std::vector<T>> result(n, std::vector<T>(n, 0));
  for (int i = 0; i < n; ++i) {
    result[0][i] = 1;
    result[n - 1][i] = 1;
    result[i][0] = 1;
    result[i][n - 1] = 1;
    result[i][i] = n;
  }
  return Matrix(result);
}

template<class T>
std::vector<std::vector<int>> Task5<T>::GetTime(T min_omega, T max_omega) {
  std::vector<T> omegas(3);
  omegas[0] = min_omega;
  omegas[1] = (min_omega + max_omega) / 2;
  omegas[2] = max_omega;
  std::vector<std::vector<int>> result;
  std::vector<int> current_result;

  for (int i = 0; i < 3; ++i) {
    for (int j = 50; j <= 4000; j += 50) {
      auto x = Task5<T>::SolveSystem(j, omegas[i]);
      current_result.push_back(x.number_of_iterations);
    }
    result.push_back(current_result);
    current_result.clear();
  }

  return result;
}

}  // namespace matrix::matrix_tasks

#endif  // TASKS_TASK5_H_

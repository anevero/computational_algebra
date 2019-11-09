#ifndef TASKS_TASK3_H_
#define TASKS_TASK3_H_

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <random>
#include <stdexcept>
#include <tuple>
#include <typeinfo>
#include <vector>

#include "../Matrix.h"
#include "../Utils/Utils.h"

namespace matrix::matrix_tasks {

template<class T>
class Task3 {
 public:
  [[nodiscard]] static auto RandomInput(int number_of_tests = 3,
                                        int max_matrix_size = 2000);
};

template<class T>
auto Task3<T>::RandomInput(int number_of_tests, int max_matrix_size) {
  static_assert(std::is_floating_point_v<T>);

  std::vector<T> result_tlu{};
  std::vector<T> result_ldl{};

  for (int size = 50; size <= max_matrix_size; size += 50) {
    std::cout << "Current size: " << size << std::endl;

    std::vector<std::vector<T>> a_vector(size, std::vector<T>(size));
    std::vector<std::vector<T>> b_vector(size, std::vector<T>(1));

    double tlu_time = 0;
    double ldl_time = 0;

    for (int i = 0; i < number_of_tests; ++i) {
      // Generating random values for the a matrix.
      for (int j = 0; j < size; ++j) {
        for (int k = 0; k <= j; ++k) {
          T value = matrix_utils::Random();
          a_vector[j][k] = value;
          a_vector[k][j] = value;
        }
      }

      // Generating random values for the b matrix.
      for (int j = 0; j < size; ++j) {
        T value = matrix_utils::Random();
        b_vector[j][0] = value;
      }

      // Creating matrices for TLU and LDL algorithms.
      Matrix a_tlu(a_vector);
      Matrix a_ldl(a_vector);
      Matrix b_tlu(b_vector);
      Matrix b_ldl(b_vector);

      auto t1_tlu = std::chrono::high_resolution_clock::now();
      void(a_tlu.SolveSystem(b_tlu));
      auto t2_tlu = std::chrono::high_resolution_clock::now();

      tlu_time += std::chrono::duration_cast<std::chrono::microseconds>(
        t2_tlu - t1_tlu).count();

      auto t1_ldl = std::chrono::high_resolution_clock::now();
      void(a_ldl.SolveSystem_Symmetric(b_ldl));
      auto t2_ldl = std::chrono::high_resolution_clock::now();

      ldl_time += std::chrono::duration_cast<std::chrono::microseconds>(
        t2_ldl - t1_ldl).count();
    }

    result_tlu.push_back(tlu_time / number_of_tests);
    result_ldl.push_back(ldl_time / number_of_tests);
  }

  return std::tuple{result_tlu, result_ldl};
}

}  // namespace matrix::matrix_tasks

#endif  // TASKS_TASK3_H_

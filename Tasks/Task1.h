#ifndef TASK1_H_
#define TASK1_H_

#include "../Matrix.h"
#include "../Utils/Utils.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <random>
#include <stdexcept>
#include <tuple>
#include <typeinfo>
#include <vector>

template<class T>
class Task1 {
 public:
  struct Task1_Data {
    long double time_common_f;
    long double time_optimized_f;
    long double time_common_d;
    long double time_optimized_d;
    long double time_common_ld;
    long double time_optimized_ld;
  };

  void FirstMatrix();
  void SecondMatrix();
  [[nodiscard]] std::vector<Task1_Data> RandomInput(
      int number_of_tests = 5, int max_matrix_size = 500);
  [[nodiscard]] auto SingleVsMultiThread(
      int number_of_tests = 5, int max_matrix_size = 2000);
  [[nodiscard]] auto TimeForAlmostTriangularMatrix(int matrix_size);
};

template<class T>
void Task1<T>::FirstMatrix() {
  static_assert(std::is_floating_point_v<T>);
  std::cout << "Task 1, first matrix." << std::endl;
  std::cout << "Type: " << typeid(T).name() << std::endl;
  std::cout << "Number of digits: " << std::numeric_limits<T>::digits10
            << std::endl;

  std::fstream in("../input_files/task1_input1");
  int size;
  in >> size;
  std::vector<std::vector<T>> a_vector(size, std::vector<T>(size));
  std::vector<std::vector<T>> e_vector(size, std::vector<T>(size, 0));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      in >> a_vector[i][j];
    }
    e_vector[i][i] = 1;
  }

  in.close();

  Matrix a(a_vector);
  Matrix e(e_vector);

  a.CountInverseMatrix_AlmostTriangular();

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "A^{-1} matrix:" << std::endl << a.GetInverseMatrix()
            << std::endl;

  std::cout << "A * A^{-1} product:" << std::endl << a * a.GetInverseMatrix()
            << std::endl;
  std::cout << "Difference with real E:" << std::endl
            << e - a * a.GetInverseMatrix() << std::endl;

  std::cout << std::endl;
}

template<class T>
void Task1<T>::SecondMatrix() {
  static_assert(std::is_floating_point_v<T>);
  std::cout << "Task 1, second matrix." << std::endl;
  std::cout << "Type: " << typeid(T).name() << std::endl;
  std::cout << "Number of digits: " << std::numeric_limits<T>::digits10
            << std::endl;

  std::fstream in("../input_files/task1_input2");
  int size;
  in >> size;
  std::vector<std::vector<T>> a_vector(size, std::vector<T>(size));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      in >> a_vector[i][j];
    }
  }

  in.close();

  Matrix a(a_vector);
  std::cout << "Trying to count an inverse matrix..." << std::endl;

  try {
    a.CountInverseMatrix_AlmostTriangular();
  } catch (const std::invalid_argument& exc) {
    std::cout << "An exception has been thrown." << std::endl;
    std::cout << "  what(): " << exc.what() << std::endl;
  }

  std::cout << std::endl;
}

template<class T>
std::vector<typename Task1<T>::Task1_Data> Task1<T>::RandomInput(
    int number_of_tests, int max_matrix_size) {
  static_assert(std::is_floating_point_v<T>);

  std::vector<Task1_Data> result{};

  for (int size = 10; size <= max_matrix_size; size += 10) {
    std::cout << "Current size: " << size << std::endl;

    // Float, double and long double matrices.
    std::vector<std::vector<float>>
        f_vector(size, std::vector<float>(size, 0));
    std::vector<std::vector<double>>
        d_vector(size, std::vector<double>(size, 0));
    std::vector<std::vector<long double>>
        ld_vector(size, std::vector<long double>(size, 0));

    std::vector<double> time_common_f;
    std::vector<double> time_optimized_f;
    std::vector<double> time_common_d;
    std::vector<double> time_optimized_d;
    std::vector<double> time_common_ld;
    std::vector<double> time_optimized_ld;

    for (int i = 0; i < number_of_tests; ++i) {
      // Generating random values for the matrix.
      for (int j = 0; j < size; ++j) {
        for (int k = 0; k < j + 2 && k < size; ++k) {
          auto value = Random();
          f_vector[j][k] = value;
          d_vector[j][k] = value;
          ld_vector[j][k] = value;
        }
      }

      // Creating pairs of matrices for common and optimized algorithm.
      Matrix f_common(f_vector);
      Matrix f_optimized(f_vector);
      Matrix d_common(d_vector);
      Matrix d_optimized(d_vector);
      Matrix ld_common(ld_vector);
      Matrix ld_optimized(ld_vector);

      auto t1_common_f = std::chrono::high_resolution_clock::now();
      f_common.CountInverseMatrix();
      auto t2_common_f = std::chrono::high_resolution_clock::now();

      auto t1_optimized_f = std::chrono::high_resolution_clock::now();
      f_optimized.CountInverseMatrix_AlmostTriangular();
      auto t2_optimized_f = std::chrono::high_resolution_clock::now();

      auto t1_common_d = std::chrono::high_resolution_clock::now();
      d_common.CountInverseMatrix();
      auto t2_common_d = std::chrono::high_resolution_clock::now();

      auto t1_optimized_d = std::chrono::high_resolution_clock::now();
      d_optimized.CountInverseMatrix_AlmostTriangular();
      auto t2_optimized_d = std::chrono::high_resolution_clock::now();

      auto t1_common_ld = std::chrono::high_resolution_clock::now();
      ld_common.CountInverseMatrix();
      auto t2_common_ld = std::chrono::high_resolution_clock::now();

      auto t1_optimized_ld = std::chrono::high_resolution_clock::now();
      ld_optimized.CountInverseMatrix_AlmostTriangular();
      auto t2_optimized_ld = std::chrono::high_resolution_clock::now();

      time_common_f.push_back(
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_common_f - t1_common_f).count());

      time_common_d.push_back(
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_common_d - t1_common_d).count());

      time_common_ld.push_back(
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_common_ld - t1_common_ld).count());

      time_optimized_f.push_back(
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_optimized_f - t1_optimized_f).count());

      time_optimized_d.push_back(
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_optimized_d - t1_optimized_d).count());

      time_optimized_ld.push_back(
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_optimized_ld - t1_optimized_ld).count());
    }

    // Sorting arrays to count proper sum later (it helps to overcome precision
    // errors).

    std::sort(time_common_f.begin(), time_common_f.end());
    std::sort(time_optimized_f.begin(), time_optimized_f.end());
    std::sort(time_common_d.begin(), time_common_d.end());
    std::sort(time_optimized_d.begin(), time_optimized_d.end());
    std::sort(time_common_ld.begin(), time_common_ld.end());
    std::sort(time_optimized_ld.begin(), time_optimized_ld.end());

    Task1_Data data;

    data.time_common_f = std::accumulate(time_common_f.begin(),
                                         time_common_f.end(),
                                         0.0L) / number_of_tests;
    data.time_optimized_f = std::accumulate(time_optimized_f.begin(),
                                            time_optimized_f.end(),
                                            0.0L) / number_of_tests;

    data.time_common_d = std::accumulate(time_common_d.begin(),
                                         time_common_d.end(),
                                         0.0L) / number_of_tests;
    data.time_optimized_d = std::accumulate(time_optimized_d.begin(),
                                            time_optimized_d.end(),
                                            0.0L) / number_of_tests;

    data.time_common_ld = std::accumulate(time_common_ld.begin(),
                                          time_common_ld.end(),
                                          0.0L) / number_of_tests;
    data.time_optimized_ld = std::accumulate(time_optimized_ld.begin(),
                                             time_optimized_ld.end(),
                                             0.0L) / number_of_tests;

    result.push_back(data);
  }

  return result;
}

template<class T>
auto Task1<T>::SingleVsMultiThread(int number_of_tests, int max_matrix_size) {
  static_assert(std::is_floating_point_v<T>);

  std::vector<double> result_single_thread{};
  std::vector<double> result_multi_thread{};

  for (int size = 25; size <= max_matrix_size; size += 25) {
    std::cout << "Current size: " << size << std::endl;

    std::vector<std::vector<double>>
        d_vector(size, std::vector<double>(size, 0));

    double time_single_thread = 0;
    double time_multi_thread = 0;

    for (int i = 0; i < number_of_tests; ++i) {
      // Generating random values for the matrix.
      for (int j = 0; j < size; ++j) {
        for (int k = 0; k < j + 2 && k < size; ++k) {
          auto value = Random();
          d_vector[j][k] = value;
        }
      }

      // Creating pairs of matrices for single and multi-threading algorithm.
      Matrix d_single_thread(d_vector);
      Matrix d_multi_thread(d_vector);

      auto t1_single_thread = std::chrono::high_resolution_clock::now();
      d_single_thread.CountInverseMatrix_AlmostTriangular_SingleThread();
      auto t2_single_thread = std::chrono::high_resolution_clock::now();

      auto t1_multi_thread = std::chrono::high_resolution_clock::now();
      d_multi_thread.CountInverseMatrix_AlmostTriangular();
      auto t2_multi_thread = std::chrono::high_resolution_clock::now();

      time_single_thread +=
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_single_thread - t1_single_thread).count();

      time_multi_thread +=
          std::chrono::duration_cast<std::chrono::microseconds>(
              t2_multi_thread - t1_multi_thread).count();
    }

    result_single_thread.push_back(time_single_thread / number_of_tests);
    result_multi_thread.push_back(time_multi_thread / number_of_tests);
  }

  return std::tuple(result_single_thread, result_multi_thread);
}

template<class T>
auto Task1<T>::TimeForAlmostTriangularMatrix(int matrix_size) {
  static_assert(std::is_floating_point_v<T>);

  std::vector<std::vector<T>>
      m_vector(matrix_size, std::vector<T>(matrix_size));

  for (int j = 0; j < matrix_size; ++j) {
    for (int k = 0; k < j + 2 && k < matrix_size; ++k) {
      m_vector[j][k] = Random();
    }
  }

  Matrix<T> m(m_vector);

  auto t1 = std::chrono::high_resolution_clock::now();
  m.CountInverseMatrix_AlmostTriangular();
  auto t2 = std::chrono::high_resolution_clock::now();

  return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
}

#endif  // TASK1_H_

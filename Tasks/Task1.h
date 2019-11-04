#ifndef TASKS_TASK1_H_
#define TASKS_TASK1_H_

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <tuple>
#include <typeinfo>
#include <vector>

#include "../Matrix.h"
#include "../Utils/Utils.h"

template<class T>
class Task1 {
 public:
  struct Task1_Data {
    double time_common_f;
    double time_tlu_f;
    double time_triangular_f;

    double time_common_d;
    double time_tlu_d;
    double time_triangular_d;

    double time_common_ld;
    double time_tlu_ld;
    double time_triangular_ld;
  };

  static void FirstMatrix();
  static void SecondMatrix();
  [[nodiscard]] static std::vector<Task1_Data> RandomInput(
      int max_matrix_size = 3000);
  [[nodiscard]] static auto SingleVsMultiThread(int max_matrix_size = 5000);
  [[nodiscard]] static auto TimeForAlmostTriangularMatrix(int matrix_size);
};

template<class T>
void Task1<T>::FirstMatrix() {
  static_assert(std::is_floating_point_v<T>);
  std::cout << "Task 1, first matrix." << std::endl;
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
  } catch (const std::runtime_error& exc) {
    std::cout << "An exception has been thrown." << std::endl;
    std::cout << "  what(): " << exc.what() << std::endl;
  }

  std::cout << std::endl;
}

template<class T>
std::vector<typename Task1<T>::Task1_Data> Task1<T>::RandomInput(
    int max_matrix_size) {
  static_assert(std::is_floating_point_v<T>);

  std::vector<Task1_Data> result{};

  for (int size = 150; size <= max_matrix_size; size += 150) {
    std::cout << "Current size: " << size << std::endl;

    // Float, double and long double matrices.
    std::vector<std::vector<float>>
        f_vector(size, std::vector<float>(size));
    std::vector<std::vector<double>>
        d_vector(size, std::vector<double>(size));
    std::vector<std::vector<long double>>
        ld_vector(size, std::vector<long double>(size));

    std::vector<double> time_common_f;
    std::vector<double> time_tlu_f;
    std::vector<double> time_triangular_f;
    std::vector<double> time_common_d;
    std::vector<double> time_tlu_d;
    std::vector<double> time_triangular_d;
    std::vector<double> time_common_ld;
    std::vector<double> time_tlu_ld;
    std::vector<double> time_triangular_ld;

    // Generating random values for the matrix.
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < j + 2 && k < size; ++k) {
        auto value = Random();
        f_vector[j][k] = value;
        d_vector[j][k] = value;
        ld_vector[j][k] = value;
      }
    }

    Matrix f_common(f_vector);
    Matrix f_tlu(f_vector);
    Matrix f_triangular(f_vector);
    Matrix d_common(f_vector);
    Matrix d_tlu(f_vector);
    Matrix d_triangular(f_vector);
    Matrix ld_common(f_vector);
    Matrix ld_tlu(f_vector);
    Matrix ld_triangular(f_vector);

    auto t1_common_f = std::chrono::high_resolution_clock::now();
    f_common.CountInverseMatrix();
    auto t2_common_f = std::chrono::high_resolution_clock::now();

    auto t1_tlu_f = std::chrono::high_resolution_clock::now();
    f_tlu.CountInverseMatrix_AlmostTriangular_Tlu();
    auto t2_tlu_f = std::chrono::high_resolution_clock::now();

    auto t1_triangular_f = std::chrono::high_resolution_clock::now();
    f_triangular.CountInverseMatrix_AlmostTriangular();
    auto t2_triangular_f = std::chrono::high_resolution_clock::now();

    auto t1_common_d = std::chrono::high_resolution_clock::now();
    d_common.CountInverseMatrix();
    auto t2_common_d = std::chrono::high_resolution_clock::now();

    auto t1_tlu_d = std::chrono::high_resolution_clock::now();
    d_tlu.CountInverseMatrix_AlmostTriangular_Tlu();
    auto t2_tlu_d = std::chrono::high_resolution_clock::now();

    auto t1_triangular_d = std::chrono::high_resolution_clock::now();
    d_triangular.CountInverseMatrix_AlmostTriangular();
    auto t2_triangular_d = std::chrono::high_resolution_clock::now();

    auto t1_common_ld = std::chrono::high_resolution_clock::now();
    ld_common.CountInverseMatrix();
    auto t2_common_ld = std::chrono::high_resolution_clock::now();

    auto t1_tlu_ld = std::chrono::high_resolution_clock::now();
    ld_tlu.CountInverseMatrix_AlmostTriangular_Tlu();
    auto t2_tlu_ld = std::chrono::high_resolution_clock::now();

    auto t1_triangular_ld = std::chrono::high_resolution_clock::now();
    ld_triangular.CountInverseMatrix_AlmostTriangular();
    auto t2_triangular_ld = std::chrono::high_resolution_clock::now();

    Task1_Data data;

    data.time_common_f = std::chrono::duration_cast<std::chrono::microseconds>(
        t2_common_f - t1_common_f).count();
    data.time_tlu_f = std::chrono::duration_cast<std::chrono::microseconds>(
        t2_tlu_f - t1_tlu_f).count();
    data.time_triangular_f =
        std::chrono::duration_cast<std::chrono::microseconds>(
            t2_triangular_f - t1_triangular_f).count();

    data.time_common_d = std::chrono::duration_cast<std::chrono::microseconds>(
        t2_common_d - t1_common_d).count();
    data.time_tlu_d = std::chrono::duration_cast<std::chrono::microseconds>(
        t2_tlu_d - t1_tlu_d).count();
    data.time_triangular_d =
        std::chrono::duration_cast<std::chrono::microseconds>(
            t2_triangular_d - t1_triangular_d).count();

    data.time_common_ld = std::chrono::duration_cast<std::chrono::microseconds>(
        t2_common_ld - t1_common_ld).count();
    data.time_tlu_ld = std::chrono::duration_cast<std::chrono::microseconds>(
        t2_tlu_ld - t1_tlu_ld).count();
    data.time_triangular_ld =
        std::chrono::duration_cast<std::chrono::microseconds>(
            t2_triangular_ld - t1_triangular_ld).count();

    result.push_back(data);
  }

  return result;
}

template<class T>
auto Task1<T>::SingleVsMultiThread(int max_matrix_size) {
  static_assert(std::is_floating_point_v<T>);

  std::vector<double> result_tlu_single_thread{};
  std::vector<double> result_tlu_multi_thread{};
  std::vector<double> result_triangular_single_thread{};
  std::vector<double> result_triangular_multi_thread{};

  for (int size = 500; size <= max_matrix_size; size += 500) {
    std::cout << "Current size: " << size << std::endl;

    std::vector<std::vector<double>>
        d_vector(size, std::vector<double>(size, 0));

    // Generating random values for the matrix.
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < j + 2 && k < size; ++k) {
        auto value = Random();
        d_vector[j][k] = value;
      }
    }

    // Creating pairs of matrices for single and multi-threading algorithm.
    Matrix d_single_thread_tlu(d_vector);
    Matrix d_multi_thread_tlu(d_vector);
    Matrix d_single_thread_triangular(d_vector);
    Matrix d_multi_thread_triangular(d_vector);

    auto t1_single_thread_tlu = std::chrono::high_resolution_clock::now();
    d_single_thread_tlu.CountInverseMatrix_AlmostTriangular_Tlu_SingleThread();
    auto t2_single_thread_tlu = std::chrono::high_resolution_clock::now();

    auto t1_multi_thread_tlu = std::chrono::high_resolution_clock::now();
    d_multi_thread_tlu.CountInverseMatrix_AlmostTriangular_Tlu();
    auto t2_multi_thread_tlu = std::chrono::high_resolution_clock::now();

    auto
        t1_single_thread_triangular = std::chrono::high_resolution_clock::now();
    d_single_thread_triangular.CountInverseMatrix_AlmostTriangular_SingleThread();
    auto
        t2_single_thread_triangular = std::chrono::high_resolution_clock::now();

    auto t1_multi_thread_triangular = std::chrono::high_resolution_clock::now();
    d_multi_thread_triangular.CountInverseMatrix_AlmostTriangular();
    auto t2_multi_thread_triangular = std::chrono::high_resolution_clock::now();

    auto time_single_thread_tlu =
        std::chrono::duration_cast<std::chrono::microseconds>(
            t2_single_thread_tlu - t1_single_thread_tlu).count();
    auto time_multi_thread_tlu =
        std::chrono::duration_cast<std::chrono::microseconds>(
            t2_multi_thread_tlu - t1_multi_thread_tlu).count();

    auto time_single_thread_triangular =
        std::chrono::duration_cast<std::chrono::microseconds>(
            t2_single_thread_triangular - t1_single_thread_triangular).count();
    auto time_multi_thread_triangular =
        std::chrono::duration_cast<std::chrono::microseconds>(
            t2_multi_thread_triangular - t1_multi_thread_triangular).count();

    result_tlu_single_thread.push_back(time_single_thread_tlu);
    result_tlu_multi_thread.push_back(time_multi_thread_tlu);
    result_triangular_single_thread.push_back(time_single_thread_triangular);
    result_triangular_multi_thread.push_back(time_multi_thread_triangular);
  }

  return std::tuple{result_tlu_single_thread,
                    result_tlu_multi_thread,
                    result_triangular_single_thread,
                    result_triangular_multi_thread};
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

#endif  // TASKS_TASK1_H_

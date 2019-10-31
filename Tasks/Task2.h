#ifndef TASKS_TASK2_H_
#define TASKS_TASK2_H_

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <tuple>
#include <typeinfo>
#include <vector>

#include "../Matrix.h"
#include "../Utils/Utils.h"

template<class T>
class Task2 {
 public:
  void FirstMatrix();
  void SecondMatrix();
  [[nodiscard]] auto RandomBForSecondMatrix(int number_of_tests = 250);
};

template<class T>
void Task2<T>::FirstMatrix() {
  static_assert(std::is_floating_point_v<T>);
  std::cout << "Task 2, first matrix." << std::endl;
  std::cout << "Type: " << typeid(T).name() << std::endl;
  std::cout << "Number of digits: " << std::numeric_limits<T>::digits10
            << std::endl;

  std::fstream in("../input_files/task2_input1");
  int size;
  in >> size;
  std::vector<std::vector<T>> a_vector(size, std::vector<T>(size));
  std::vector<std::vector<T>> b_vector(size, std::vector<T>(1));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      in >> a_vector[i][j];
    }
  }
  for (int i = 0; i < size; ++i) {
    in >> b_vector[i][0];
  }

  in.close();

  Matrix a(a_vector);
  Matrix b(b_vector);

  a.CountTluDecomposition();
  auto x = a.SolveSystem(b);
  a.CountConditionNumber();
  auto condition_number = a.GetConditionNumber().value();

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "A matrix condition number: " << condition_number << std::endl;
  std::cout << "B matrix:" << std::endl << b << std::endl;
  std::cout << "X matrix:" << std::endl << x << std::endl;

  std::cout << "Condition number is not too big (but also not small enough)."
            << std::endl << "The result should be mostly right." << std::endl;

  std::cout << "AX product:" << std::endl << a * x << std::endl;
  std::cout << "Difference with real B:" << std::endl << b - a * x << std::endl;

  std::cout << std::endl;
}

template<class T>
void Task2<T>::SecondMatrix() {
  static_assert(std::is_floating_point_v<T>);
  std::cout << "Task 2, second matrix." << std::endl;
  std::cout << "Type: " << typeid(T).name() << std::endl;
  std::cout << "Number of digits: " << std::numeric_limits<T>::digits10
            << std::endl;

  std::fstream in("../input_files/task2_input2");
  int size;
  in >> size;
  std::vector<std::vector<T>> a_vector(size, std::vector<T>(size));
  std::vector<std::vector<T>> b_vector(size, std::vector<T>(1));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      in >> a_vector[i][j];
    }
  }
  for (int i = 0; i < size; ++i) {
    in >> b_vector[i][0];
  }

  in.close();

  Matrix a(a_vector);
  Matrix b(b_vector);

  a.CountTluDecomposition();
  auto x = a.SolveSystem(b);
  a.CountConditionNumber();
  auto condition_number = a.GetConditionNumber().value();

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "A matrix condition number: " << condition_number << std::endl;
  std::cout << "B matrix:" << std::endl << b << std::endl;
  std::cout << "X matrix:" << std::endl << x << std::endl;

  std::cout << "Condition number is too big." << std::endl
            << "Results can be far from right." << std::endl;

  std::cout << "AX product:" << std::endl << a * x << std::endl;
  std::cout << "Difference with real B:" << std::endl << b - a * x << std::endl;

  std::cout << std::endl;
}

template<class T>
auto Task2<T>::RandomBForSecondMatrix(int number_of_tests) {
  static_assert(std::is_floating_point_v<T>);

  std::fstream in("../input_files/task2_input2");
  int size;
  in >> size;

  std::vector<std::vector<float>> a_vector_f(size, std::vector<float>(size));
  std::vector<std::vector<double>> a_vector_d(size, std::vector<double>(size));
  std::vector<std::vector<long double>>
      a_vector_ld(size, std::vector<long double>(size));

  std::vector<std::vector<float>> b_vector_f(size, std::vector<float>(1));
  std::vector<std::vector<double>> b_vector_d(size, std::vector<double>(1));
  std::vector<std::vector<long double>>
      b_vector_ld(size, std::vector<long double>(1));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      in >> a_vector_ld[i][j];
      a_vector_f[i][j] = a_vector_ld[i][j];
      a_vector_d[i][j] = a_vector_ld[i][j];
    }
  }

  in.close();

  Matrix a_f(a_vector_f);
  Matrix a_d(a_vector_d);
  Matrix a_ld(a_vector_ld);

  a_f.CountTluDecomposition();
  a_d.CountTluDecomposition();
  a_ld.CountTluDecomposition();

  std::vector<T> result_f{};
  std::vector<T> result_d{};
  std::vector<T> result_ld{};

  int seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937_64 random_generator(seed);

  for (int i = 0; i < number_of_tests; ++i) {
    std::cout << "Test number: " << i + 1 << std::endl;
    for (int j = 0; j < size; ++j) {
      b_vector_ld[j][0] = Random();
      b_vector_d[j][0] = b_vector_ld[j][0];
      b_vector_f[j][0] = b_vector_ld[j][0];
    }

    Matrix b_f(b_vector_f);
    Matrix b_d(b_vector_d);
    Matrix b_ld(b_vector_ld);

    auto x_f = a_f.SolveSystem(b_f);
    auto ax_f = a_f * x_f;
    auto diff_f = b_f - ax_f;

    auto x_d = a_d.SolveSystem(b_d);
    auto ax_d = a_d * x_d;
    auto diff_d = b_d - ax_d;

    auto x_ld = a_ld.SolveSystem(b_ld);
    auto ax_ld = a_ld * x_ld;
    auto diff_ld = b_ld - ax_ld;

    diff_f.CountNorm();
    diff_d.CountNorm();
    diff_ld.CountNorm();

    result_f.push_back(diff_f.GetNorm().value());
    result_d.push_back(diff_d.GetNorm().value());
    result_ld.push_back(diff_ld.GetNorm().value());
  }

  std::sort(result_f.begin(), result_f.end());
  std::sort(result_d.begin(), result_d.end());
  std::sort(result_ld.begin(), result_ld.end());

  return std::tuple{result_f, result_d, result_ld};
}

#endif  // TASKS_TASK2_H_

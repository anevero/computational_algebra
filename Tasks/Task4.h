#ifndef TASKS_TASK4_H_
#define TASKS_TASK4_H_

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <tuple>
#include <typeinfo>
#include <vector>

#include "../Matrix.h"
#include "../Utils/Utils.h"

namespace matrix::matrix_tasks {

template<class T>
class Task4 {
 public:
  static void FirstMatrix();
  static void SecondMatrix();
};

template<class T>
void Task4<T>::FirstMatrix() {
  static_assert(std::is_floating_point_v<T>);
  std::cout << "Task 4, first matrix." << std::endl;
  std::fstream in("../input_files/task4_input1");
  int size;
  in >> size;
  std::vector<std::vector<T>> a_vector(size, std::vector<T>(4));
  std::vector<std::vector<T>> b_vector(size, std::vector<T>(1));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < 3; ++j) {
      in >> a_vector[i][j];
    }
    a_vector[i][3] = 0;
  }

  for (int i = 0; i < size; ++i) {
    in >> b_vector[i][0];
  }

  in.close();

  Matrix a(a_vector);
  Matrix b(b_vector);

  auto sol = a.SolveSystem_Tridiagonal(b);

  std::cout << "A matrix:" << std::endl << a.GetTridiagonalMatrixAsNormal()
            << std::endl;
  std::cout << "b matrix:" << std::endl << b << std::endl;
  std::cout << "x matrix:" << std::endl << sol << std::endl;

  std::cout << "A * x product:" << std::endl
            << a.GetTridiagonalMatrixAsNormal() * sol << std::endl;
  std::cout << "Difference with real b:" << std::endl
            << b - a.GetTridiagonalMatrixAsNormal() * sol << std::endl;

  std::cout << std::endl;
}

template<class T>
void Task4<T>::SecondMatrix() {
  static_assert(std::is_floating_point_v<T>);
  std::cout << "Task 4, second matrix." << std::endl;
  std::fstream in("../input_files/task4_input2");
  int size;
  in >> size;
  std::vector<std::vector<T>> a_vector(size, std::vector<T>(4));
  std::vector<std::vector<T>> b_vector(size, std::vector<T>(1));

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < 3; ++j) {
      in >> a_vector[i][j];
    }
    a_vector[i][3] = 0;
  }

  for (int i = 0; i < size; ++i) {
    in >> b_vector[i][0];
  }

  in.close();

  Matrix a(a_vector);
  Matrix b(b_vector);

  auto sol = a.SolveSystem_Tridiagonal(b);

  std::cout << "A matrix:" << std::endl << a.GetTridiagonalMatrixAsNormal()
            << std::endl;
  std::cout << "b matrix:" << std::endl << b << std::endl;
  std::cout << "x matrix:" << std::endl << sol << std::endl;

  std::cout << "A * x product:" << std::endl
            << a.GetTridiagonalMatrixAsNormal() * sol << std::endl;
  std::cout << "Difference with real b:" << std::endl
            << b - a.GetTridiagonalMatrixAsNormal() * sol << std::endl;

  std::cout << std::endl;
}

}  // namespace matrix::matrix_tasks

#endif  // TASKS_TASK4_H_

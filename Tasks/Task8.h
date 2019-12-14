#ifndef TASKS_TASK8_H_
#define TASKS_TASK8_H_

#include <limits>
#include <vector>

#include "../Matrix.h"
#include "../Utils/Utils.h"

namespace matrix::matrix_tasks {

template<class T> requires std::is_floating_point_v<T>
class Task8 {
 public:
  static void FirstMatrix();
  static void SecondMatrix();
  [[nodiscard]] static std::vector<int> RandomInput(
      int max_matrix_size = 100, T epsilon = 0.0000000001,
      int max_number_of_iterations = 100000);
  [[nodiscard]] static int TimeForMatrix(
      int matrix_size, T epsilon = 0.0000000001,
      int max_number_of_iterations = 100000);
};

template<class T>
requires std::is_floating_point_v<T>
void Task8<T>::FirstMatrix() {
  std::cout << "Task 8, first matrix." << std::endl;
  std::fstream in("../input_files/task6_input1");
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
  a.RunQrAlgorithm(0.000001);
  auto result = a.GetEigenvectorsFromHessenbergMatrix(0.000001);

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "Eigenvalues:" << std::endl;
  std::cout << "(" << result.size() << " items)" << std::endl;

  for (const auto& item : result) {
    std::cout << item.value << std::endl;
  }

  std::cout << std::endl;
}

template<class T>
requires std::is_floating_point_v<T>
void Task8<T>::SecondMatrix() {
  std::cout << "Task 8, second matrix." << std::endl;
  std::fstream in("../input_files/task6_input2");
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
  a.RunQrAlgorithm(0.000001);
  auto result = a.GetEigenvectorsFromHessenbergMatrix(0.000001);

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "Eigenvalues:" << std::endl;
  std::cout << "(" << result.size() << " items)" << std::endl;

  for (const auto& item : result) {
    std::cout << item.value << std::endl;
  }

  std::cout << std::endl;
}

template<class T>
requires std::is_floating_point_v<T>
std::vector<int> Task8<T>::RandomInput(int max_matrix_size, T epsilon,
                                       int max_number_of_iterations) {
  matrix_utils::SetRandomSeed(42);

  std::vector<int> result;

  for (int size = 5; size <= max_matrix_size; size += 5) {
    std::cout << "Current size: " << size << std::endl;
    int time = 0;
    int number_of_tests = 3;

    for (int j = 0; j < number_of_tests; ++j) {
      auto diagonal_matrix_vector = matrix_utils::RandomDiagonalMatrix<T>(size);
      auto square_matrix_vector = matrix_utils::RandomSquareMatrix<T>(size);

      auto square_matrix = Matrix(square_matrix_vector);
      square_matrix.CountInverseMatrix();
      auto matrix = square_matrix * Matrix(diagonal_matrix_vector)
          * square_matrix.GetInverseMatrix();

      auto t1 = std::chrono::high_resolution_clock::now();
      matrix.RunQrAlgorithm(epsilon, max_number_of_iterations);
      auto t2 = std::chrono::high_resolution_clock::now();
      time += std::chrono::duration_cast<std::chrono::microseconds>(
          t2 - t1).count();
    }

    result.push_back(time / number_of_tests);
  }

  return result;
}

template<class T>
requires std::is_floating_point_v<T>
int Task8<T>::TimeForMatrix(int size, T epsilon, int max_number_of_iterations) {
  matrix_utils::SetRandomSeed(42);

  auto diagonal_matrix_vector = matrix_utils::RandomDiagonalMatrix<T>(size);
  auto square_matrix_vector = matrix_utils::RandomSquareMatrix<T>(size);

  auto square_matrix = Matrix(square_matrix_vector);
  square_matrix.CountInverseMatrix();
  auto matrix = square_matrix * Matrix(diagonal_matrix_vector)
      * square_matrix.GetInverseMatrix();

  auto t1 = std::chrono::high_resolution_clock::now();
  matrix.RunQrAlgorithm(epsilon, max_number_of_iterations);
  auto t2 = std::chrono::high_resolution_clock::now();

  return std::chrono::duration_cast<std::chrono::microseconds>(
      t2 - t1).count();
}

}  // namespace matrix::matrix_tasks

#endif  // TASKS_TASK8_H_

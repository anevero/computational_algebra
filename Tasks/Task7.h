#ifndef TASKS_TASK7_H_
#define TASKS_TASK7_H_

#include <limits>
#include <vector>

#include "../Matrix.h"
#include "../Utils/Utils.h"

namespace matrix::matrix_tasks {

template<class T> requires std::is_floating_point_v<T>
class Task7 {
 public:
  static void FirstMatrix();
  static void SecondMatrix();
};

template<class T>
requires std::is_floating_point_v<T>
void Task7<T>::FirstMatrix() {
  std::cout << "Task 7, first matrix." << std::endl;
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
  a.CountFrobeniusMatrix();
  a.CountCharacteristicPolynomial();
  a.FindCharacteristicPolynomialRoots();

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "Characteristic polynomial:" << std::endl
            << a.GetCharacteristicPolynomial() << std::endl;

  std::cout << "Characteristic polynomial real roots:" << std::endl;
  for (const auto& item : a.GetEigenvectorsFromFrobeniusMatrix()) {
    std::cout << item.value << std::endl;
  }

  std::cout << "Characteristic polynomial real roots (eigenvalues) and "
               "corresponding eigenvectors:" << std::endl;

  for (const auto& item : a.GetEigenvectorsFromFrobeniusMatrix()) {
    std::cout << item.value << std::endl << item.vector << std::endl;
  }

  std::cout << std::endl;
}

template<class T>
requires std::is_floating_point_v<T>
void Task7<T>::SecondMatrix() {
  std::cout << "Task 7, second matrix." << std::endl;
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
  a.CountFrobeniusMatrix();
  a.CountCharacteristicPolynomial();
  a.FindCharacteristicPolynomialRoots();

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "Characteristic polynomial:" << std::endl
            << a.GetCharacteristicPolynomial() << std::endl;

  std::cout << "Characteristic polynomial real roots:" << std::endl;
  for (const auto& item : a.GetEigenvectorsFromFrobeniusMatrix()) {
    std::cout << item.value << std::endl;
  }

  std::cout << "Characteristic polynomial real roots (eigenvalues) and "
               "corresponding eigenvectors:" << std::endl;

  for (const auto& item : a.GetEigenvectorsFromFrobeniusMatrix()) {
    std::cout << item.value << std::endl << item.vector << std::endl;
  }

  std::cout << std::endl;
}

}  // namespace matrix::matrix_tasks

#endif  // TASKS_TASK7_H_

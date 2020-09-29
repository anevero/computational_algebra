#ifndef TASKS_TASK6_H_
#define TASKS_TASK6_H_

#include <limits>
#include <vector>

#include "../matrix.h"
#include "../utils/utils.h"

namespace matrix::matrix_tasks {

template<class T> requires std::is_floating_point_v<T>
class Task6 {
 public:
  static void FirstMatrix();
  static void SecondMatrix();
};

template<class T>
requires std::is_floating_point_v<T>
void Task6<T>::FirstMatrix() {
  std::cout << "Task 6, first matrix." << std::endl;
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
  auto result =
      a.GetPowerIterationResults(std::numeric_limits<T>::epsilon(), 2000);

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "Real maximum module eigenvalues and "
               "corresponding eigenvectors:" << std::endl;

  for (const auto& item : result) {
    std::cout << item.value << std::endl << item.vector << std::endl;
  }

  if (result.empty()) {
    std::cout << "- Maximum eigenvalue is not real." << std::endl;
  }

  std::cout << std::endl;
}

template<class T>
requires std::is_floating_point_v<T>
void Task6<T>::SecondMatrix() {
  std::cout << "Task 6, second matrix." << std::endl;
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
  auto result =
      a.GetPowerIterationResults(std::numeric_limits<double>::epsilon(), 2000);

  std::cout << "A matrix:" << std::endl << a << std::endl;
  std::cout << "Real maximum module eigenvalues and "
               "corresponding eigenvectors:" << std::endl;

  for (const auto& item : result) {
    std::cout << item.value << std::endl << item.vector << std::endl;
  }

  if (result.empty()) {
    std::cout << "- Maximum eigenvalue is not real." << std::endl;
  }

  std::cout << std::endl;
}

}  // namespace matrix::matrix_tasks

#endif  // TASKS_TASK6_H_

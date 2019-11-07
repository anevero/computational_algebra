#include <iomanip>
#include <iostream>

#include "Tests/Tests.h"
#include "Utils/DataOutput.h"

int main() {
  matrix::matrix_tests::TestEverything();
  matrix::matrix_utils::DataOutput::Task1_RandomInput();
  matrix::matrix_utils::DataOutput::Task1_SingleVsMultiThread();
  matrix::matrix_utils::DataOutput::Task2_RandomBForSecondMatrix();
  matrix::matrix_utils::DataOutput::Task3_RandomInput();
  matrix::matrix_utils::DataOutput::Task5_DifferentOmegas();

  int matrix_size;
  while (true) {
    std::cout << "Enter matrix size:" << std::endl;
    std::cin >> matrix_size;
    if (matrix_size <= 0) {
      break;
    }
    std::cout << "Time:" << std::endl;
    auto time_us = static_cast<double>(
        matrix::matrix_tasks::Task1<double>::TimeForAlmostTriangularMatrix(
            matrix_size));
    std::cout << time_us / 1000000 << "seconds" << std::endl;
  }

  return 0;
}

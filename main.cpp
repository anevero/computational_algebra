#include <iomanip>
#include <iostream>

#include "Tests/Tests.h"
#include "Utils/DataOutput.h"

using namespace matrix;

int main() {
//  matrix::matrix_tests::TestEverything();
//  matrix::matrix_utils::DataOutput::Task1_RandomInput();
//  matrix::matrix_utils::DataOutput::Task1_SingleVsMultiThread();
//  matrix::matrix_utils::DataOutput::Task2_RandomBForSecondMatrix();
//  matrix::matrix_utils::DataOutput::Task3_RandomInput();
//  matrix::matrix_utils::DataOutput::Task5_DifferentOmegas();

//  matrix_tests::TestQr();

  std::vector<std::vector<double>> v = {{5, 3, -18, -3},
                                        {1, -82.7, 101.7, 81.7},
                                        {4, -146.2, 168.2, 146.2},
                                        {-8, -8, 0, 7}};
  Matrix m(v);
  for (int i = 0; i < 100; ++i) {
    m.RunQRAlgorithmIteration();
    std::cout << m.GetUpperHessenbergMatrix() << std::endl;
  }

//  int matrix_size;
//  while (true) {
//    std::cout << "Enter matrix size:" << std::endl;
//    std::cin >> matrix_size;
//    if (matrix_size <= 0) {
//      break;
//    }
//    std::cout << "Time:" << std::endl;
//    auto time_us = static_cast<double>(
//      matrix::matrix_tasks::Task1<double>::TimeForAlmostTriangularMatrix(
//        matrix_size));
//    std::cout << time_us / 1000000 << " seconds" << std::endl;
//  }

  return 0;
}

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

//  std::vector<std::vector<double>> v = {{5, 3, -18, -3},
//                                        {1, -82.7, 101.7, 81.7},
//                                        {4, -146.2, 168.2, 146.2},
//                                        {-8, -8, 0, 7}};
  std::vector<std::vector<double>>
      v = {{90, 72, -36},
           {-91, -73, 36},
           {34, 34, -18}};

  Matrix m(v);
  auto t = m.GetPowerIterationResults(0.000001);
  for (auto item : t) {
    std::cout << std::get<0>(item) << std::endl << std::get<1>(item)
              << std::endl;
  }

  return 0;
}

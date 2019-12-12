#include <iomanip>
#include <iostream>

#include "Tests/Tests.h"
#include "Utils/DataOutput.h"

using namespace matrix;
using namespace polynomial;

int main() {
//  matrix::matrix_tests::TestEverything();

//  std::vector<std::vector<double>> v = {{5, 3, -18, -3},
//                                        {1, -82.7, 101.7, 81.7},
//                                        {4, -146.2, 168.2, 146.2},
//                                        {-8, -8, 0, 7}};

//  std::vector<std::vector<double>> v = {{90, 72, -36},
//                                        {-91, -73, 36},
//                                        {34, 34, -18}};

//  std::vector<std::vector<double>> v = {{7, -12, -2},
//                                        {3, -4, 0},
//                                        {-2, 0, -2}};

//  std::vector<std::vector<double>> v = {{-2, 2, 4},
//                                        {2, 3, -1},
//                                        {4, -1, 1}};

//  std::vector<std::vector<double>> v = {{2, 2, 4},
//                                        {2, 3, -1},
//                                        {4, -1, 1}};

//  Matrix m(v);
//  m.CountFrobeniusMatrix();
//  std::cout << m.GetFrobeniusMatrix() << std::endl;

//  std::fstream in("../input_files/task6_input2");
//  double a;
//  in >> a;
//  std::vector<std::vector<double>> v(20, std::vector<double>(20));
//  for (int i = 0; i < 20; ++i) {
//    for (int j = 0; j < 20; ++j) {
//      in >> v[i][j];
//    }
//  }

//  Matrix m(v);
//  m.RunQrAlgorithm(0.0000001);
//  auto t = m.GetEigenvectorsFromHessenbergMatrix(0.001);
//  for (const auto& item : t) {
//    std::cout << item.value << std::endl << item.vector << std::endl;
//  }

//  Matrix m(v);
//  std::cout << m << std::endl;
//  auto t = m.GetPowerIterationResults(0.000001);
//  for (auto& item : t) {
//    std::cout << item.value << std::endl << item.vector << std::endl;
//  }

  Polynomial<int> p({1, 0, 1, 2});
  Polynomial<int> q({-1, 0, 1});

  std::cout << 5 * p * q;


  return 0;
}

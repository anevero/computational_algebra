#include <iomanip>
#include <iostream>

#include "Tests/Tests.h"
#include "Utils/DataOutput.h"

using namespace matrix;
using namespace polynomial;

int main() {
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

//  std::vector<std::vector<double>> v = {{1, 5, 10, 15},
//                                        {0, 54, 45, 32},
//                                        {0, 42, 42, 1},
//                                        {0, 0, 0, 4}};

//  std::vector<std::vector<double>> v = {{25, 20, -10},
//                                        {-26, -21, 10},
//                                        {8, 8, -5}};

  std::fstream in("../input_files/task6_input1");
  double a;
  in >> a;
  std::vector<std::vector<double>> v(a, std::vector<double>(a));
  for (int i = 0; i < a; ++i) {
    for (int j = 0; j < a; ++j) {
      in >> v[i][j];
    }
  }

  Matrix m(v);
  m.CountFrobeniusMatrix();
  m.CountCharacteristicPolynomial();
  m.FindCharacteristicPolynomialRoots();
  std::cout << m.GetFrobeniusMatrix() << std::endl;

  for (const auto& item : m.GetEigenvectorsFromFrobeniusMatrix()) {
    std::cout << item.value << std::endl << item.vector << std::endl;
  }

//  Matrix m(v);
//  m.RunQrAlgorithm(0.00000001);
//  auto t = m.GetEigenvectorsFromHessenbergMatrix(0.01);
//  for (const auto& item : t) {
//    std::cout << item.value << std::endl << item.vector << std::endl;
//  }
//  std::cout << m.GetUpperHessenbergMatrix();

//  Matrix m(v);
//  std::cout << m << std::endl;
//  auto t = m.GetPowerIterationResults(0.000001);
//  for (auto& item : t) {
//    std::cout << item.value << std::endl << item.vector << std::endl;
//  }

  return 0;
}

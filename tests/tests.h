#ifndef TESTS_TESTS_H_
#define TESTS_TESTS_H_

#include <chrono>
#include <random>
#include <vector>

#include "../matrix.h"
#include "../tasks/task5.h"
#include "../utils/utils.h"

// These tests should be run manually. Some of them can fail (it's expected
// behaviour), so you'll need to check the failure message.

namespace matrix::matrix_tests {

void TestEverything();

void TestTlu(int number_of_tests = 500, int matrix_size = 100);
void TestSystemSolution(int number_of_tests = 500, int matrix_size = 100);
void TestInverse(int number_of_tests = 500, int matrix_size = 100);

void TestTlu_AlmostTriangular(int number_of_tests = 500, int matrix_size = 100);
void TestSystemSolution_AlmostTriangular(int number_of_tests = 500,
                                         int matrix_size = 100);
void TestInverse_AlmostTriangular(int number_of_tests = 500,
                                  int matrix_size = 100);
void TestInverse_TLU_AlmostTriangular(int number_of_tests = 500,
                                      int matrix_size = 100);

void TestLdl(int number_of_tests = 500, int matrix_size = 100);
void TestSystemSolution_Symmetric(int number_of_tests = 500,
                                  int matrix_size = 100);

void TestSystemSolution_Tridiagonal(int number_of_tests = 1000,
                                    int matrix_size = 200);

void TestSystemSolution_Sor(int number_of_tests = 200);

void TestQr(int number_of_tests = 500, int matrix_size = 200);

}  // namespace matrix::matrix_tests

#endif  // TESTS_TESTS_H_

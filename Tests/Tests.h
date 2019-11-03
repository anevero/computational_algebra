#ifndef TESTS_TESTS_H_
#define TESTS_TESTS_H_

#include <vector>

#include "../Matrix.h"
#include "../Utils/Utils.h"

// These tests should be run manually. Some of them can fail (it's expected
// behaviour), so you'll need to check the failure message.

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

#endif  // TESTS_TESTS_H_

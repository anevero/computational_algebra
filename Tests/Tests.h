#ifndef TESTS_H_
#define TESTS_H_

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

void TestLdl(int number_of_tests = 500, int matrix_size = 100);
void TestSystemSolution_Symmetric(int number_of_tests = 500,
                                  int matrix_size = 100);

#endif  // TESTS_H_

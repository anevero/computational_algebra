#include "Tests.h"

void TestEverything() {
  TestTlu();
  TestSystemSolution();
  TestInverse();

  TestTlu_AlmostTriangular();
  TestSystemSolution_AlmostTriangular();
  TestInverse_TLU_AlmostTriangular();
  TestInverse_AlmostTriangular();

  TestLdl();
  TestSystemSolution_Symmetric();

  TestSystemSolution_Tridiagonal();

  TestSystemSolution_Sor();
}

void TestTlu(int number_of_tests, int matrix_size) {
  std::vector<std::vector<double>>
      v(matrix_size, std::vector<double>(matrix_size));

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < matrix_size; ++j) {
        v[i][j] = Random();
      }
    }

    Matrix<double> m(v, 0.0001);
    m.CountTluDecomposition();
    auto mdlu =
        m.GetTInverseMatrix_TLU() * m.GetLMatrix_TLU() * m.GetUMatrix_TLU();

    if (m != mdlu) {
      std::cout << "TestTlu: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << m << std::endl;
      std::cout << "D^{-1} L U:" << std::endl << mdlu << std::endl;
      return;
    }
  }

  std::cout << "TestTlu: completed." << std::endl;
}

void TestSystemSolution(int number_of_tests, int matrix_size) {
  std::vector<std::vector<double>>
      v1(matrix_size, std::vector<double>(matrix_size));
  std::vector<std::vector<double>>
      v2(matrix_size, std::vector<double>(1));

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < matrix_size; ++j) {
        v1[i][j] = Random();
      }
      v2[i][0] = Random();
    }

    Matrix<double> a(v1, 0.0001);
    Matrix<double> b(v2, 0.0001);
    auto sol = a.SolveSystem(b);

    if (a * sol != b) {
      std::cout << "TestSystemSolution: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << a << std::endl;
      std::cout << "Solution:" << std::endl << sol << std::endl;
      std::cout << "B:" << std::endl << b << std::endl;
      std::cout << "Corresponding B:" << std::endl << a * sol << std::endl;
      std::cout << "B - AX:" << std::endl << b - a * sol << std::endl;
      return;
    }
  }

  std::cout << "TestSystemSolution: completed." << std::endl;
}

void TestInverse(int number_of_tests, int matrix_size) {
  std::vector<std::vector<double>>
      v1(matrix_size, std::vector<double>(matrix_size));

  std::vector<std::vector<double>>
      v2(matrix_size, std::vector<double>(matrix_size, 0));
  for (int i = 0; i < matrix_size; ++i) {
    v2[i][i] = 1;
  }

  Matrix<double> id(v2);

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < matrix_size; ++j) {
        v1[i][j] = Random();
      }
    }

    Matrix<double> m(v1, 0.0001);
    m.CountInverseMatrix();
    auto inverse = m.GetInverseMatrix();

    if (m * inverse != id) {
      std::cout << "TestInverse: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << m << std::endl;
      std::cout << "Inverse:" << std::endl << inverse << std::endl;
      std::cout << "Product:" << std::endl << m * inverse << std::endl;
      return;
    }
  }

  std::cout << "TestInverse: completed." << std::endl;
}

void TestTlu_AlmostTriangular(int number_of_tests, int matrix_size) {
  std::vector<std::vector<long double>>
      v(matrix_size, std::vector<long double>(matrix_size, 0));

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < i + 2 && j < matrix_size; ++j) {
        v[i][j] = Random();
      }
    }

    Matrix<long double> m(v, 0.00001);
    m.CountTluDecomposition_AlmostTriangular();
    auto mdlu =
        m.GetTInverseMatrix_TLU() * m.GetLMatrix_TLU() * m.GetUMatrix_TLU();

    if (m != mdlu) {
      std::cout << "TestTlu_AlmostTriangular: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << m << std::endl;
      std::cout << "D^{-1}LU:" << std::endl << mdlu << std::endl;
      return;
    }
  }

  std::cout << "TestTlu_AlmostTriangular: completed." << std::endl;
}

void TestSystemSolution_AlmostTriangular(int number_of_tests, int matrix_size) {
  std::vector<std::vector<long double>>
      v1(matrix_size, std::vector<long double>(matrix_size, 0));
  std::vector<std::vector<long double>>
      v2(matrix_size, std::vector<long double>(1));

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < i + 2 && j < matrix_size; ++j) {
        v1[i][j] = Random();
      }
      v2[i][0] = Random();
    }

    Matrix<long double> a(v1, 0.00001);
    Matrix<long double> b(v2, 0.00001);
    auto sol = a.SolveSystem_AlmostTriangular(b);

    if (a * sol != b) {
      std::cout << "TestSystemSolution_AlmostTriangular: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << a << std::endl;
      std::cout << "Solution:" << std::endl << sol << std::endl;
      std::cout << "B:" << std::endl << b << std::endl;
      std::cout << "Corresponding B:" << std::endl << a * sol << std::endl;
      std::cout << "B - AX:" << std::endl << b - a * sol << std::endl;
      return;
    }
  }

  std::cout << "TestSystemSolution_AlmostTriangular: completed." << std::endl;
}

void TestInverse_TLU_AlmostTriangular(int number_of_tests, int matrix_size) {
  std::vector<std::vector<long double>>
      v1(matrix_size, std::vector<long double>(matrix_size, 0));

  std::vector<std::vector<long double>>
      v2(matrix_size, std::vector<long double>(matrix_size, 0));
  for (int i = 0; i < matrix_size; ++i) {
    v2[i][i] = 1;
  }

  Matrix<long double> id(v2, 0.000001);

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < i + 2 && j < matrix_size; ++j) {
        v1[i][j] = Random();
      }
    }

    Matrix<long double> m(v1, 0.000001);
    m.CountInverseMatrix_AlmostTriangular_Tlu();
    auto inverse = m.GetInverseMatrix();

    if (m * inverse != id) {
      std::cout << "TestInverse_TLU_AlmostTriangular: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << m << std::endl;
      std::cout << "Inverse:" << std::endl << inverse << std::endl;
      std::cout << "Product:" << std::endl << m * inverse << std::endl;
      return;
    }
  }

  std::cout << "TestInverse_TLU_AlmostTriangular: completed." << std::endl;
}

void TestInverse_AlmostTriangular(int number_of_tests, int matrix_size) {
  std::vector<std::vector<double>>
      v1(matrix_size, std::vector<double>(matrix_size, 0));

  std::vector<std::vector<double>>
      v2(matrix_size, std::vector<double>(matrix_size, 0));
  for (int i = 0; i < matrix_size; ++i) {
    v2[i][i] = 1;
  }

  Matrix<double> id(v2, 0.0001);

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < i + 2 && j < matrix_size; ++j) {
        v1[i][j] = Random();
      }
    }

    Matrix<double> m(v1, 0.0001);
    m.CountInverseMatrix_AlmostTriangular();
    auto inverse = m.GetInverseMatrix();

    if (m * inverse != id) {
      std::cout << "TestInverse_AlmostTriangular: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << m << std::endl;
      std::cout << "Inverse:" << std::endl << inverse << std::endl;
      std::cout << "Product:" << std::endl << m * inverse << std::endl;
      return;
    }
  }

  std::cout << "TestInverse_AlmostTriangular: completed." << std::endl;
}

void TestLdl(int number_of_tests, int matrix_size) {
  std::vector<std::vector<double>>
      v(matrix_size, std::vector<double>(matrix_size));

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j <= i; ++j) {
        v[i][j] = Random();
        v[j][i] = v[i][j];
      }
    }

    Matrix<double> m(v, 0.0001);
    m.CountLdlDecomposition_Symmetric();

    auto mldl = m.GetLTMatrix_LDL().GetTranspose() * m.GetDMatrix_LDL()
        * m.GetLTMatrix_LDL();

    if (m != mldl) {
      std::cout << "TestLdl: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << m << std::endl;
      std::cout << "L D L^T:" << std::endl << mldl << std::endl;
      std::cout << "L D L^T - m" << std::endl << mldl - m << std::endl;
      return;
    }
  }

  std::cout << "TestLdl: completed." << std::endl;
}

void TestSystemSolution_Symmetric(int number_of_tests, int matrix_size) {
  std::vector<std::vector<double>>
      v1(matrix_size, std::vector<double>(matrix_size));
  std::vector<std::vector<double>>
      v2(matrix_size, std::vector<double>(1));

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j <= i; ++j) {
        v1[i][j] = Random();
        v1[j][i] = v1[i][j];
      }
      v2[i][0] = Random();
    }

    Matrix<double> a(v1, 0.0001);
    Matrix<double> b(v2, 0.0001);
    auto sol = a.SolveSystem_Symmetric(b);

    if (a * sol != b) {
      std::cout << "TestSystemSolution_Symmetric: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << a << std::endl;
      std::cout << "Solution:" << std::endl << sol << std::endl;
      std::cout << "B:" << std::endl << b << std::endl;
      std::cout << "Corresponding B:" << std::endl << a * sol << std::endl;
      std::cout << "B - AX:" << std::endl << b - a * sol << std::endl;
      return;
    }
  }

  std::cout << "TestSystemSolution_Symmetric: completed." << std::endl;
}

void TestSystemSolution_Tridiagonal(int number_of_tests, int matrix_size) {
  std::vector<std::vector<double>> v1(matrix_size, std::vector<double>(4));
  std::vector<std::vector<double>> v2(matrix_size, std::vector<double>(1));

  for (int k = 0; k < number_of_tests; ++k) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < 3; ++j) {
        v1[i][j] = Random();
      }
      v1[i][3] = 0;
      v2[i][0] = Random();
    }

    Matrix<double> a(v1, 0.00001);
    Matrix<double> b(v2, 0.00001);
    auto sol = a.SolveSystem_Tridiagonal(b);

    if (a.GetTridiagonalMatrixAsNormal() * sol != b) {
      std::cout << "TestSystemSolution_Tridiagonal: FAIL!" << std::endl;
      std::cout << "Matrix:" << std::endl << a.GetTridiagonalMatrixAsNormal()
                << std::endl;
      std::cout << "Solution:" << std::endl << sol << std::endl;
      std::cout << "B:" << std::endl << b << std::endl;
      std::cout << "Corresponding B:" << std::endl
                << a.GetTridiagonalMatrixAsNormal() * sol << std::endl;
      std::cout << "B - AX:" << std::endl
                << b - a.GetTridiagonalMatrixAsNormal() * sol << std::endl;
      return;
    }
  }

  std::cout << "TestSystemSolution_Tridiagonal: completed." << std::endl;
}

void TestSystemSolution_Sor(int number_of_tests) {
  std::mt19937_64 random_generator
      (std::chrono::system_clock::now().time_since_epoch().count());

  for (int i = 0; i < number_of_tests; ++i) {
    int n = 1 + static_cast<int>(random_generator() % 3000);
    auto x = Task5<double>::SolveSystem(n);
    auto b = Task5<double>::GetAMatrix(n) * x.matrix;

    std::vector<std::vector<double>> id_vector(n, std::vector<double>(1, 1));
    Matrix id(id_vector);

    if (b != id) {
      std::cout << "TestSystemSolution_Sor: FAIL!" << std::endl;
      std::cout << "n:" << std::endl << n << std::endl;
      std::cout << "Solution:" << std::endl << x.matrix << std::endl;
      std::cout << "b:" << std::endl << b << std::endl;
      std::cout << "id - AX:" << std::endl << id - b << std::endl;
      return;
    }
  }

  std::cout << "TestSystemSolution_Sor: completed." << std::endl;
}

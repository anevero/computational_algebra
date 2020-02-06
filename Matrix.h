#ifndef MATRIX_H_
#define MATRIX_H_

// Base Matrix template class. Requires template arguments satisfying
// std::is_floating_point<T> condition.
// A constructor from the two-dimensional vector is provided, as well as
// comparison operators (equal and not equal), sum, difference and product
// operators. << operator is overloaded for printing matrices.
// Methods for counting TLU / LDL / QR decomposition, solving systems of linear
// equations, counting an inverse matrix and counting a condition number,
// searching for eigenvalues and eigenvectors are implemented.
// Their descriptions (often with time asymptotics estimate) are
// located below.

#include <algorithm>
#include <cmath>
#include <complex>
#include <condition_variable>
#include <iostream>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "Polynomial.h"
#include "Utils/ThreadPool.h"
#include "Utils/Utils.h"

namespace matrix {

using polynomial::Polynomial;
using matrix_utils::ThreadPool;
using matrix_utils::Equal;

template<class T>
concept MatrixNumber = std::is_floating_point_v<T>;

template<MatrixNumber T>
class Matrix {
 public:
  struct Eigenvector {
    Matrix<T> vector;
    std::complex<T> value;
  };

 public:
  explicit Matrix(std::vector<std::vector<T>> matrix,
                  T epsilon = std::numeric_limits<T>::epsilon());
  // Default constructor initializes the matrix with [0].
  Matrix();
  ~Matrix() = default;

// ---------------------------------------------------------------------------
// Comparison operators. Use predefined epsilons and Equal function from Utils
// to compare small floating point values properly.

  bool operator==(const Matrix& other) const;
  bool operator!=(const Matrix& other) const;

// ---------------------------------------------------------------------------
// Arithmetic operators.

  Matrix operator+(const Matrix& other) const;
  Matrix operator-(const Matrix& other) const;
  Matrix operator*(const Matrix& other) const;
  Matrix operator*(T number) const;
  Matrix operator/(T number) const;
  template<class U>
  friend Matrix<U> operator*(U number, const Matrix<U>& matrix);

// ---------------------------------------------------------------------------
// Printing the matrix.

  [[nodiscard]] std::string ToString() const;

  // Operator << is overloaded outside the class.

// ---------------------------------------------------------------------------
// Matrix norm functions. Matrix norm, induced by the vector max-norm, is used.

  void CountNorm();
  std::optional<T> GetNorm() const;
  T CountAndGetNorm();

// ---------------------------------------------------------------------------
// Matrix getters.

  int GetNumberOfRows() const;
  int GetNumberOfColumns() const;
  T GetEpsilon() const;
  [[nodiscard]] std::vector<std::vector<T>> GetMatrixAsVector() const;

// ---------------------------------------------------------------------------
// TLU decomposition, solving systems of linear equations, counting the inverse
// matrix and the condition number for ANY matrix.

  // Counts TLU decomposition with choosing a main element in the column (so
  // rows can be permutated). If the matrix is singular, throws an exception.
  // Time complexity: O(n^3).
  void CountTluDecomposition();
  // Returns N*1 matrix, which is the solution of AX = b system. b should be
  // N*1 matrix, where N is equal to the number of columns in A matrix.
  // 'this' is used as A matrix.
  // Uses TLU decomposition for solving the system. If it hasn't been
  // counted yet, runs CountTluDecomposition method.
  // Time complexity (without counting TLU): O(n^2).
  [[nodiscard]] Matrix<T> SolveSystem(Matrix<T> b);
  // Counts an inverse matrix using TLU decomposition and SolveSystem methods.
  // Takes O(n^3) time.
  void CountInverseMatrix();
  // Counts a condition number using matrix norm, induced by the vector
  // max-norm.
  // If an inverse matrix hasn't been counted yet, runs CountInverseMatrix
  // function. Time complexity (without counting the inverse matrix): O(n^2).
  void CountConditionNumber();

// ---------------------------------------------------------------------------
// TLU decomposition, solving systems of linear equations, counting the inverse
// matrix for the ALMOST TRIANGULAR matrix (look at the task 1 conditions).
// Almost triangulat matrix in this task is a lower Hessenberg matrix.
// There're two ways for counting the inverse matrix. The both have pros and
// cons, which are described in the report.
// CountInverseMatrix_AlmostTriangular_TLU counts the inverse matrix using
// the optimized algorithm for TLU decomposition.
// CountInverseMatrix_AlmostTriangular counts the inverse matrix making
// the matrix lower-triangular and using this property.

  // Checks if the matrix is satisfying Task 1 conditions (is almost triangular
  // or lower Hessenberg, i.e. has only one non-zero element to the right
  // of the diagonal).
  bool IsAlmostTriangular() const;

  // Counts TLU decomposition of an almost triangular matrix. Uses this fact
  // to reduce the number of subtraction operations.
  // No elements are chosen as main, because it will break the structure of the
  // matrix. If the current element is equal to zero, and it's not possible to
  // continue applying this algorithm, the function throws an exception.
  // Overall time complexity is O(n^2).
  void CountTluDecomposition_AlmostTriangular();
  // Works just as common SolveSystem function, but uses the fact, that
  // U matrix from the decomposition of the almost triangular matrix is almost
  // diagonal. So we can reduce the number of subtraction operations.
  [[nodiscard]] Matrix<T> SolveSystem_AlmostTriangular(Matrix<T> b);
  // Counts an inverse matrix using corresponding methods for the almost
  // triangular matrices. Time complexity: O(n^3).
  void CountInverseMatrix_AlmostTriangular_Tlu();
  // Counts inverse matrix using corresponding methods for the almost
  // triangular matrices and not using multithreading.
  // Time complexity: O(n^3).
  void CountInverseMatrix_AlmostTriangular_Tlu_SingleThread();

  // Using the fact, that the matrix is almost triangular, makes it
  // lower-triangular and counts the inverse. Overall time complexity is
  // O(n^3).
  void CountInverseMatrix_AlmostTriangular();
  // The same as previous, but uses only one thread.
  void CountInverseMatrix_AlmostTriangular_SingleThread();

// ---------------------------------------------------------------------------
// LDL decomposition, solving systems of linear equations for the SYMMETRIC
// matrix.

  // Checks if the matrix is symmetric.
  bool IsSymmetric() const;
  // Returns the transpose of the current matrix. Its generation takes O(n^2)
  // time.
  [[nodiscard]] Matrix<T> GetTranspose() const;

  // Counts LDL decomposition of the symmetric matrix (using square roots
  // algorithm).
  // No elements are chosen as main, because it will break the structure of the
  // matrix. If a current element is equal to zero, and it's not possible to
  // continue applying this algorithm, the function throws an exception.
  // Time complexity is O(n^3).
  void CountLdlDecomposition_Symmetric();
  // Returns N*1 matrix, which is the solution of AX = b system. b should be
  // N*1 matrix, where N is equal to the number of columns in A matrix.
  // 'this' is used as A matrix.
  // Uses LDL decomposition for solving the system. If it hasn't been
  // counted yet, runs CountLdlDecomposition method.
  // Time complexity (without counting LDL): O(n^2).
  [[nodiscard]] Matrix<T> SolveSystem_Symmetric(Matrix<T> b);

// ---------------------------------------------------------------------------
// Solving systems of linear equations for the TRIDIAGONAL matrix.
// Tridiagonal matrix representation must include three columns (the first
// column is the element to the left from the main diagonal, and so on).

  // Checks if the current matrix can represent the tridiagonal matrix (i.e
  // the current matrix has three columns).
  bool IsTridiagonal() const;
  // Returns the tridiagonal matrix represented as three columns of values
  // in the normal square matrix representation.
  [[nodiscard]] Matrix<T> GetTridiagonalMatrixAsNormal() const;

  // Returns N*1 matrix, which is the solution of AX = b system. b should be
  // N*1 matrix, where N is equal to the number of columns in A matrix.
  // 'this' is used as A matrix.
  // Uses tridiagonal matrices properties for solving the system.
  // Time complexity: O(n).
  [[nodiscard]] Matrix<T> SolveSystem_Tridiagonal(Matrix<T> b) const;

// ---------------------------------------------------------------------------
// QR decomposition for any matrix.

 private:
  // A function to multiply the matrix by the rotation matrix in O(n)
  // time (matrix = rotation_matrix * matrix). Some tricks are used to make
  // compiler apply vectorization here.
  // This function is used in QR decomposition algorithm and in QR algorithm
  // for getting all the eigenvalues of the matrix.
  void MultiplyMatrixByRotation_Left(std::vector<std::vector<T>>* matrix,
                                     T sin, T cos,
                                     int i, int j, int size) const;

  // The same as previous, but performs multiplication like:
  // matrix = matrix * rotation_matrix.
  void MultiplyMatrixByRotation_Right(std::vector<std::vector<T>>* matrix,
                                      T sin, T cos,
                                      int i, int j, int size) const;

 public:
  // Counts QR decomposition of the matrix (using Givens method). Internal
  // orthogonal matrices are multiplied efficiently (only two rows and two
  // columns are changed) to improve time complexity.
  // Time complexity is O(n^3).
  void CountQrDecomposition();

// ---------------------------------------------------------------------------
// QR algorithm for any matrix (based on upper Hessenberg matrices). If the
// matrix is suitable for the algorithm, and it is able to converge properly,
// it allows to get all the eigenvalues of the matrix. Otherwise eigenvalues,
// returned by the algorithms below, can be incorrrect.
// If the matrix is symmetric, the algorithm is also able to return
// eigenvectors.
// If the matrix is not symmetric, returned eigenvectors will be incorrect
// (they will be orthogonal vectors, which are formed from the eigenvectors
// by the Gram-Schmidt orthogonalization process).

  // Counts an upper Hessenberg (almost triangular) matrix, which is similar
  // to this. This matrix is saved in the internal field of the class.
  // Time complexity is O(n^3).
  void CountUpperHessenbergMatrix();

 private:
  // Runs one iteration of QR algorithm. One iteration updates current
  // Hessenberg matrix (with O(n^2) time).
  // Saves the information about rotation matrices in the corresponding private
  // field of the class. This info is necessary for restoring the eigenvectors
  // later.
  // This method assumes, that the Hessenberg form of the matrix has been
  // already counted.
  void RunQRAlgorithmIteration();

 public:
  // Runs QR algorithm iterations, until the Hessenberg matrix converges to the
  // matrix with eigenvalues on the diagonal.
  // Convergence condition: for every item on the diagonal and to the left of
  // diagonal the difference between it and its version at the previous
  // iteration should be <= epsilon.
  // If the algorithm converges slowly (or even diverges), it will be stopped
  // after some number of iterations. The results will be unpredictable, the
  // information about this situation will be printed to std::cerr.
  // If the Hessenberg form of the matrix hasn't been counted yet, runs
  // the corresponding method.
  void RunQrAlgorithm(T epsilon = 0.00001,
                      int max_number_of_iterations = 1000000);

  // Returns a vector with all the eigenvalues of the matrix. These eigenvalues
  // are extracted from the Hessenberg matrix after applying QR algorithm to
  // it. This method assumes, that the Hessenberg form of the matrix has been
  // already counted.
  // Returned vectors are correct eigenvectors for the symmetric matrices.
  // Otherwise they will be incorrect (they will be orthogonal vectors, which
  // are formed from the eigenvectors by the Gram-Schmidt orthogonalization
  // process).
  // NOTE: it's recommended to use the same epsilon, as with RunQrAlgorithm()
  // method, or even smaller. Otherwise the results can be incorrect.
  // Time complexity is O(n^2).
  std::vector<Eigenvector>
  GetEigenvectorsFromHessenbergMatrix(T epsilon = 0.00001) const;

// ---------------------------------------------------------------------------
// Danilevsky algorithm implementation for any matrix (based on upper
// Frobenius matrices). This algorithm is used to count the characteristic
// polynomial of the matrix. Then at least real roots can be extracted from it
// (look at Polynomial class methods), and these roots will be the eigenvalues.
// If the Frobenius matrix, similar to this, isn't a block matrix, the
// transition matrix S (F = S^{-1} A S) can be used to get the eigenvectors,
// corresponding to the eigenvalues.
// Note that due to restrictions of our methods to find the roots of the
// polynomials, behaviour of the algorithms below can be unpredictable, if
// the characteristic polynomial isn't suitable (for example, has a lot of
// roots with degree > 1).

  // Counts the Frobenius (possibly block) matrix, which is similar to this.
  // Time complexity is O(n^3). This matrix is stored in the internal field
  // of the class.
  void CountFrobeniusMatrix();

  // Counts characteristic polynomial for the matrix, using its Frobenius form.
  // Saves it in the internal characteristic polynomial field.
  // If the Frobenius matrix hasn't been counted yet, runs CountFrobeniusMatrix
  // method.
  // Time complexity: O(n).
  void CountCharacteristicPolynomial();

  // Calls FindRoots() method for the characteristic polynomial. This method
  // assumes that the characteristic polynomial has been already counted.
  void FindCharacteristicPolynomialRoots();

  // This method uses eigenvalues (the roots of the characteristic polynomial)
  // and Frobenius transition matrix to count the eigenvectors,
  // corresponding to eigenvalues.
  // The method assumes that the Frobenius form of the matrix and the roots
  // of the characteristic polynomial have been already counted.
  // If some eigenvalues have degree greater than 1, behaviour is undefined:
  // the method can return several similar eigenvectors or just one eigenvector
  // for it.
  // If the Frobenius matrix is a block matrix, returned eigenvectors will be
  // incorrect.
  std::vector<Eigenvector> GetEigenvectorsFromFrobeniusMatrix() const;

// ---------------------------------------------------------------------------
// Power iteration algorithm for any matrix. Allows to get the eigenvalue with
// maximum module and the corresponding eigenvector. Works only with real
// numbers.

  // Runs Power Iterations method to count the eigencalue and the eigenvector.
  // Returned vector can consist of 0, 1 or 2 eigenvectors and corresponding
  // values.
  // 0: the method hasn't converged (it's possible that eigenvalues are too
  // similar, or we just need even more iterations, or the eigenvalues are
  // not real).
  // 1: the matrix has an eigenvector with eigenvalue, which module is greater
  // or equal than others. This eigenvalue and eigenvector are returned.
  // 2: the matrix has two eigenvectors with opposite eigenvalues, which
  // modules are greater or equal than others, OR one of the returned vectors
  // will be the vector with maximum module, the second will be some another
  // vector.
  // Note that complex eigenvalues and corresponding eigenvectors aren't taken
  // into account.
  std::vector<Eigenvector> GetPowerIterationResults(
      T epsilon = 0.00001, int max_number_of_iterations = 1500) const;

// ---------------------------------------------------------------------------
// Getters for the results of TLU decomposition, LDL decomposition, QR
// decomposition, counting the inverse matrix and the condition number,
// counting the Hessenberg and Frobenius forms of the matrix, etc.
// All the methods below assume that necessary matrices have been already
// counted. Otherwise the returned value is undefined.

  [[nodiscard]] Matrix<T> GetLMatrix_TLU() const;
  [[nodiscard]] Matrix<T> GetUMatrix_TLU() const;
  // Returns T matrix in normal form (not in the internal form of a
  // one-dimensional vector).
  [[nodiscard]] Matrix<T> GetTMatrix_TLU() const;
  // Returns T^{-1} matrix in normal form (not in the internal form of a
  // one-dimensional vector).
  [[nodiscard]] Matrix<T> GetTInverseMatrix_TLU() const;

  [[nodiscard]] Matrix<T> GetLTMatrix_LDL() const;
  [[nodiscard]] Matrix<T> GetDMatrix_LDL() const;

  // GetQMatrix returns Q matrix in normal (not transpose) state.
  [[nodiscard]] Matrix<T> GetQMatrix_QR() const;
  [[nodiscard]] Matrix<T> GetRMatrix_QR() const;

  // Returns an upper Hessenberg matrix, similar to this.
  [[nodiscard]] Matrix<T> GetUpperHessenbergMatrix() const;

  // Returns the Frobenius matrix (possibly, a block matrix), similar to this.
  [[nodiscard]] Matrix<T> GetFrobeniusMatrix() const;

  // Returns the characteristic polynomial of the matrix.
  [[nodiscard]] Polynomial<T> GetCharacteristicPolynomial() const;

  // Runs characteristic_polynomial.GetRoots() method (so this method can
  // actually be replaced by GetCharacteristicPolynomial() and manual
  // searching for the roots).
  [[nodiscard]] std::vector<T> GetCharacteristicPolynomialRealRoots() const;

  [[nodiscard]] Matrix<T> GetInverseMatrix() const;
  std::optional<T> GetConditionNumber() const;

 private:
  std::vector<std::vector<T>> matrix_;
  std::optional<T> norm_ = std::nullopt;

  int rows_;
  int columns_;
  T epsilon_;

// ---------------------------------------------------------------------------
// Variables connected with TLU decomposition.

  std::vector<std::vector<T>> L_matrix_TLU_{};
  std::vector<std::vector<T>> U_matrix_TLU_{};
  // Can be used while solving systems: AX = B will be equal to LU = TB.
  // To save the memory, this matrix is stored as a single-dimensional vector.
  // Number on the ith position points to the coordinate of 1 in the ith row.
  std::vector<int> T_matrix_TLU_{};
  // Can be used to express A matrix as the product: A = T^{-1} LU.
  // To save the memory, this matrix is stored as a single-dimensional vector.
  // Number on the ith position points to the coordinate of 1 in the ith column.
  std::vector<int> T_inverse_matrix_TLU_{};

// ---------------------------------------------------------------------------
// Variables connected with LDL decomposition.

  std::vector<std::vector<T>> LT_matrix_LDL_{};
  // To save the memory, this matrix is stored as a single-dimensional vector.
  // Number on the ith position points to the coordinate of 1 in the ith row.
  std::vector<int> D_matrix_LDL_{};

// ---------------------------------------------------------------------------
// Variables connected with QR decomposition.

  // Q matrix is stored in transpose state (so to get Q matrix for A = QR
  // equality you need to transpose Q_matrix_QR_).
  std::vector<std::vector<T>> Q_matrix_QR_{};
  std::vector<std::vector<T>> R_matrix_QR_{};

// ---------------------------------------------------------------------------
// Variables connected with QR algorithm.

  // Upper Hessenberg matrix, similar to 'this', is stored here. This variable
  // is updated every iteration of the QR algorithm, and finally it will
  // have all the eigenvalues at the cells of the main diagonal (it the method
  // converges for the given matrix).
  std::vector<std::vector<T>> hessenberg_matrix_{};

  // Stored Q matrix, which is got after applying QR algorithm (this Q
  // matrix will consist of eigenvectors in case of the symmetric matrix).
  std::vector<std::vector<T>> hessenberg_rotation_matrix_{};

// ---------------------------------------------------------------------------
// Variables connected with Danilevsky algorithm.

  // Frobenius (it's possible, a block) matrix, which is similar to this
  // matrix. Counting this matrix is a part of the Danilevsky algorithm for
  // counting the characteristic polynomial of the matrix.
  std::vector<std::vector<T>> frobenius_matrix_;

  // If the Frobenius matrix is found as F = S^{-1} A S, then this matrix is S.
  // It's necessary for restoring eigenvectors after applying Danilevsky
  // algorithm and counting eigenvalues.
  std::vector<std::vector<T>> frobenius_transition_matrix_;

  Polynomial<T> characteristic_polynomial_;

// ---------------------------------------------------------------------------
// Results of applying different algorithms.

  std::vector<std::vector<T>> inverse_matrix_{};
  std::optional<T> condition_number_ = std::nullopt;
};

template <MatrixNumber T>
Matrix<T>::Matrix(std::vector<std::vector<T>> matrix, T epsilon)
    : matrix_(std::move(matrix)),
      epsilon_(epsilon) {
  if (matrix_.empty() || matrix_[0].empty()) {
    throw std::invalid_argument("Empty matrix can't be constructed.");
  }
  if (epsilon_ <= 0) {
    throw std::invalid_argument("Epsilon must be positive.");
  }
  rows_ = matrix_.size();
  columns_ = matrix_[0].size();
}

template <MatrixNumber T>
Matrix<T>::Matrix() : Matrix({{0}}) {
}

// ---------------------------------------------------------------------------
// Comparison operators.

template <MatrixNumber T>
bool Matrix<T>::operator==(const Matrix& other) const {
  if (rows_ != other.rows_ || columns_ != other.columns_) return false;
  auto max_epsilon = std::max(epsilon_, other.epsilon_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < columns_; ++j) {
      if (std::abs(matrix_[i][j] - other.matrix_[i][j]) >= max_epsilon) {
        return false;
      }
    }
  }
  return true;
}

template <MatrixNumber T>
bool Matrix<T>::operator!=(const Matrix& other) const {
  return !operator==(other);
}

// ---------------------------------------------------------------------------
// Arithmetic operators.

template <MatrixNumber T>
Matrix<T> Matrix<T>::operator+(const Matrix& other) const {
  if (rows_ != other.rows_ || columns_ != other.columns_) {
    throw std::runtime_error(
        "Matrices have different sizes; the sum can't be counted.");
  }

  Matrix<T> result(matrix_, std::max(epsilon_, other.epsilon_));
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < columns_; ++j) {
      result.matrix_[i][j] += other.matrix_[i][j];
    }
  }
  return result;
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::operator-(const Matrix& other) const {
  if (rows_ != other.rows_ || columns_ != other.columns_) {
    throw std::runtime_error(
        "Matrices have different sizes; the difference can't be counted.");
  }

  Matrix<T> result(matrix_, std::max(epsilon_, other.epsilon_));
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < columns_; ++j) {
      result.matrix_[i][j] -= other.matrix_[i][j];
    }
  }
  return result;
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::operator*(const Matrix& other) const {
  if (columns_ != other.rows_) {
    throw std::runtime_error(
        "Matrices are not consistent; the product can't be counted.");
  }

  std::vector<std::vector<T>> result(rows_, std::vector<T>(other.columns_, 0));

  int number_of_threads = static_cast<int>(std::thread::hardware_concurrency());
  number_of_threads = (number_of_threads == 0) ? 8 : number_of_threads;
  std::vector<std::thread> threads(number_of_threads);
  int rows_per_thread;
  int current_row = 0;
  int rows_remaining = rows_;

  for (int i = 0; i < number_of_threads; ++i) {
    rows_per_thread = rows_remaining / (number_of_threads - i);
    threads[i] = std::thread(
        [this, &result, &other, current_row, rows_per_thread]() {
          for (int r = current_row; r < current_row + rows_per_thread; ++r) {
            for (int j = 0; j < other.columns_; ++j) {
              for (int k = 0; k < columns_; ++k) {
                result[r][j] += matrix_[r][k] * other.matrix_[k][j];
              }
            }
          }
        });

    current_row += rows_per_thread;
    rows_remaining -= rows_per_thread;
  }

  for (int i = 0; i < number_of_threads; ++i) {
    threads[i].join();
  }

  return Matrix<T>(result, std::max(epsilon_, other.epsilon_));
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::operator*(T number) const {
  Matrix<T> result(matrix_, epsilon_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < columns_; ++j) {
      result.matrix_[i][j] *= number;
    }
  }
  return result;
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::operator/(T number) const {
  Matrix<T> result(matrix_, epsilon_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < columns_; ++j) {
      result.matrix_[i][j] /= number;
    }
  }
  return result;
}

template<class U>
Matrix<U> operator*(U number, const Matrix<U>& matrix) {
  return matrix.operator*(number);
}

// ---------------------------------------------------------------------------
// Printing the matrix.

template <MatrixNumber T>
std::string Matrix<T>::ToString() const {
  std::stringstream stream;
  stream << '[';
  for (int i = 0; i < rows_; ++i) {
    if (i != 0) {
      stream << ' ';
    }
    stream << '[';
    for (int j = 0; j < columns_; ++j) {
      stream << matrix_[i][j];
      if (j != columns_ - 1) {
        stream << ", ";
      }
    }
    stream << ']';
    if (i != rows_ - 1) {
      stream << ',' << std::endl;
    }
  }
  stream << ']';
  return stream.str();
}

template<class U>
std::ostream& operator<<(std::ostream& out, const Matrix<U>& matrix) {
  return out << matrix.ToString();
}

// ---------------------------------------------------------------------------
// Matrix norm functions.

template <MatrixNumber T>
void Matrix<T>::CountNorm() {
  norm_ = 0;
  for (int i = 0; i < rows_; ++i) {
    T current_sum = 0;
    for (int j = 0; j < columns_; ++j) {
      current_sum += std::abs(matrix_[i][j]);
    }
    norm_ = std::max(norm_.value(), current_sum);
  }
}

template <MatrixNumber T>
std::optional<T> Matrix<T>::GetNorm() const {
  return norm_;
}

template <MatrixNumber T>
T Matrix<T>::CountAndGetNorm() {
  if (!norm_.has_value()) {
    CountNorm();
  }
  return norm_.value();
}

// ---------------------------------------------------------------------------
// Getters.

template <MatrixNumber T>
int Matrix<T>::GetNumberOfRows() const {
  return rows_;
}

template <MatrixNumber T>
int Matrix<T>::GetNumberOfColumns() const {
  return columns_;
}

template <MatrixNumber T>
T Matrix<T>::GetEpsilon() const {
  return epsilon_;
}

template <MatrixNumber T>
std::vector<std::vector<T>> Matrix<T>::GetMatrixAsVector() const {
  return matrix_;
}

// ---------------------------------------------------------------------------
// TLU decomposition, solving systems of linear equations, counting the inverse
// matrix and the condition number for ANY matrix.

template <MatrixNumber T>
void Matrix<T>::CountTluDecomposition() {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }
  int size = rows_;

  // Initializing L with zero matrix, U with main matrix, T and T^{-1} with
  // identity matrices.
  L_matrix_TLU_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  U_matrix_TLU_ = matrix_;

  T_matrix_TLU_ = std::vector<int>(size, 0);
  T_inverse_matrix_TLU_ = std::vector<int>(size, 0);
  std::iota(T_matrix_TLU_.begin(), T_matrix_TLU_.end(), 0);
  std::iota(T_inverse_matrix_TLU_.begin(), T_inverse_matrix_TLU_.end(), 0);

  // Filling L, U and T matrices out. Takes O(n^3) time.
  for (int i = 0; i < size; ++i) {
    // Searching for the biggest element in the column. Takes O(n) time.
    int index_of_biggest = i;
    for (int j = i; j < size; ++j) {
      if (std::abs(U_matrix_TLU_[j][i])
          > std::abs(U_matrix_TLU_[index_of_biggest][i])) {
        index_of_biggest = j;
      }
    }

    // Swapping rows in L and U matrices. Swapping vectors works in O(1)
    // time (std::vector::swap just swaps internal pointers), so updating
    // these matrices will take O(1) time.
    U_matrix_TLU_[i].swap(U_matrix_TLU_[index_of_biggest]);
    L_matrix_TLU_[i].swap(L_matrix_TLU_[index_of_biggest]);

    // Swapping elements in T and T^{-1} matrices. Takes O(1) time.
    std::swap(T_matrix_TLU_[i], T_matrix_TLU_[index_of_biggest]);
    std::swap(T_inverse_matrix_TLU_[i],
              T_inverse_matrix_TLU_[index_of_biggest]);

    if (Equal(U_matrix_TLU_[i][i], T(), epsilon_)) {
      L_matrix_TLU_.clear();
      U_matrix_TLU_.clear();
      T_matrix_TLU_.clear();
      T_inverse_matrix_TLU_.clear();
      throw std::runtime_error("Matrix is singular.");
    }

    // Updating current column in L matrix. Takes O(n) time.
    auto multiplier = U_matrix_TLU_[i][i];
    for (int j = i; j < size; ++j) {
      L_matrix_TLU_[j][i] = U_matrix_TLU_[j][i] / multiplier;
    }

    // Updating current column (filling it with zeros) in U matrix. Takes
    // O(n^2) time.
    for (int j = i + 1; j < size; ++j) {
      multiplier = (-1) * U_matrix_TLU_[j][i] / U_matrix_TLU_[i][i];
      for (int k = i; k < size; ++k) {
        U_matrix_TLU_[j][k] += multiplier * U_matrix_TLU_[i][k];
      }
    }
  }
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::SolveSystem(Matrix<T> b) {
  if (columns_ != b.GetNumberOfRows() || b.GetNumberOfColumns() != 1) {
    throw std::invalid_argument("B is not a proper vector for this matrix.");
  }

  if (L_matrix_TLU_.empty()) {
    CountTluDecomposition();
  }
  int size = rows_;

  // Counting the product of T matrix and b matrix efficiently (in linear time).
  std::vector<std::vector<T>> new_b(size, std::vector<T>{});
  for (int i = 0; i < size; ++i) {
    new_b[i] = b.matrix_[T_matrix_TLU_[i]];
  }
  b.matrix_ = std::move(new_b);

  // Solving L(UX) = b system using the fact that L is lower-triangular.
  // Takes O(n^2) time.
  for (int i = 0; i < size; ++i) {
    b.matrix_[i][0] /= L_matrix_TLU_[i][i];
    auto multiplier = (-1) * b.matrix_[i][0];
    for (int j = i + 1; j < size; ++j) {
      b.matrix_[j][0] += multiplier * L_matrix_TLU_[j][i];
    }
  }

  // Solving UX = b system using the fact that U is upper-triangular.
  // Takes O(n^2) time.
  for (int i = size - 1; i >= 0; --i) {
    b.matrix_[i][0] /= U_matrix_TLU_[i][i];
    auto multiplier = (-1) * b.matrix_[i][0];
    for (int j = i - 1; j >= 0; --j) {
      b.matrix_[j][0] += multiplier * U_matrix_TLU_[j][i];
    }
  }

  return b;
}

template <MatrixNumber T>
void Matrix<T>::CountInverseMatrix() {
  if (L_matrix_TLU_.empty()) {
    CountTluDecomposition();
  }
  int size = rows_;

  inverse_matrix_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  std::vector<Matrix<T>> identity_columns(
      size, Matrix<T>(std::vector<std::vector<T>>(size, std::vector<T>(1, 0))));
  for (int i = 0; i < size; ++i) {
    identity_columns[i].matrix_[i][0] = 1;
  }

  // Solving Ax = e systems, where A is 'this' matrix, x is the column of the
  // inverse matrix, e is the column of the identity matrix.
  // These systems can be solved independently for all the columns of the
  // inverse matrix, so we will use multithreading here.
  // Threads are created only once, so no need for the thread pools here.

  int number_of_threads = static_cast<int>(std::thread::hardware_concurrency());
  number_of_threads = (number_of_threads == 0) ? 8 : number_of_threads;
  std::vector<std::thread> threads(number_of_threads);
  int columns_per_thread;
  int current_column = 0;
  int columns_remaining = size;

  for (int i = 0; i < number_of_threads; ++i) {
    columns_per_thread = columns_remaining / (number_of_threads - i);
    threads[i] = std::thread(
        [this, &identity_columns, current_column, columns_per_thread, size]() {
          for (int k = current_column;
               k < current_column + columns_per_thread;
               ++k) {
            auto inverse_column = SolveSystem(identity_columns[k]).matrix_;
            for (int j = 0; j < size; ++j) {
              inverse_matrix_[j][k] = inverse_column[j][0];
            }
          }
        });

    current_column += columns_per_thread;
    columns_remaining -= columns_per_thread;
  }

  for (int i = 0; i < number_of_threads; ++i) {
    threads[i].join();
  }
}

template <MatrixNumber T>
void Matrix<T>::CountConditionNumber() {
  if (inverse_matrix_.empty()) {
    CountInverseMatrix();
  }
  if (!norm_.has_value()) {
    CountNorm();
  }

  int size = rows_;
  T inverse_norm = 0;

  for (int i = 0; i < size; ++i) {
    T current_sum = 0;
    for (int j = 0; j < size; ++j) {
      current_sum += std::abs(inverse_matrix_[i][j]);
    }
    inverse_norm = std::max(inverse_norm, current_sum);
  }

  condition_number_ = norm_.value() * inverse_norm;
}

// ---------------------------------------------------------------------------
// TLU decomposition, solving systems of linear equations, counting the inverse
// matrix for the ALMOST TRIANGULAR matrix (look at the task 1 conditions).

template <MatrixNumber T>
bool Matrix<T>::IsAlmostTriangular() const {
  if (rows_ != columns_) return false;
  for (int i = 0; i < rows_; ++i) {
    for (int j = i + 2; j < columns_; ++j) {
      if (matrix_[i][j] != 0) return false;
    }
  }
  return true;
}

template <MatrixNumber T>
void Matrix<T>::CountTluDecomposition_AlmostTriangular() {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }
  int size = rows_;

  // Initializing L with zero matrix, U with main matrix. T and T^{-1} matrices
  // are initialized, but not used after, because no permutations are allowed.
  L_matrix_TLU_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  U_matrix_TLU_ = matrix_;

  T_matrix_TLU_ = std::vector<int>(size, 0);
  T_inverse_matrix_TLU_ = std::vector<int>(size, 0);
  std::iota(T_matrix_TLU_.begin(), T_matrix_TLU_.end(), 0);
  std::iota(T_inverse_matrix_TLU_.begin(), T_inverse_matrix_TLU_.end(), 0);

  // Filling L and U matrices out. Takes O(n^2) time.
  for (int i = 0; i < size; ++i) {
    if (Equal(U_matrix_TLU_[i][i], T(), epsilon_)) {
      L_matrix_TLU_.clear();
      U_matrix_TLU_.clear();
      T_matrix_TLU_.clear();
      T_inverse_matrix_TLU_.clear();
      throw std::runtime_error(
          "This algorithm isn't suitable for this matrix.");
    }

    // Updating the current column in L matrix. Takes O(n) time.
    auto multiplier = U_matrix_TLU_[i][i];
    for (int j = i; j < size; ++j) {
      L_matrix_TLU_[j][i] = U_matrix_TLU_[j][i] / multiplier;
    }

    // Updating current column (filling it with zeros) in U matrix. Takes
    // O(n) time, because we need to perform subtraction only on two columns
    // (all the next elements of the current row are zeros).
    for (int j = i + 1; j < size; ++j) {
      multiplier = (-1) * U_matrix_TLU_[j][i] / U_matrix_TLU_[i][i];
      for (int k = i; k < i + 2 && k < size; ++k) {
        U_matrix_TLU_[j][k] += multiplier * U_matrix_TLU_[i][k];
      }
    }
  }
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::SolveSystem_AlmostTriangular(Matrix<T> b) {
  if (columns_ != b.GetNumberOfRows() || b.GetNumberOfColumns() != 1) {
    throw std::invalid_argument("B is not a proper vector for this matrix.");
  }

  if (L_matrix_TLU_.empty()) {
    CountTluDecomposition_AlmostTriangular();
  }
  int size = rows_;

  // Solving L(UX) = b system using the fact that L is lower-triangular.
  // Takes O(n^2) time.
  for (int i = 0; i < size; ++i) {
    b.matrix_[i][0] /= L_matrix_TLU_[i][i];
    auto multiplier = (-1) * b.matrix_[i][0];
    for (int j = i + 1; j < size; ++j) {
      b.matrix_[j][0] += multiplier * L_matrix_TLU_[j][i];
    }
  }

  // Solving UX = b system using the fact that U is upper-triangular. More,
  // U is almost diagonal now, so we can reduce the number of subtractions and
  // get O(n) time complexity.
  for (int i = size - 1; i >= 0; --i) {
    b.matrix_[i][0] /= U_matrix_TLU_[i][i];
    auto multiplier = (-1) * b.matrix_[i][0];
    for (int j = i - 1; j > i - 3 && j >= 0; --j) {
      b.matrix_[j][0] += multiplier * U_matrix_TLU_[j][i];
    }
  }

  return b;
}

template <MatrixNumber T>
void Matrix<T>::CountInverseMatrix_AlmostTriangular_Tlu() {
  if (L_matrix_TLU_.empty()) {
    CountTluDecomposition_AlmostTriangular();
  }
  int size = rows_;

  inverse_matrix_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  std::vector<Matrix<T>> identity_columns(
      size, Matrix<T>(std::vector<std::vector<T>>(size, std::vector<T>(1, 0))));
  for (int i = 0; i < size; ++i) {
    identity_columns[i].matrix_[i][0] = 1;
  }

  // Solving Ax = e systems, where A is 'this' matrix, x is the column of the
  // inverse matrix, e is the column of the identity matrix.
  // These systems can be solved independently for all the columns of the
  // inverse matrix, so we will use multithreading here.
  // Threads are created only once, so no need for the thread pools here.

  int number_of_threads = static_cast<int>(std::thread::hardware_concurrency());
  number_of_threads = (number_of_threads == 0) ? 8 : number_of_threads;
  std::vector<std::thread> threads(number_of_threads);
  int columns_per_thread;
  int current_column = 0;
  int columns_remaining = size;

  for (int i = 0; i < number_of_threads; ++i) {
    columns_per_thread = columns_remaining / (number_of_threads - i);
    threads[i] = std::thread(
        [this, &identity_columns, current_column, columns_per_thread, size]() {
          for (int k = current_column;
               k < current_column + columns_per_thread;
               ++k) {
            auto inverse_column =
                SolveSystem_AlmostTriangular(identity_columns[k]).matrix_;
            for (int j = 0; j < size; ++j) {
              inverse_matrix_[j][k] = inverse_column[j][0];
            }
          }
        });

    current_column += columns_per_thread;
    columns_remaining -= columns_per_thread;
  }

  for (int i = 0; i < number_of_threads; ++i) {
    threads[i].join();
  }
}

template <MatrixNumber T>
void Matrix<T>::CountInverseMatrix_AlmostTriangular_Tlu_SingleThread() {
  if (L_matrix_TLU_.empty()) {
    CountTluDecomposition_AlmostTriangular();
  }
  int size = rows_;
  inverse_matrix_ = std::vector<std::vector<T>>(size, std::vector<T>(size));
  Matrix<T> identity_column(
      std::vector<std::vector<T>>(size, std::vector<T>(1, 0)));

  for (int i = 0; i < size; ++i) {
    identity_column.matrix_[i][0] = 1;
    auto inverse_column = SolveSystem_AlmostTriangular(identity_column).matrix_;
    for (int j = 0; j < size; ++j) {
      inverse_matrix_[j][i] = inverse_column[j][0];
    }
    identity_column.matrix_[i][0] = 0;
  }
}

template <MatrixNumber T>
void Matrix<T>::CountInverseMatrix_AlmostTriangular() {
  int size = rows_;
  auto a_matrix = matrix_;
  inverse_matrix_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    inverse_matrix_[i][i] = 1;
  }

  // Making the matrix lower-triangular.
  for (int i = size - 1; i > 0; --i) {
    if (Equal(a_matrix[i][i], T(), epsilon_)) {
      inverse_matrix_.clear();
      throw std::runtime_error("Optimized algorithm cannot be applied to "
                               "this matrix.");
    }

    auto multiplier = (-1) * a_matrix[i - 1][i] / a_matrix[i][i];
    for (int j = 0; j <= i; ++j) {
      a_matrix[i - 1][j] += multiplier * a_matrix[i][j];
    }
    for (int j = i; j < size; ++j) {
      inverse_matrix_[i - 1][j] += multiplier * inverse_matrix_[i][j];
    }
  }

  if (Equal(a_matrix[0][0], T(), epsilon_)) {
    inverse_matrix_.clear();
    throw std::runtime_error("Optimized algorithm cannot be applied to "
                             "this matrix.");
  }

  // Creating the pool of threads. We will use different threads to parallelize
  // the subtraction operations. We will have a lot of 'groups' of operations
  // and will need to use the pool of threads, because creating new threads
  // every time is inefficient.
  int number_of_threads = static_cast<int>(std::thread::hardware_concurrency());
  number_of_threads = (number_of_threads == 0) ? 8 : number_of_threads;
  auto thread_pool = ThreadPool(number_of_threads);
  thread_pool.StartWorkers();

  int rows_completed = rows_;
  std::mutex mutex;
  std::condition_variable rows_completed_cv;

  // Subtracting the rows from the every row below using multithreading
  // (different rows of the matrix are processed by different threads).
  for (int i = 0; i < size; ++i) {
    // Lock the thread until all the subtraction operations from the previous
    // iterations won't be completed.
    {
      std::unique_lock<std::mutex> locker(mutex);
      rows_completed_cv.wait(locker, [this, &rows_completed]() {
        return (rows_completed == rows_);
      });
    }

    rows_completed = i + 1;
    int current_row = i + 1;
    int rows_remaining = size - i - 1;
    int rows_per_thread;

    for (int j = 0; j < number_of_threads; ++j) {
      rows_per_thread = rows_remaining / (number_of_threads - j);
      if (rows_per_thread == 0) continue;
      thread_pool.Schedule(
          [this, &a_matrix, &rows_completed, &mutex, &rows_completed_cv,
              current_row, rows_per_thread, size, i]() {
            for (int k = current_row; k < current_row + rows_per_thread; ++k) {
              if (Equal(a_matrix[i][i], T(), epsilon_)) {
                throw std::runtime_error(
                    "Optimized algorithm cannot be applied to "
                    "this matrix.");
              }
              auto multiplier = a_matrix[k][i] / a_matrix[i][i];
              for (int l = 0; l < size; ++l) {
                inverse_matrix_[k][l] -= multiplier * inverse_matrix_[i][l];
              }
            }
            {
              std::lock_guard<std::mutex> locker(mutex);
              rows_completed += rows_per_thread;
              rows_completed_cv.notify_one();
            }
          });
      current_row += rows_per_thread;
      rows_remaining -= rows_per_thread;
    }
  }

  // Dividing the rows by the diagonal elements.
  int rows_per_thread;
  int current_row = 0;
  int rows_remaining = size;

  for (int i = 0; i < number_of_threads; ++i) {
    rows_per_thread = rows_remaining / (number_of_threads - i);
    if (rows_per_thread == 0) continue;
    thread_pool.Schedule(
        [this, &a_matrix, current_row, rows_per_thread, size]() {
          for (int k = current_row; k < current_row + rows_per_thread; ++k) {
            auto multiplier = a_matrix[k][k];
            for (int j = 0; j < size; ++j) {
              inverse_matrix_[k][j] /= multiplier;
            }
          }
        });
    current_row += rows_per_thread;
    rows_remaining -= rows_per_thread;
  }
}

template <MatrixNumber T>
void Matrix<T>::CountInverseMatrix_AlmostTriangular_SingleThread() {
  int size = rows_;
  auto a_matrix = matrix_;
  inverse_matrix_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    inverse_matrix_[i][i] = 1;
  }

  // Making the matrix lower-triangular.
  for (int i = size - 1; i > 0; --i) {
    if (Equal(a_matrix[i][i], T(), epsilon_)) {
      inverse_matrix_.clear();
      throw std::runtime_error("Optimized algorithm cannot be applied to "
                               "this matrix.");
    }

    auto multiplier = (-1) * a_matrix[i - 1][i] / a_matrix[i][i];
    for (int j = 0; j <= i; ++j) {
      a_matrix[i - 1][j] += multiplier * a_matrix[i][j];
    }
    for (int j = i; j < size; ++j) {
      inverse_matrix_[i - 1][j] += multiplier * inverse_matrix_[i][j];
    }
  }

  if (Equal(a_matrix[0][0], T(), epsilon_)) {
    inverse_matrix_.clear();
    throw std::runtime_error("Optimized algorithm cannot be applied to "
                             "this matrix.");
  }

  for (int i = 0; i < size; ++i) {
    // Dividing the current row by the diagonal element.
    auto multiplier = a_matrix[i][i];
    for (int j = 0; j < size; ++j) {
      inverse_matrix_[i][j] /= multiplier;
    }

    // Subtracting this row from the every row below.
    for (int j = i + 1; j < size; ++j) {
      multiplier = a_matrix[j][i];
      for (int k = 0; k < size; ++k) {
        inverse_matrix_[j][k] -= multiplier * inverse_matrix_[i][k];
      }
    }
  }
}

// ---------------------------------------------------------------------------
// LDL decomposition, solving systems of linear equations for the SYMMETRIC
// matrix.

template <MatrixNumber T>
bool Matrix<T>::IsSymmetric() const {
  if (rows_ != columns_) return false;
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j <= i; ++j) {
      if (matrix_[i][j] != matrix_[j][i]) return false;
    }
  }
  return true;
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetTranspose() const {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }

  std::vector<std::vector<T>> transpose_matrix = matrix_;
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j <= i; ++j) {
      std::swap(transpose_matrix[i][j], transpose_matrix[j][i]);
    }
  }

  return Matrix(transpose_matrix, epsilon_);
}

template <MatrixNumber T>
void Matrix<T>::CountLdlDecomposition_Symmetric() {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }
  int size = rows_;

  // Initializing LT with 'this' matrix, D with 'identity' matrix.
  LT_matrix_LDL_ = matrix_;
  D_matrix_LDL_ = std::vector<int>(size, 1);

  // Filling LT and D matrices out. Takes O(n^3) time.
  for (int i = 0; i < size; ++i) {
    // Checking if the current diagonal element is equal to zero.
    if (Equal(LT_matrix_LDL_[i][i], T(), epsilon_)) {
      LT_matrix_LDL_.clear();
      D_matrix_LDL_.clear();
      throw std::runtime_error(
          "LDL decomposition algorithm can't be applied to this matrix.");
    }

    // Subtracting the current row from all the rows below (and filling the
    // current column with zeros).
    // We use the symmetric property of the matrix here. Every lower right
    // submatrix we get at every step is symmetric, and we don't need to count
    // anything for the lower left part of it. So we can reduce the number of
    // subtractions.
    // Look at LDL algorithm description for the more detailed information.
    for (int j = i + 1; j < size; ++j) {
      auto multiplier = (-1) * LT_matrix_LDL_[i][j] / LT_matrix_LDL_[i][i];
      LT_matrix_LDL_[j][i] = 0;
      for (int k = j; k < size; ++k) {
        LT_matrix_LDL_[j][k] += multiplier * LT_matrix_LDL_[i][k];
      }
    }

    // Dividing the current row by the square root of the diagonal element.
    T multiplier;

    // Checking if the current diagonal element is negative.
    if (LT_matrix_LDL_[i][i] < 0) {
      multiplier = -std::sqrt(-LT_matrix_LDL_[i][i]);
      D_matrix_LDL_[i] = -1;
    } else {
      multiplier = std::sqrt(LT_matrix_LDL_[i][i]);
    }

    // Dividing the row.
    for (int j = i; j < size; ++j) {
      LT_matrix_LDL_[i][j] /= multiplier;
    }
  }
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::SolveSystem_Symmetric(Matrix<T> b) {
  if (columns_ != b.GetNumberOfRows() || b.GetNumberOfColumns() != 1) {
    throw std::invalid_argument("B is not a proper vector for this matrix.");
  }

  if (LT_matrix_LDL_.empty()) {
    CountLdlDecomposition_Symmetric();
  }
  int size = rows_;

  // Solving L(D L^T X) = b system, using the fact, that L is lower-triangular.
  // Takes O(n^2) time.
  for (int i = 0; i < size; ++i) {
    b.matrix_[i][0] /= LT_matrix_LDL_[i][i];
    auto multiplier = (-1) * b.matrix_[i][0];
    for (int j = i + 1; j < size; ++j) {
      b.matrix_[j][0] += multiplier * LT_matrix_LDL_[i][j];
    }
  }

  // Working with D matrix.
  for (int i = 0; i < size; ++i) {
    b.matrix_[i][0] *= D_matrix_LDL_[i];
  }

  // Solving L^T X = b system using the fact that L^T is upper-triangular.
  // Takes O(n^2) time.
  for (int i = size - 1; i >= 0; --i) {
    b.matrix_[i][0] /= LT_matrix_LDL_[i][i];
    auto multiplier = (-1) * b.matrix_[i][0];
    for (int j = i - 1; j >= 0; --j) {
      b.matrix_[j][0] += multiplier * LT_matrix_LDL_[j][i];
    }
  }

  return b;
}

// ---------------------------------------------------------------------------
// Solving systems of linear equations for the TRIDIAGONAL matrix.

template <MatrixNumber T>
bool Matrix<T>::IsTridiagonal() const {
  return (columns_ == 3);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetTridiagonalMatrixAsNormal() const {
  std::vector<std::vector<T>> result(rows_, std::vector<T>(rows_, 0));
  for (int i = 0; i < rows_; ++i) {
    if (i != 0) {
      result[i][i - 1] = matrix_[i][0];
    }
    result[i][i] = matrix_[i][1];
    if (i != rows_ - 1) {
      result[i][i + 1] = matrix_[i][2];
    }
  }
  return Matrix(result);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::SolveSystem_Tridiagonal(Matrix<T> b) const {
  if (rows_ != b.GetNumberOfRows() || b.GetNumberOfColumns() != 1) {
    throw std::invalid_argument("B is not a proper vector for this matrix.");
  }

  // Creating current matrix copy for storing temporary values.
  auto matrix = matrix_;

  // Adding the fourth column.
  for (int i = 0; i < rows_; ++i) {
    matrix[i].push_back(0);
  }

  // Going right and down.
  for (int i = 0; i < rows_ - 1; ++i) {
    // Swapping rows. Assuming that at the moment the fourth element of the
    // (i+1)th row is equal to zero.
    if (std::abs(matrix[i][1]) < std::abs(matrix[i + 1][0])) {
      std::swap(matrix[i][1], matrix[i + 1][0]);
      std::swap(matrix[i][2], matrix[i + 1][1]);
      std::swap(matrix[i][3], matrix[i + 1][2]);
      b.matrix_[i].swap(b.matrix_[i + 1]);
    }

    if (Equal(matrix[i][1], T(), epsilon_)) {
      throw std::runtime_error(
          "This algorithm cannot be applied to this matrix.");
    }

    auto multiplier = (-1) * matrix[i + 1][0] / matrix[i][1];
    matrix[i + 1][0] = 0;
    matrix[i + 1][1] += multiplier * matrix[i][2];
    matrix[i + 1][2] += multiplier * matrix[i][3];
    b.matrix_[i + 1][0] += multiplier * b.matrix_[i][0];
  }

  if (Equal(matrix[rows_ - 1][1], T(), epsilon_)) {
    throw std::runtime_error(
        "This algorithm cannot be applied to this matrix.");
  }

  // Going left and up.
  for (int i = rows_ - 1; i >= 0; --i) {
    b.matrix_[i][0] /= matrix[i][1];
    if (i > 0) {
      b.matrix_[i - 1][0] -= matrix[i - 1][2] * b.matrix_[i][0];
    }
    if (i > 1) {
      b.matrix_[i - 2][0] -= matrix[i - 2][3] * b.matrix_[i][0];
    }
  }

  return b;
}

// ---------------------------------------------------------------------------
// QR decomposition for any matrix.

template <MatrixNumber T>
void Matrix<T>::MultiplyMatrixByRotation_Left(
    std::vector<std::vector<T>>* matrix, T sin, T cos,
    int i, int j, int size) const {
  std::vector<T> new_j_row = (*matrix)[j];
  std::vector<T> new_i_row = (*matrix)[j];

  for (int k = 0; k < size; ++k) {
    new_j_row[k] *= cos;
  }
  for (int k = 0; k < size; ++k) {
    new_j_row[k] -= sin * (*matrix)[i][k];
  }
  for (int k = 0; k < size; ++k) {
    new_i_row[k] *= sin;
  }
  for (int k = 0; k < size; ++k) {
    new_i_row[k] += cos * (*matrix)[i][k];
  }

  (*matrix)[j] = std::move(new_j_row);
  (*matrix)[i] = std::move(new_i_row);
}

template <MatrixNumber T>
void Matrix<T>::MultiplyMatrixByRotation_Right(
    std::vector<std::vector<T>>* matrix, T sin, T cos,
    int i, int j, int size) const {
  std::vector<T> new_j_column = std::vector<T>(size, 0);
  std::vector<T> new_i_column = std::vector<T>(size, 0);

  for (int k = 0; k < size; ++k) {
    new_j_column[k] = cos * (*matrix)[k][j] + sin * (*matrix)[k][i];
  }
  for (int k = 0; k < size; ++k) {
    new_i_column[k] = -sin * (*matrix)[k][j] + cos * (*matrix)[k][i];
  }

  for (int k = 0; k < size; ++k) {
    (*matrix)[k][j] = new_j_column[k];
  }
  for (int k = 0; k < size; ++k) {
    (*matrix)[k][i] = new_i_column[k];
  }
}

template <MatrixNumber T>
void Matrix<T>::CountQrDecomposition() {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }
  int size = rows_;

  R_matrix_QR_ = matrix_;
  Q_matrix_QR_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    Q_matrix_QR_[i][i] = 1;
  }

  for (int j = 0; j < size - 1; ++j) {
    for (int i = j + 1; i < size; ++i) {
      // Making matrix[i][j] element equal to zero.

      if (Equal(R_matrix_QR_[i][j], T(), epsilon_)) {
        continue;
      }

      T divisor = std::sqrt(R_matrix_QR_[i][j] * R_matrix_QR_[i][j] +
          R_matrix_QR_[j][j] * R_matrix_QR_[j][j]);

      if (Equal(divisor, T(), epsilon_)) {
        Q_matrix_QR_.clear();
        R_matrix_QR_.clear();
        throw std::runtime_error(
            "Givens algorithm can't be applied to this matrix.");
      }

      T cos = R_matrix_QR_[j][j] / divisor;
      T sin = -1 * R_matrix_QR_[i][j] / divisor;

      MultiplyMatrixByRotation_Left(&Q_matrix_QR_, sin, cos, i, j, size);
      MultiplyMatrixByRotation_Left(&R_matrix_QR_, sin, cos, i, j, size);

      R_matrix_QR_[i][j] = 0;
    }
  }
}

// ---------------------------------------------------------------------------
// QR algorithm for any matrix (based on upper Hessenberg matrices).

template <MatrixNumber T>
void Matrix<T>::CountUpperHessenbergMatrix() {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }
  int size = rows_;

  hessenberg_matrix_ = matrix_;

  for (int j = 0; j < size - 2; ++j) {
    for (int i = j + 2; i < size; ++i) {
      // Making matrix[i][j] element equal to zero. Multiplying current matrix
      // from the left and from the right (so the new matrix will be similar
      // to it and will have all the eigenvalues saved).

      if (Equal(hessenberg_matrix_[i][j], T(), epsilon_)) {
        continue;
      }

      T divisor = std::sqrt(
          hessenberg_matrix_[j + 1][j] * hessenberg_matrix_[j + 1][j] +
              hessenberg_matrix_[i][j] * hessenberg_matrix_[i][j]);

      if (Equal(divisor, T(), epsilon_)) {
        hessenberg_matrix_.clear();
        throw std::runtime_error(
            "Givens algorithm can't be applied to this matrix.");
      }

      T cos = hessenberg_matrix_[j + 1][j] / divisor;
      T sin = -1 * hessenberg_matrix_[i][j] / divisor;

      MultiplyMatrixByRotation_Left(
          &hessenberg_matrix_, sin, cos, i, j + 1, size);
      MultiplyMatrixByRotation_Right(
          &hessenberg_matrix_, -sin, cos, i, j + 1, size);
      MultiplyMatrixByRotation_Right(
          &hessenberg_rotation_matrix_, -sin, cos, i, j + 1, size);

      hessenberg_matrix_[i][j] = 0;
    }
  }
}

template <MatrixNumber T>
void Matrix<T>::RunQRAlgorithmIteration() {
  int size = rows_;
  std::vector<std::tuple<T, T, int, int>> rotation_matrices_;

  for (int i = 1; i < size; ++i) {
    // Making matrix[i][i - 1] element equal to zero.

    T divisor = std::sqrt(
        hessenberg_matrix_[i][i - 1] * hessenberg_matrix_[i][i - 1] +
            hessenberg_matrix_[i - 1][i - 1]
                * hessenberg_matrix_[i - 1][i - 1]);

    if (Equal(divisor, T(), epsilon_)) {
      throw std::runtime_error(
          "This algorithm can't be applied to this matrix.");
    }

    T cos = hessenberg_matrix_[i - 1][i - 1] / divisor;
    T sin = -1 * hessenberg_matrix_[i][i - 1] / divisor;

    MultiplyMatrixByRotation_Left(&hessenberg_matrix_,
                                  sin, cos, i, i - 1, size);
    rotation_matrices_.emplace_back(sin, cos, i, i - 1);

    hessenberg_matrix_[i][i - 1] = 0;
  }

  int number_of_rotation_matrices = rotation_matrices_.size();
  for (int i = 0; i < number_of_rotation_matrices; ++i) {
    MultiplyMatrixByRotation_Right(
        &hessenberg_matrix_,
        -std::get<0>(rotation_matrices_[i]),
        std::get<1>(rotation_matrices_[i]),
        std::get<2>(rotation_matrices_[i]),
        std::get<3>(rotation_matrices_[i]),
        size);

    MultiplyMatrixByRotation_Right(
        &hessenberg_rotation_matrix_,
        -std::get<0>(rotation_matrices_[i]),
        std::get<1>(rotation_matrices_[i]),
        std::get<2>(rotation_matrices_[i]),
        std::get<3>(rotation_matrices_[i]),
        size);
  }
}

template <MatrixNumber T>
void Matrix<T>::RunQrAlgorithm(T epsilon, int max_number_of_iterations) {
  int size = rows_;

  if (hessenberg_matrix_.empty()) {
    hessenberg_rotation_matrix_ =
        std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
    for (int i = 0; i < size; ++i) {
      hessenberg_rotation_matrix_[i][i] = 1;
    }

    CountUpperHessenbergMatrix();
  }

  std::vector<T> previous_diagonal_values(size, 0);
  std::vector<T> previous_left_diagonal_values(size - 1, 0);

  bool one_more_iteration = true;
  int number_of_iterations = 0;

  while (one_more_iteration
      && number_of_iterations < max_number_of_iterations) {
    RunQRAlgorithmIteration();
    ++number_of_iterations;

    one_more_iteration = false;
    for (int i = 0; i < size; ++i) {
      if (!one_more_iteration &&
          std::abs(previous_diagonal_values[i] - hessenberg_matrix_[i][i])
              > epsilon) {
        one_more_iteration = true;
      }
      previous_diagonal_values[i] = hessenberg_matrix_[i][i];
    }
    for (int i = 1; i < size; ++i) {
      if (!one_more_iteration && std::abs(
          previous_left_diagonal_values[i - 1] - hessenberg_matrix_[i][i - 1])
          > epsilon) {
        one_more_iteration = true;
      }
      previous_left_diagonal_values[i - 1] = hessenberg_matrix_[i][i - 1];
    }
  }

  if (number_of_iterations == max_number_of_iterations) {
    std::cerr << "QR algorithm: stopped after the maximum number of iterations"
                 " reached." << std::endl;
  }
}

template <MatrixNumber T>
std::vector<typename Matrix<T>::Eigenvector>
Matrix<T>::GetEigenvectorsFromHessenbergMatrix(T epsilon) const {
  int size = rows_;
  std::vector<Eigenvector> result;

  auto get_ith_column = [&size](
      const std::vector<std::vector<T>>& matrix, int i) -> Matrix<T> {
    std::vector<std::vector<T>> result(size, std::vector<T>(1));
    for (int j = 0; j < size; ++j) {
      result[j][0] = matrix[j][i];
    }
    return Matrix(result);
  };

  for (int i = 0; i < size; ++i) {
    if (i < size - 1 && std::abs(hessenberg_matrix_[i + 1][i]) > epsilon) {
      T a = hessenberg_matrix_[i][i];
      T b = hessenberg_matrix_[i][i + 1];
      T c = hessenberg_matrix_[i + 1][i];
      T d = hessenberg_matrix_[i + 1][i + 1];

      T real_part = (a + d) / 2;
      T imaginary_part =
          std::sqrt(4 * a * d - 4 * b * c - (a + d) * (a + d)) / 2;

      result.push_back(
          {get_ith_column(hessenberg_rotation_matrix_, i),
           {real_part, imaginary_part}});
      result.push_back(
          {get_ith_column(hessenberg_rotation_matrix_, i + 1),
           {real_part, -imaginary_part}});
      ++i;
    } else {
      result.push_back(
          {get_ith_column(hessenberg_rotation_matrix_, i),
           {hessenberg_matrix_[i][i], 0}});
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// Danilevsky algorithm for any matrix.

template <MatrixNumber T>
void Matrix<T>::CountFrobeniusMatrix() {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }

  int size = rows_;
  frobenius_matrix_ = matrix_;

  frobenius_transition_matrix_ =
      std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    frobenius_transition_matrix_[i][i] = 1;
  }

  int current_row = size - 1;
  int current_column = size - 2;

  while (current_row > 0) {
    int non_zero_column = current_column;
    while (non_zero_column >= 0 && Equal(
        frobenius_matrix_[current_row][non_zero_column], T(), epsilon_)) {
      --non_zero_column;
    }
    if (non_zero_column == -1) {
      --current_row;
      current_column = current_row - 1;
      continue;
    }
    if (non_zero_column != current_column) {
      for (int i = 0; i < size; ++i) {
        std::swap(frobenius_matrix_[i][non_zero_column],
                  frobenius_matrix_[i][current_column]);
        std::swap(frobenius_transition_matrix_[i][non_zero_column],
                  frobenius_transition_matrix_[i][current_column]);
      }
      frobenius_matrix_[non_zero_column].swap(
          frobenius_matrix_[current_column]);
    }

    auto divisor = frobenius_matrix_[current_row][current_column];
    for (int i = 0; i < current_row; ++i) {
      frobenius_matrix_[i][current_column] /= divisor;
    }
    for (int i = 0; i < size; ++i) {
      frobenius_transition_matrix_[i][current_column] /= divisor;
    }
    for (int i = 0; i < size; ++i) {
      frobenius_matrix_[current_column][i] *= divisor;
    }
    frobenius_matrix_[current_row][current_column] = 1;

    for (int i = 0; i < size; ++i) {
      if (i == current_column) continue;
      auto multiplier = frobenius_matrix_[current_row][i];
      frobenius_matrix_[current_row][i] = 0;

      for (int j = 0; j < current_row; ++j) {
        frobenius_matrix_[j][i] -=
            multiplier * frobenius_matrix_[j][current_column];
      }
      for (int j = 0; j < size; ++j) {
        frobenius_transition_matrix_[j][i] -=
            multiplier * frobenius_transition_matrix_[j][current_column];
      }
      for (int j = 0; j < size; ++j) {
        frobenius_matrix_[current_column][j] +=
            multiplier * frobenius_matrix_[i][j];
      }
    }

    --current_row;
    --current_column;
  }
}

template <MatrixNumber T>
void Matrix<T>::CountCharacteristicPolynomial() {
  if (frobenius_matrix_.empty()) {
    CountFrobeniusMatrix();
  }

  characteristic_polynomial_ = Polynomial({1}, epsilon_);
  int size = rows_;

  int current_row = size - 1;
  int current_column = current_row - 1;

  while (current_row >= 0) {
    if (current_column == -1 ||
        Equal(frobenius_matrix_[current_row][current_column], T(), epsilon_)) {
      std::vector<T> current_coefficients;
      for (int i = current_column + 1; i < size; ++i) {
        current_coefficients.push_back(-frobenius_matrix_[current_row][i]);
      }
      std::reverse(current_coefficients.begin(), current_coefficients.end());
      current_coefficients.push_back(1);
      characteristic_polynomial_ =
          characteristic_polynomial_ * Polynomial(current_coefficients);
      if ((size - current_column - 1) % 2 == 1) {
        characteristic_polynomial_ =
            static_cast<T>(-1) * characteristic_polynomial_;
      }

      size = current_row;
    }
    --current_row;
    --current_column;
  }
}

template <MatrixNumber T>
void Matrix<T>::FindCharacteristicPolynomialRoots() {
  characteristic_polynomial_.FindRoots();
}

template <MatrixNumber T>
std::vector<typename Matrix<T>::Eigenvector>
Matrix<T>::GetEigenvectorsFromFrobeniusMatrix() const {
  std::vector<Eigenvector> result;
  int size = rows_;
  auto transition_matrix = Matrix(frobenius_transition_matrix_);

  for (auto eigenvalue : characteristic_polynomial_.GetRoots()) {
    std::vector<std::vector<T>> eigenvector(size, std::vector<T>(1));
    T current_eigenvalue_degree = 1;
    for (int i = size - 1; i >= 0; --i) {
      eigenvector[i][0] = current_eigenvalue_degree;
      current_eigenvalue_degree *= eigenvalue;
    }
    result.push_back({transition_matrix * Matrix(eigenvector),
                      eigenvalue});
  }

  return result;
}

// ---------------------------------------------------------------------------
// Power iteration algorithm for any matrix.

template <MatrixNumber T>
std::vector<typename Matrix<T>::Eigenvector>
Matrix<T>::GetPowerIterationResults(T epsilon,
                                    int max_number_of_iterations) const {
  if (rows_ != columns_) {
    throw std::runtime_error("Matrix is not square.");
  }

  std::vector<std::vector<T>> u_vector(rows_, std::vector<T>(1, 0));
  u_vector[0][0] = 1;

  auto vector_max_value = [](const Matrix<T>& vector) -> T {
    int size = vector.matrix_.size();
    T value = vector.matrix_[0][0];
    for (int i = 0; i < size; ++i) {
      value = (std::abs(value) < std::abs(vector.matrix_[i][0]))
              ? vector.matrix_[i][0] : value;
    }
    return value;
  };

  auto iterations =
      [this, &u_vector, vector_max_value, max_number_of_iterations, epsilon](
          bool sign) -> std::optional<Eigenvector> {
        T lambda = 0;
        Matrix u(u_vector);
        Matrix v(u_vector);
        Matrix this_multiply_v = operator*(v);

        int number_of_iterations = 0;

        while ((this_multiply_v - lambda * v).CountAndGetNorm() > epsilon) {
          v = this_multiply_v;
          u = operator*(v);

          lambda = std::sqrt(std::abs(vector_max_value(u)));
          if (sign) {
            lambda = -lambda;
          }

          v = lambda * v + u;
          v = v / vector_max_value(v);

          this_multiply_v = operator*(v);

          ++number_of_iterations;
          if (number_of_iterations > max_number_of_iterations) {
            return std::nullopt;
          }
        }

        return Eigenvector{v, lambda};
      };

  std::promise<std::optional<Eigenvector>> first_completed;
  std::promise<std::optional<Eigenvector>> second_completed;

  auto first_future = first_completed.get_future();
  auto second_future = second_completed.get_future();

  std::thread([&iterations](
      std::promise<std::optional<Eigenvector>>& first_completed) {
    auto result = iterations(false);
    first_completed.set_value(result);
  }, std::ref(first_completed)).join();

  std::thread([&iterations](
      std::promise<std::optional<Eigenvector>>& second_completed) {
    auto result = iterations(true);
    second_completed.set_value(result);
  }, std::ref(second_completed)).join();

  std::vector<Eigenvector> result = {};

  auto first_future_result = first_future.get();
  if (first_future_result.has_value()) {
    result.push_back(std::move(first_future_result.value()));
  }

  auto second_future_result = second_future.get();
  if (second_future_result.has_value()) {
    result.push_back(std::move(second_future_result.value()));
  }

  return result;
}

// ---------------------------------------------------------------------------
// Getters for the results of TLU, LDL, QR decompositions, counting the
// inverse matrix and the condition number, etc.

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetLMatrix_TLU() const {
  return Matrix(L_matrix_TLU_, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetUMatrix_TLU() const {
  return Matrix(U_matrix_TLU_, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetTMatrix_TLU() const {
  int size = rows_;
  std::vector<std::vector<T>> t_matrix(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    t_matrix[i][T_matrix_TLU_[i]] = 1;
  }
  return Matrix(t_matrix, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetTInverseMatrix_TLU() const {
  int size = rows_;
  std::vector<std::vector<T>> d_inverse_matrix(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    d_inverse_matrix[T_inverse_matrix_TLU_[i]][i] = 1;
  }
  return Matrix(d_inverse_matrix, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetLTMatrix_LDL() const {
  return Matrix(LT_matrix_LDL_, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetDMatrix_LDL() const {
  int size = rows_;
  std::vector<std::vector<T>> d_matrix(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    d_matrix[i][i] = D_matrix_LDL_[i];
  }
  return Matrix(d_matrix, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetQMatrix_QR() const {
  return Matrix(Q_matrix_QR_, epsilon_).GetTranspose();
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetRMatrix_QR() const {
  return Matrix(R_matrix_QR_, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetUpperHessenbergMatrix() const {
  return Matrix(hessenberg_matrix_, epsilon_);
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetFrobeniusMatrix() const {
  return Matrix(frobenius_matrix_, epsilon_);
}

template <MatrixNumber T>
Polynomial<T> Matrix<T>::GetCharacteristicPolynomial() const {
  return characteristic_polynomial_;
}

template <MatrixNumber T>
std::vector<T> Matrix<T>::GetCharacteristicPolynomialRealRoots() const {
  return characteristic_polynomial_.GetRoots();
}

template <MatrixNumber T>
Matrix<T> Matrix<T>::GetInverseMatrix() const {
  return Matrix(inverse_matrix_, epsilon_);
}

template <MatrixNumber T>
std::optional<T> Matrix<T>::GetConditionNumber() const {
  return condition_number_;
}

}  // namespace matrix

#endif  // MATRIX_H_

#ifndef MATRIX_H_
#define MATRIX_H_

// Base Matrix template class. Requires template arguments satisfying
// std::is_floating_point<T> condition.
// A constructor from the two-dimensional vector is provided, as well as
// comparison operators (equal and not equal), sum, difference and product
// operators. << operator is overloaded for printing matrices.
// Methods for counting TLU / LDL decomposition, solving systems of linear
// equations, counting an inverse matrix and counting a condition number are
// implemented. Their descriptions (with time asymptotics estimate) are
// located below.

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "Utils/ThreadPool.h"

namespace matrix {

using matrix_utils::ThreadPool;

template<class T>
class Matrix {
 public:
  explicit Matrix(std::vector<std::vector<T>> matrix,
                  T epsilon = std::numeric_limits<T>::epsilon());
  ~Matrix() = default;

// ---------------------------------------------------------------------------
// Comparison operators. Use predefined epsilons to compare floating point
// values properly.

  bool operator==(const Matrix& other) const;
  bool operator!=(const Matrix& other) const;

// ---------------------------------------------------------------------------
// Arithmetic operators.

  Matrix operator+(const Matrix& other) const;
  Matrix operator-(const Matrix& other) const;
  Matrix operator*(const Matrix& other) const;
  Matrix operator*(T number) const;
  template<class U>
  friend Matrix<U> operator*(U number, const Matrix<U>& matrix);

// ---------------------------------------------------------------------------
// Printing the matrix.

  [[nodiscard]] std::string ToString() const;

// ---------------------------------------------------------------------------
// Matrix norm functions. Matrix norm, induced by the vector max-norm, is used.

  void CountNorm();
  std::optional<T> GetNorm() const;

// ---------------------------------------------------------------------------
// Getters.

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
// There're two ways for counting the inverse matrix. The both have pros and
// cons, which are described in the report.
// CountInverseMatrix_AlmostTriangular_TLU counts the inverse matrix using
// the optimized algorithm for TLU decomposition.
// CountInverseMatrix_AlmostTriangular counts the inverse matrix making
// the matrix lower-triangular and using this property.

  // Checks if the matrix is satisfying Task 1 conditions (is almost triangular,
  // i.e. has only one non-zero element to the right of the diagonal).
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

  // Counts QR decomposition of the matrix (using Givens method). Internal
  // orthogonal matrices are multiplied efficiently (only two rows and two
  // columns are changed) to improve time complexity.
  // Time complexity is O(n^3).
  void CountQrDecomposition();

// ---------------------------------------------------------------------------
// Getters for the results of TLU decomposition, LDL decomposition, QR
// decomposition, counting the inverse matrix and the condition number.

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
// Results of applying different algorithms.

  std::vector<std::vector<T>> inverse_matrix_{};
  std::optional<T> condition_number_ = std::nullopt;
};

template<class T>
Matrix<T>::Matrix(std::vector<std::vector<T>> matrix, T epsilon)
    : matrix_(std::move(matrix)),
      epsilon_(epsilon) {
  static_assert(std::is_floating_point_v<T>,
                "Matrix must contain floating-point numerical values.");
  if (matrix_.empty() || matrix_[0].empty()) {
    throw std::invalid_argument("Empty matrix can't be constructed.");
  }
  if (epsilon_ <= 0) {
    throw std::invalid_argument("Epsilon must be positive.");
  }
  rows_ = matrix_.size();
  columns_ = matrix_[0].size();
}

// ---------------------------------------------------------------------------
// Comparison operators.

template<class T>
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

template<class T>
bool Matrix<T>::operator!=(const Matrix& other) const {
  return !operator==(other);
}

// ---------------------------------------------------------------------------
// Arithmetic operators.

template<class T>
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

template<class T>
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

template<class T>
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

template<class T>
Matrix<T> Matrix<T>::operator*(T number) const {
  Matrix<T> result(matrix_, epsilon_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < columns_; ++j) {
      result.matrix_[i][j] *= number;
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

template<class T>
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

template<class T>
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

template<class T>
std::optional<T> Matrix<T>::GetNorm() const {
  return norm_;
}

// ---------------------------------------------------------------------------
// Getters.

template<class T>
int Matrix<T>::GetNumberOfRows() const {
  return rows_;
}

template<class T>
int Matrix<T>::GetNumberOfColumns() const {
  return columns_;
}

template<class T>
T Matrix<T>::GetEpsilon() const {
  return epsilon_;
}

template<class T>
std::vector<std::vector<T>> Matrix<T>::GetMatrixAsVector() const {
  return matrix_;
}

// ---------------------------------------------------------------------------
// TLU decomposition, solving systems of linear equations, counting the inverse
// matrix and the condition number for ANY matrix.

template<class T>
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

    if (std::abs(U_matrix_TLU_[i][i]) < epsilon_) {
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

template<class T>
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

template<class T>
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

template<class T>
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

template<class T>
bool Matrix<T>::IsAlmostTriangular() const {
  if (rows_ != columns_) return false;
  for (int i = 0; i < rows_; ++i) {
    for (int j = i + 2; j < columns_; ++j) {
      if (matrix_[i][j] != 0) return false;
    }
  }
  return true;
}

template<class T>
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
    if (std::abs(U_matrix_TLU_[i][i]) < epsilon_) {
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

template<class T>
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

template<class T>
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

template<class T>
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

template<class T>
void Matrix<T>::CountInverseMatrix_AlmostTriangular() {
  int size = rows_;
  auto a_matrix = matrix_;
  inverse_matrix_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    inverse_matrix_[i][i] = 1;
  }

  // Making the matrix lower-triangular.
  for (int i = size - 1; i > 0; --i) {
    if (std::abs(a_matrix[i][i]) < epsilon_) {
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

  if (std::abs(a_matrix[0][0]) < epsilon_) {
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
              if (std::abs(a_matrix[i][i]) < epsilon_) {
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

template<class T>
void Matrix<T>::CountInverseMatrix_AlmostTriangular_SingleThread() {
  int size = rows_;
  auto a_matrix = matrix_;
  inverse_matrix_ = std::vector<std::vector<T>>(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    inverse_matrix_[i][i] = 1;
  }

  // Making the matrix lower-triangular.
  for (int i = size - 1; i > 0; --i) {
    if (std::abs(a_matrix[i][i]) < epsilon_) {
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

  if (std::abs(a_matrix[0][0]) < epsilon_) {
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

template<class T>
bool Matrix<T>::IsSymmetric() const {
  if (rows_ != columns_) return false;
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j <= i; ++j) {
      if (matrix_[i][j] != matrix_[j][i]) return false;
    }
  }
  return true;
}

template<class T>
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

template<class T>
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
    if (std::abs(LT_matrix_LDL_[i][i]) < epsilon_) {
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

template<class T>
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

template<class T>
bool Matrix<T>::IsTridiagonal() const {
  return (columns_ == 3);
}

template<class T>
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

template<class T>
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

    if (std::abs(matrix[i][1]) < epsilon_) {
      throw std::runtime_error(
          "This algorithm cannot be applied to this matrix.");
    }

    auto multiplier = (-1) * matrix[i + 1][0] / matrix[i][1];
    matrix[i + 1][0] = 0;
    matrix[i + 1][1] += multiplier * matrix[i][2];
    matrix[i + 1][2] += multiplier * matrix[i][3];
    b.matrix_[i + 1][0] += multiplier * b.matrix_[i][0];
  }

  if (std::abs(matrix[rows_ - 1][1]) < epsilon_) {
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

template<class T>
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

  // A function to multiply the matrix by the rotation matrix in O(n)
  // time (matrix = rotation_matrix * matrix). We also try to make compiler
  // apply vectorization here.
  auto multiply_rotation =
      [](std::vector<std::vector<T>>& matrix,
         T sin, T cos,
         int i, int j, int size) {
        std::vector<T> new_j_row = matrix[j];
        std::vector<T> new_i_row = matrix[j];

        for (int k = 0; k < size; ++k) {
          new_j_row[k] *= cos;
        }
        for (int k = 0; k < size; ++k) {
          new_j_row[k] -= sin * matrix[i][k];
        }
        for (int k = 0; k < size; ++k) {
          new_i_row[k] *= sin;
        }
        for (int k = 0; k < size; ++k) {
          new_i_row[k] += cos * matrix[i][k];
        }

        matrix[j] = std::move(new_j_row);
        matrix[i] = std::move(new_i_row);
      };

  for (int j = 0; j < size - 1; ++j) {
    for (int i = j + 1; i < size; ++i) {
      // Making matrix[i][j] element equal to zero.

      T divisor = std::sqrt(R_matrix_QR_[i][j] * R_matrix_QR_[i][j] +
          R_matrix_QR_[j][j] * R_matrix_QR_[j][j]);

      if (divisor < epsilon_) {
        throw std::runtime_error(
            "Givens algorithm can't be applied to this matrix.");
      }

      T cos = R_matrix_QR_[j][j] / divisor;
      T sin = -1 * R_matrix_QR_[i][j] / divisor;

      multiply_rotation(Q_matrix_QR_, sin, cos, i, j, size);
      multiply_rotation(R_matrix_QR_, sin, cos, i, j, size);

      R_matrix_QR_[i][j] = 0;
    }
  }
}

// ---------------------------------------------------------------------------
// Getters for the results of TLU, LDL, QR decompositions, counting the
// inverse matrix and the condition number.

template<class T>
Matrix<T> Matrix<T>::GetLMatrix_TLU() const {
  return Matrix(L_matrix_TLU_, epsilon_);
}

template<class T>
Matrix<T> Matrix<T>::GetUMatrix_TLU() const {
  return Matrix(U_matrix_TLU_, epsilon_);
}

template<class T>
Matrix<T> Matrix<T>::GetTMatrix_TLU() const {
  int size = rows_;
  std::vector<std::vector<T>> t_matrix(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    t_matrix[i][T_matrix_TLU_[i]] = 1;
  }
  return Matrix(t_matrix, epsilon_);
}

template<class T>
Matrix<T> Matrix<T>::GetTInverseMatrix_TLU() const {
  int size = rows_;
  std::vector<std::vector<T>> d_inverse_matrix(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    d_inverse_matrix[T_inverse_matrix_TLU_[i]][i] = 1;
  }
  return Matrix(d_inverse_matrix, epsilon_);
}

template<class T>
Matrix<T> Matrix<T>::GetLTMatrix_LDL() const {
  return Matrix(LT_matrix_LDL_, epsilon_);
}

template<class T>
Matrix<T> Matrix<T>::GetDMatrix_LDL() const {
  int size = rows_;
  std::vector<std::vector<T>> d_matrix(size, std::vector<T>(size, 0));
  for (int i = 0; i < size; ++i) {
    d_matrix[i][i] = D_matrix_LDL_[i];
  }
  return Matrix(d_matrix, epsilon_);
}

template<class T>
Matrix<T> Matrix<T>::GetQMatrix_QR() const {
  return Matrix(Q_matrix_QR_, epsilon_).GetTranspose();
}

template<class T>
Matrix<T> Matrix<T>::GetRMatrix_QR() const {
  return Matrix(R_matrix_QR_, epsilon_);
}

template<class T>
Matrix<T> Matrix<T>::GetInverseMatrix() const {
  return Matrix(inverse_matrix_, epsilon_);
}

template<class T>
std::optional<T> Matrix<T>::GetConditionNumber() const {
  return condition_number_;
}

}  // namespace matrix

#endif  // MATRIX_H_

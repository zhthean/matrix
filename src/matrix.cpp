#include <cstdint>
#include <stdexcept>

#include "exceptions/matrix_exceptions.h"
#include "matrix.h"

template <typename T>
Matrix<T>::Matrix(const unsigned int row_size, const unsigned int col_size)
    : m_row_size(row_size), m_col_size(col_size), m_elements(row_size, std::vector<T>(col_size, 0)) {}

template <typename T> Matrix<T>::Matrix(const std::vector<std::vector<T>> &elements) {
  unsigned int row_size = elements.size();
  if (elements.size() <= 0) {
    throw InvalidMatrixDimensionException("The row of the Matrix should not be empty.");
  }

  unsigned int col_size = elements[0].size();
  for (const auto &row : elements) {
    if (row.size() != col_size) {
      throw InvalidMatrixDimensionException("All rows must have the same number of columns.");
    }
  }

  m_col_size = col_size;
  m_row_size = row_size;
  m_elements = elements;
}

template <typename T>
template <typename U>
Matrix<T>::Matrix(const Matrix<U> &matrix)
    : m_row_size(matrix.get_row_size()), m_col_size(matrix.get_col_size()),
      m_elements(matrix.get_row_size(), std::vector<T>(matrix.get_col_size())) {

  for (unsigned int row = 0; row < matrix.get_row_size(); row++) {
    for (unsigned int col = 0; col < matrix.get_col_size(); col++) {
      m_elements[row][col] = static_cast<T>(matrix[row][col]);
    }
  }
}

template <typename T> unsigned int Matrix<T>::get_row_size() const { return m_row_size; }

template <typename T> unsigned int Matrix<T>::get_col_size() const { return m_col_size; }

template <typename T> std::vector<T> &Matrix<T>::operator[](const unsigned int index) {
  if (index >= m_row_size) {
    throw std::out_of_range(
        "The matrix has " + std::to_string(m_row_size) + " rows, but tried to access row at index " +
        std::to_string(index)
    );
  }
  return m_elements[index];
}

template <typename T> const std::vector<T> &Matrix<T>::operator[](const unsigned int index) const {
  if (index >= m_row_size) {
    throw std::out_of_range(
        "The matrix has " + std::to_string(m_row_size) + " rows, but tried to access row at index " +
        std::to_string(index)
    );
  }
  return m_elements[index];
}

template <typename T>
template <typename U>
auto Matrix<T>::operator+(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>> {
  if (m_row_size != matrix.get_row_size() || m_col_size != matrix.get_col_size()) {
    throw MatrixMismatchException(
        "A " + std::to_string(m_row_size) + "x" + std::to_string(m_col_size) + " matrix is trying to add with a " +
        std::to_string(matrix.get_row_size()) + "x" + std::to_string(matrix.get_col_size()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = m_elements[row][col] + matrix[row][col];
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto Matrix<T>::operator-(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>> {
  if (m_row_size != matrix.get_row_size() || m_col_size != matrix.get_col_size()) {
    throw MatrixMismatchException(
        "A " + std::to_string(m_row_size) + "x" + std::to_string(m_col_size) + " matrix is trying to subtract with a " +
        std::to_string(matrix.get_row_size()) + "x" + std::to_string(matrix.get_col_size()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = m_elements[row][col] - matrix[row][col];
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto Matrix<T>::operator*(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>> {
  if (m_col_size != matrix.get_row_size()) {
    throw MatrixMismatchException(
        "The first Matrix has a column size of " + std::to_string(m_col_size) +
        ", and the second Matrix has a row size of " + std::to_string(matrix.get_row_size())
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(m_row_size, matrix.get_col_size());

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < matrix.get_col_size(); col++) {
      for (unsigned int factor_row = 0; factor_row < matrix.get_row_size(); factor_row++) {
        result[row][col]    += 
            m_elements[row][factor_row] * matrix[factor_row][col];
      }
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto Matrix<T>::operator*(const U &factor) const -> Matrix<std::common_type_t<T, U>> {
  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = factor * m_elements[row][col];
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto Matrix<T>::element_wise_product(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>> {
  if (m_row_size != matrix.get_row_size() || m_col_size != matrix.get_col_size()) {
    throw MatrixMismatchException(
        "A " + std::to_string(m_row_size) + "x" + std::to_string(m_col_size) +
        " matrix is trying to element-wise product with a " + std::to_string(matrix.get_row_size()) + "x" +
        std::to_string(matrix.get_col_size()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = m_elements[row][col] * matrix[row][col];
    }
  }

  return result;
}

template <typename T> Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> transposed(m_col_size, m_row_size);

  for (unsigned int row = 0; row < m_col_size; row++) {
    for (unsigned int col = 0; col < m_row_size; col++) {
      transposed[row][col] = m_elements[col][row];
    }
  }

  return transposed;
}

template <typename T, typename U>
auto operator*(const U &factor, const Matrix<T> &matrix) -> Matrix<std::common_type_t<T, U>> {
  using CommonType = std::common_type_t<T, U>;
  Matrix<CommonType> result(matrix.get_row_size(), matrix.get_col_size());

  for (unsigned int row = 0; row < matrix.get_row_size(); row++) {
    for (unsigned int col = 0; col < matrix.get_col_size(); col++) {
      result[row][col] = factor * matrix[row][col];
    }
  }

  return result;
}

// Explicit Instantiation
template class Matrix<short>;
template class Matrix<int>;
template class Matrix<long>;
template class Matrix<long long>;
template class Matrix<float>;
template class Matrix<double>;
template class Matrix<long double>;

INSTANTIATION_MATRIX_CONSTRUCTOR(short, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(short, long double)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(int, long double)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(long, long double)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(long long, long double)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, double)
INSTANTIATION_MATRIX_CONSTRUCTOR(float, long double)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(double, long double)
INSTANTIATION_MATRIX_CONSTRUCTOR(long double, short)
INSTANTIATION_MATRIX_CONSTRUCTOR(long double, int)
INSTANTIATION_MATRIX_CONSTRUCTOR(long double, long)
INSTANTIATION_MATRIX_CONSTRUCTOR(long double, long long)
INSTANTIATION_MATRIX_CONSTRUCTOR(long double, float)
INSTANTIATION_MATRIX_CONSTRUCTOR(long double, double)

INSTANTIATION_MATRIX_FUNCTIONS(short, short)
INSTANTIATION_MATRIX_FUNCTIONS(short, int)
INSTANTIATION_MATRIX_FUNCTIONS(short, long)
INSTANTIATION_MATRIX_FUNCTIONS(short, long long)
INSTANTIATION_MATRIX_FUNCTIONS(short, float)
INSTANTIATION_MATRIX_FUNCTIONS(short, double)
INSTANTIATION_MATRIX_FUNCTIONS(short, long double)
INSTANTIATION_MATRIX_FUNCTIONS(int, short)
INSTANTIATION_MATRIX_FUNCTIONS(int, int)
INSTANTIATION_MATRIX_FUNCTIONS(int, long)
INSTANTIATION_MATRIX_FUNCTIONS(int, long long)
INSTANTIATION_MATRIX_FUNCTIONS(int, float)
INSTANTIATION_MATRIX_FUNCTIONS(int, double)
INSTANTIATION_MATRIX_FUNCTIONS(int, long double)
INSTANTIATION_MATRIX_FUNCTIONS(long, short)
INSTANTIATION_MATRIX_FUNCTIONS(long, int)
INSTANTIATION_MATRIX_FUNCTIONS(long, long)
INSTANTIATION_MATRIX_FUNCTIONS(long, long long)
INSTANTIATION_MATRIX_FUNCTIONS(long, float)
INSTANTIATION_MATRIX_FUNCTIONS(long, double)
INSTANTIATION_MATRIX_FUNCTIONS(long, long double)
INSTANTIATION_MATRIX_FUNCTIONS(long long, short)
INSTANTIATION_MATRIX_FUNCTIONS(long long, int)
INSTANTIATION_MATRIX_FUNCTIONS(long long, long)
INSTANTIATION_MATRIX_FUNCTIONS(long long, long long)
INSTANTIATION_MATRIX_FUNCTIONS(long long, float)
INSTANTIATION_MATRIX_FUNCTIONS(long long, double)
INSTANTIATION_MATRIX_FUNCTIONS(long long, long double)
INSTANTIATION_MATRIX_FUNCTIONS(float, short)
INSTANTIATION_MATRIX_FUNCTIONS(float, int)
INSTANTIATION_MATRIX_FUNCTIONS(float, long)
INSTANTIATION_MATRIX_FUNCTIONS(float, long long)
INSTANTIATION_MATRIX_FUNCTIONS(float, float)
INSTANTIATION_MATRIX_FUNCTIONS(float, double)
INSTANTIATION_MATRIX_FUNCTIONS(float, long double)
INSTANTIATION_MATRIX_FUNCTIONS(double, short)
INSTANTIATION_MATRIX_FUNCTIONS(double, int)
INSTANTIATION_MATRIX_FUNCTIONS(double, long)
INSTANTIATION_MATRIX_FUNCTIONS(double, long long)
INSTANTIATION_MATRIX_FUNCTIONS(double, float)
INSTANTIATION_MATRIX_FUNCTIONS(double, double)
INSTANTIATION_MATRIX_FUNCTIONS(double, long double)
INSTANTIATION_MATRIX_FUNCTIONS(long double, short)
INSTANTIATION_MATRIX_FUNCTIONS(long double, int)
INSTANTIATION_MATRIX_FUNCTIONS(long double, long)
INSTANTIATION_MATRIX_FUNCTIONS(long double, long long)
INSTANTIATION_MATRIX_FUNCTIONS(long double, float)
INSTANTIATION_MATRIX_FUNCTIONS(long double, double)
INSTANTIATION_MATRIX_FUNCTIONS(long double, long double)

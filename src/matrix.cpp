﻿#include <cstdint>
#include <stdexcept>

#include "comparison.h"
#include "exceptions/matrix_exceptions.h"
#include "matrix.h"
#include "utils/matrix_func.h"

template <typename T> supp_math::Matrix<T>::Matrix() : m_row_size(0), m_col_size(0) {}

template <typename T>
supp_math::Matrix<T>::Matrix(const unsigned int row_size, const unsigned int col_size)
    : m_row_size(row_size), m_col_size(col_size), m_elements(row_size, std::vector<T>(col_size, 0)) {}

template <typename T> supp_math::Matrix<T>::Matrix(const std::vector<std::vector<T>> &elements) {
  unsigned int row_size = static_cast<unsigned int>(elements.size());
  if (elements.size() <= 0) {
    throw supp_math::InvalidMatrixDimensionException("The row of the Matrix should not be empty.");
  }

  unsigned int col_size = static_cast<unsigned int>(elements[0].size());
  for (const auto &row : elements) {
    if (row.size() != col_size) {
      throw supp_math::InvalidMatrixDimensionException("All rows must have the same number of columns.");
    }
  }

  m_col_size = col_size;
  m_row_size = row_size;
  m_elements = elements;
}

template <typename T>
template <typename U>
supp_math::Matrix<T>::Matrix(const supp_math::Matrix<U> &other)
    : m_row_size(other.get_row_size()), m_col_size(other.get_col_size()),
      m_elements(other.get_row_size(), std::vector<T>(other.get_col_size())) {

  for (unsigned int row = 0; row < other.get_row_size(); row++) {
    for (unsigned int col = 0; col < other.get_col_size(); col++) {
      m_elements[row][col] = static_cast<T>(other[row][col]);
    }
  }
}

template <typename T> unsigned int supp_math::Matrix<T>::get_row_size() const { return m_row_size; }

template <typename T> unsigned int supp_math::Matrix<T>::get_col_size() const { return m_col_size; }

template <typename T> std::vector<T> &supp_math::Matrix<T>::operator[](const unsigned int index) {
  if (index >= m_row_size) {
    throw std::out_of_range(
        "The matrix has " + std::to_string(m_row_size) + " rows, but tried to access row at index " +
        std::to_string(index)
    );
  }
  return m_elements[index];
}

template <typename T> const std::vector<T> &supp_math::Matrix<T>::operator[](const unsigned int index) const {
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
auto supp_math::Matrix<T>::operator+(const supp_math::Matrix<U> &matrix) const
    -> supp_math::Matrix<std::common_type_t<T, U>> {
  if (m_row_size != matrix.get_row_size() || m_col_size != matrix.get_col_size()) {
    throw supp_math::MatrixMismatchException(
        "A " + std::to_string(m_row_size) + "x" + std::to_string(m_col_size) + " matrix is trying to add with a " +
        std::to_string(matrix.get_row_size()) + "x" + std::to_string(matrix.get_col_size()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  supp_math::Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = m_elements[row][col] + matrix[row][col];
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto supp_math::Matrix<T>::operator-(const supp_math::Matrix<U> &matrix) const
    -> supp_math::Matrix<std::common_type_t<T, U>> {
  if (m_row_size != matrix.get_row_size() || m_col_size != matrix.get_col_size()) {
    throw MatrixMismatchException(
        "A " + std::to_string(m_row_size) + "x" + std::to_string(m_col_size) + " matrix is trying to subtract with a " +
        std::to_string(matrix.get_row_size()) + "x" + std::to_string(matrix.get_col_size()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  supp_math::Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = m_elements[row][col] - matrix[row][col];
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto supp_math::Matrix<T>::operator*(const supp_math::Matrix<U> &matrix) const
    -> supp_math::Matrix<std::common_type_t<T, U>> {
  if (m_col_size != matrix.get_row_size()) {
    throw supp_math::MatrixMismatchException(
        "The first Matrix has a column size of " + std::to_string(m_col_size) +
        ", and the second Matrix has a row size of " + std::to_string(matrix.get_row_size())
    );
  }

  using CommonType = std::common_type_t<T, U>;
  supp_math::Matrix<CommonType> result(m_row_size, matrix.get_col_size());

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < matrix.get_col_size(); col++) {
      for (unsigned int factor_row = 0; factor_row < matrix.get_row_size(); factor_row++) {
        result[row][col] += m_elements[row][factor_row] * matrix[factor_row][col];
      }
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto supp_math::Matrix<T>::operator*(const U &factor) const -> supp_math::Matrix<std::common_type_t<T, U>> {
  using CommonType = std::common_type_t<T, U>;
  supp_math::Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = factor * m_elements[row][col];
    }
  }

  return result;
}

template <typename T>
template <typename U>
auto supp_math::Matrix<T>::element_wise_product(const supp_math::Matrix<U> &matrix) const
    -> supp_math::Matrix<std::common_type_t<T, U>> {
  if (m_row_size != matrix.get_row_size() || m_col_size != matrix.get_col_size()) {
    throw supp_math::MatrixMismatchException(
        "A " + std::to_string(m_row_size) + "x" + std::to_string(m_col_size) +
        " matrix is trying to element-wise product with a " + std::to_string(matrix.get_row_size()) + "x" +
        std::to_string(matrix.get_col_size()) + " matrix"
    );
  }

  using CommonType = std::common_type_t<T, U>;
  supp_math::Matrix<CommonType> result(m_row_size, m_col_size);

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      result[row][col] = m_elements[row][col] * matrix[row][col];
    }
  }

  return result;
}

template <typename T> supp_math::Matrix<T> supp_math::Matrix<T>::transpose() const {
  supp_math::Matrix<T> transposed(m_col_size, m_row_size);

  for (unsigned int row = 0; row < m_col_size; row++) {
    for (unsigned int col = 0; col < m_row_size; col++) {
      transposed[row][col] = m_elements[col][row];
    }
  }

  return transposed;
}

template <typename T> double supp_math::Matrix<T>::determinant() const {
  if (m_row_size != m_col_size) {
    throw supp_math::MatrixNonSquareException("The matrix is not a square matrix.");
  }

  if (m_row_size == 1) {
    return static_cast<double>(m_elements[0][0]);
  }

  if (m_row_size == 2) {
    return static_cast<double>(m_elements[0][0] * m_elements[1][1] - m_elements[0][1] * m_elements[1][0]);
  }

  if (m_row_size == 3) {
    return static_cast<double>(
        m_elements[0][0] * (m_elements[1][1] * m_elements[2][2] - m_elements[1][2] * m_elements[2][1]) -
        m_elements[0][1] * (m_elements[1][0] * m_elements[2][2] - m_elements[1][2] * m_elements[2][0]) +
        m_elements[0][2] * (m_elements[1][0] * m_elements[2][1] - m_elements[1][1] * m_elements[2][0])
    );
  }

  try {
    auto [eliminated_matrix, num_swaps] = gaussian_elimination(*this);
    std::cout << eliminated_matrix << std::endl;
    double det = 1;
    for (unsigned int pivot = 0; pivot < m_row_size; pivot++) {
      det *= eliminated_matrix[pivot][pivot];
    }

    return std::pow(-1, num_swaps) * det;
  } catch (supp_math::MatrixSingularException &) {
    return 0;
  }
}

template <typename T> supp_math::Matrix<double> supp_math::Matrix<T>::inverse() const {
  if (m_row_size != m_col_size) {
    throw supp_math::MatrixNonSquareException("The matrix is not a square matrix.");
  }

  const double det = determinant();

  if (det == 0) {
    throw supp_math::MatrixNonInvertibleException("The matrix is singular and hence is not invertible.");
  }

  supp_math::Matrix<double> inversed(m_row_size, m_col_size);

  if (m_row_size == 1) {
    inversed[0][0] = 1.0 / static_cast<double>(m_elements[0][0]);
    return inversed;
  }

  if (m_row_size == 2) {
    inversed[0][0] = m_elements[1][1] / det;
    inversed[0][1] = -m_elements[0][1] / det;
    inversed[1][0] = -m_elements[1][0] / det;
    inversed[1][1] = m_elements[0][0] / det;
    return inversed;
  }

  return gaussian_elimination_inverse(*this);
}

template <typename T>
void supp_math::Matrix<T>::swap(const unsigned int chosen_index, const unsigned int swapped_index, const char axis) {
  if (axis == 'r') {
    if (chosen_index >= m_row_size || swapped_index >= m_row_size) {
      throw std::out_of_range(
          "The matrix has " + std::to_string(m_row_size) + " rows, but tried to swap row at index " +
          std::to_string(chosen_index) + " and " + std::to_string(swapped_index)
      );
    }
    m_elements[chosen_index].swap(m_elements[swapped_index]);
  } else if (axis == 'c') {
    if (chosen_index >= m_col_size || swapped_index >= m_col_size) {
      throw std::out_of_range(
          "The matrix has " + std::to_string(m_col_size) + " columns, but tried to swap column at index " +
          std::to_string(chosen_index) + " and " + std::to_string(swapped_index)
      );
    }

    for (unsigned int row = 0; row < m_row_size; row++) {
      std::swap(m_elements[row][chosen_index], m_elements[row][swapped_index]);
    }
  } else {
    throw std::invalid_argument("The axis should be either 'r'(row) or 'c'(column).");
  }
}

template <typename T>
template <typename U>
bool supp_math::Matrix<T>::operator==(const supp_math::Matrix<U> &matrix) const {
  if (m_row_size != matrix.get_row_size() || m_col_size != matrix.get_col_size()) {
    return false;
  }

  double epsilon = std::is_same<T, float>::value ? 1e-5 : 1e-9;

  for (unsigned int row = 0; row < m_row_size; row++) {
    for (unsigned int col = 0; col < m_col_size; col++) {
      if (!numeric_almost_equal(matrix[row][col], m_elements[row][col], epsilon)) {
        return false;
      }
    }
  }

  return true;
}

template <typename T, typename U>
auto supp_math::operator*(const U &factor, const supp_math::Matrix<T> &matrix)
    -> supp_math::Matrix<std::common_type_t<T, U>> {
  using CommonType = std::common_type_t<T, U>;
  supp_math::Matrix<CommonType> result(matrix.get_row_size(), matrix.get_col_size());

  for (unsigned int row = 0; row < matrix.get_row_size(); row++) {
    for (unsigned int col = 0; col < matrix.get_col_size(); col++) {
      result[row][col] = factor * matrix[row][col];
    }
  }

  return result;
}

// Explicit Instantiation
template class supp_math::Matrix<short>;
template class supp_math::Matrix<int>;
template class supp_math::Matrix<long>;
template class supp_math::Matrix<long long>;
template class supp_math::Matrix<float>;
template class supp_math::Matrix<double>;
template class supp_math::Matrix<long double>;

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

// template std::ostream &operator<<(std::ostream&, const Matrix<short> &);
// template std::ostream &operator<<(std::ostream&, const Matrix<int> &);
// template std::ostream &operator<<(std::ostream &, const Matrix<long> &);
// template std::ostream &operator<<(std::ostream &, const Matrix<long long> &);
// template std::ostream &operator<<(std::ostream &, const Matrix<float> &);
// template std::ostream &operator<<(std::ostream &, const Matrix<double> &);
// template std::ostream &operator<<(std::ostream &, const Matrix<long double> &);

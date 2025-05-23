﻿#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "matrix_enums.h"
#include "types.h"

#define INSTANTIATION_MATRIX_CONSTRUCTOR(Type1, Type2) template Matrix<Type1>::Matrix(const Matrix<Type2> &);

#define INSTANTIATION_MATRIX_FUNCTIONS(Type1, Type2)                                                                   \
  template auto supmath::Matrix<Type1>::operator+(const Matrix<Type2> &) const                                         \
      -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
  template auto supmath::Matrix<Type1>::operator-(const Matrix<Type2> &) const                                         \
      -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
  template auto supmath::Matrix<Type1>::operator*(const Matrix<Type2> &) const                                         \
      -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
  template auto supmath::Matrix<Type1>::operator*(const Type2) const -> Matrix<std::common_type_t<Type1, Type2>>;      \
  template auto supmath::Matrix<Type1>::elementWiseProduct(const Matrix<Type2> &) const                                \
      -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
  template bool supmath::Matrix<Type1>::operator==(const Matrix<Type2> &) const;                                       \
  template auto supmath::operator*(const Type1, const Matrix<Type2> &) -> Matrix<std::common_type_t<Type1, Type2>>;

namespace supmath
{
template<Numerical T = double>
class Matrix
{
private:
  size_t                      _row_size;
  size_t                      _col_size;
  std::vector<std::vector<T>> _elements;

public:
  Matrix();

  Matrix(const Matrix<T> &other)                   = default;
  Matrix<T> &operator=(const Matrix<T> &other)     = default;
  Matrix(Matrix<T> &&other) noexcept               = default;
  Matrix<T> &operator=(Matrix<T> &&other) noexcept = default;

  ~Matrix()                                        = default;

  explicit Matrix(const size_t row_size, const size_t col_size);
  explicit Matrix(const std::vector<std::vector<T>> &elements);

  template<Numerical U>
  Matrix(const Matrix<U> &other);

  inline size_t getRowSize() const
  {
    return _row_size;
  };
  inline size_t getColSize() const
  {
    return _col_size;
  };

  std::vector<T>       &operator[](const size_t index);
  const std::vector<T> &operator[](const size_t index) const;

  template<Numerical U>
  auto operator+(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  template<Numerical U>
  auto operator-(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  template<Numerical U>
  auto operator*(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  template<Numerical U>
  auto operator*(const U factor) const -> Matrix<std::common_type_t<T, U>>;

  template<Numerical U>
  auto           elementWiseProduct(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  Matrix<T>      transpose() const;

  double         determinant() const;
  Matrix<double> inverse() const;

  void           swap(const size_t chosen_index, const size_t swapped_index, supmath::MatrixOrder order);

  template<Numerical U>
  bool operator==(const Matrix<U> &matrix) const;

  operator std::string() const;
};

template<Numerical T, Numerical U>
auto operator*(const T factor, const Matrix<U> &matrix) -> Matrix<std::common_type_t<T, U>>;

template<Numerical T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix);
}    // namespace supmath

#endif    // MATRIX_H

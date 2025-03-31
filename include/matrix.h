#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#define INSTANTIATION_MATRIX_CONSTRUCTOR(Type1, Type2)                                                                 \
  template supp_math::Matrix<Type1>::Matrix(const supp_math::Matrix<Type2> &);

#define INSTANTIATION_MATRIX_FUNCTIONS(Type1, Type2)                                                                   \
  template auto supp_math::Matrix<Type1>::operator+(const supp_math::Matrix<Type2> &) const                            \
      -> supp_math::Matrix<std::common_type_t<Type1, Type2>>;                                                          \
  template auto supp_math::Matrix<Type1>::operator-(const supp_math::Matrix<Type2> &) const                            \
      -> supp_math::Matrix<std::common_type_t<Type1, Type2>>;                                                          \
  template auto supp_math::Matrix<Type1>::operator*(const supp_math::Matrix<Type2> &) const                            \
      -> supp_math::Matrix<std::common_type_t<Type1, Type2>>;                                                          \
  template auto supp_math::Matrix<Type1>::operator*(const Type2 &) const                                               \
      -> supp_math::Matrix<std::common_type_t<Type1, Type2>>;                                                          \
  template auto supp_math::Matrix<Type1>::element_wise_product(const supp_math::Matrix<Type2> &) const                 \
      -> supp_math::Matrix<std::common_type_t<Type1, Type2>>;                                                          \
  template auto supp_math::operator*(const Type1 &, const supp_math::Matrix<Type2> &)                                  \
      -> supp_math::Matrix<std::common_type_t<Type1, Type2>>;                                                          \
  template bool supp_math::Matrix<Type1>::operator==(const supp_math::Matrix<Type2> &) const;

namespace supp_math {
template <typename T = double> class Matrix {
private:
  unsigned int m_row_size;
  unsigned int m_col_size;
  std::vector<std::vector<T>> m_elements;

public:
  Matrix();

  Matrix(const Matrix<T> &other) = default;
  Matrix<T> &operator=(const Matrix<T> &other) = default;
  Matrix(Matrix<T> &&other) noexcept = default;
  Matrix<T> &operator=(Matrix<T> &&other) noexcept = default;

  ~Matrix() = default;

  explicit Matrix(const unsigned int row_size, const unsigned int col_size);
  explicit Matrix(const std::vector<std::vector<T>> &elements);

  template <typename U> Matrix(const Matrix<U> &other);

  unsigned int get_row_size() const;
  unsigned int get_col_size() const;

  std::vector<T> &operator[](const unsigned int index);
  const std::vector<T> &operator[](const unsigned int index) const;

  template <typename U> auto operator+(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  template <typename U> auto operator-(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  template <typename U> auto operator*(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  template <typename U> auto operator*(const U &factor) const -> Matrix<std::common_type_t<T, U>>;

  template <typename U> auto element_wise_product(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

  Matrix<T> transpose() const;

  double determinant() const;
  Matrix<double> inverse() const;

  void swap(const unsigned int chosen_index, const unsigned int swapped_index, const char axis = "r");

  template <typename U> bool operator==(const Matrix<U> &matrix) const;
};

template <typename T, typename U>
auto operator*(const U &l_factor, const supp_math::Matrix<T> &r_factor) -> supp_math::Matrix<std::common_type_t<T, U>>;

template <typename T> std::ostream &operator<<(std::ostream &os, const supp_math::Matrix<T> &matrix) {
  std::string message = "Matrix([";

  for (unsigned int row = 0; row < matrix.get_row_size(); row++) {
    message += "[";
    for (unsigned int col = 0; col < matrix.get_col_size(); col++) {
      message += std::to_string(matrix[row][col]);
      if (col != matrix.get_col_size() - 1) {
        message += ", ";
      }
    }

    message += row == matrix.get_row_size() - 1 ? "]" : "], ";
  }

  message += "])";
  os << message;

  return os;
}
} // namespace supp_math

#endif // MATRIX_H

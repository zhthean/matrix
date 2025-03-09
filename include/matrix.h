#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <type_traits>
#include <vector>

#define INSTANTIATION_MATRIX_CONSTRUCTOR(Type1, Type2) template Matrix<Type1>::Matrix(const Matrix<Type2> &);

#define INSTANTIATION_MATRIX_FUNCTIONS(Type1, Type2)                                                                   \
  template auto Matrix<Type1>::operator+(const Matrix<Type2> &) const -> Matrix<std::common_type_t<Type1, Type2>>;     \
  template auto Matrix<Type1>::operator-(const Matrix<Type2> &) const -> Matrix<std::common_type_t<Type1, Type2>>;     \
  template auto Matrix<Type1>::operator*(const Matrix<Type2> &) const -> Matrix<std::common_type_t<Type1, Type2>>;     \
  template auto Matrix<Type1>::operator*(const Type2 &) const -> Matrix<std::common_type_t<Type1, Type2>>;             \
  template auto Matrix<Type1>::element_wise_product(const Matrix<Type2> &) const                                       \
      -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
  template auto operator*(const Type1 &, const Matrix<Type2> &)->Matrix<std::common_type_t<Type1, Type2>>;

template <typename T = double> class Matrix {
private:
  unsigned int m_row_size;
  unsigned int m_col_size;
  std::vector<std::vector<T>> m_elements;

public:
  explicit Matrix(const unsigned int row_size, const unsigned int col_size);
  explicit Matrix(const std::vector<std::vector<T>> &elements);

  template <typename U> Matrix(const Matrix<U> &matrix);

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

  // Matrix<double> inverse();
  // double determinant();

  void p() {
    for (auto &row : this->m_elements) {
      for (auto &col : row) {
        std::cout << col << " ";
      }
      std::cout << "\n";
    }
  }
};

template <typename T, typename U>
auto operator*(const U &l_factor, const Matrix<T> &r_factor) -> Matrix<std::common_type_t<T, U>>;

#endif // MATRIX_H

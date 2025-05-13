#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
#include "types.h"

#define INSTANTIATION_MATRIX_CONSTRUCTOR(Type1, Type2) template Matrix<Type1>::Matrix(const Matrix<Type2> &);

#define INSTANTIATION_MATRIX_FUNCTIONS(Type1, Type2)                                                                   \
	template auto supmath::Matrix<Type1>::operator+(const Matrix<Type2> &) const                                         \
	    -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
	template auto supmath::Matrix<Type1>::operator-(const Matrix<Type2> &) const                                         \
	    -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
	template auto supmath::Matrix<Type1>::operator*(const Matrix<Type2> &) const                                         \
	    -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
	template auto supmath::Matrix<Type1>::operator*(const Type2 &) const -> Matrix<std::common_type_t<Type1, Type2>>;    \
	template auto supmath::Matrix<Type1>::element_wise_product(const Matrix<Type2> &) const                              \
	    -> Matrix<std::common_type_t<Type1, Type2>>;                                                                     \
	template bool supmath::Matrix<Type1>::operator==(const Matrix<Type2> &) const;                                       \
	template auto supmath::operator*(const Type1 &factor, const Matrix<Type2> &matrix)                                   \
	    -> Matrix<std::common_type_t<Type1, Type2>>;

namespace supmath
{
template<Numerical T = double>
class Matrix
{
private:
	unsigned int                m_row_size;
	unsigned int                m_col_size;
	std::vector<std::vector<T>> m_elements;

public:
	Matrix();

	Matrix(const Matrix<T> &other)                   = default;
	Matrix<T> &operator=(const Matrix<T> &other)     = default;
	Matrix(Matrix<T> &&other) noexcept               = default;
	Matrix<T> &operator=(Matrix<T> &&other) noexcept = default;

	~Matrix()                                        = default;

	explicit Matrix(const unsigned int row_size, const unsigned int col_size);
	explicit Matrix(const std::vector<std::vector<T>> &elements);

	template<Numerical U>
	Matrix(const Matrix<U> &other);

	unsigned int          get_row_size() const;
	unsigned int          get_col_size() const;

	std::vector<T>       &operator[](const unsigned int index);
	const std::vector<T> &operator[](const unsigned int index) const;

	template<Numerical U>
	auto operator+(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

	template<Numerical U>
	auto operator-(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

	template<Numerical U>
	auto operator*(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

	template<Numerical U>
	auto operator*(const U &factor) const -> Matrix<std::common_type_t<T, U>>;

	template<Numerical U>
	auto           element_wise_product(const Matrix<U> &matrix) const -> Matrix<std::common_type_t<T, U>>;

	Matrix<T>      transpose() const;

	double         determinant() const;
	Matrix<double> inverse() const;

	void           swap(const unsigned int chosen_index, const unsigned int swapped_index, const char axis = 'r');

	template<Numerical U>
	bool operator==(const Matrix<U> &matrix) const;
};

template<Numerical T, Numerical U>
auto operator*(const U &factor, const Matrix<T> &matrix) -> Matrix<std::common_type_t<T, U>>;

template<Numerical T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix);
}    // namespace supmath

#endif    // MATRIX_H

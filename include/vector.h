#ifndef VECTOR_H
#define VECTOR_H

#include <initializer_list>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "matrix.h"
#include "types.h"

#define INSTANTIATION_VECTOR_CONSTRUCTOR(Type1, Type2)                                                                 \
  template supmath::Vector<Type1>::Vector(const Vector<Type2> &);                                                      \
  template supmath::Vector<Type1>::Vector(const Vector<Type2> &, double);

#define INSTANTIATION_VECTOR_FUNCTIONS(Type1, Type2)                                                                   \
  template auto supmath::Vector<Type1>::operator+(const Vector<Type2> &) const                                         \
      -> Vector<std::common_type_t<Type1, Type2>>;                                                                     \
  template auto supmath::Vector<Type1>::operator-(const Vector<Type2> &) const                                         \
      -> Vector<std::common_type_t<Type1, Type2>>;                                                                     \
  template auto supmath::Vector<Type1>::operator*(const Type2) const -> Vector<std::common_type_t<Type1, Type2>>;      \
  template auto supmath::Vector<Type1>::operator*(const Matrix<Type2> &) const                                         \
      -> Vector<std::common_type_t<Type1, Type2>>;                                                                     \
  template auto supmath::Vector<Type1>::crossProduct(const Vector<Type2> &) const                                      \
      -> Vector<std::common_type_t<Type1, Type2>>;                                                                     \
  template double         supmath::Vector<Type1>::dotProduct(const Vector<Type2> &) const;                             \
  template double         supmath::Vector<Type1>::angleBetween(const Vector<Type2> &) const;                           \
  template double         supmath::Vector<Type1>::distanceBetween(const Vector<Type2> &) const;                        \
  template Vector<double> supmath::Vector<Type1>::projectionOnto(const Vector<Type2> &) const;                         \
  template bool           supmath::Vector<Type1>::operator==(const Vector<Type2> &) const;                             \
  template auto supmath::operator*(const Type1, const Vector<Type2> &) -> Vector<std::common_type_t<Type1, Type2>>;    \
  template auto supmath::operator*(const Matrix<Type1> &, const Vector<Type2> &)                                       \
      -> Vector<std::common_type_t<Type1, Type2>>;

namespace supmath
{
template<Numerical T = double>
class Vector
{
private:
  std::vector<T> _elements;

public:
  Vector() = default;
  explicit Vector(size_t size);
  explicit Vector(const std::initializer_list<T> &list);

  template<Numerical U>
  explicit Vector(const Vector<U> &other);

  template<Numerical U>
  explicit Vector(const Vector<U> &other, double number);

  Vector(const Vector<T> &other)            = default;
  Vector<T> &operator=(const Vector &other) = default;
  Vector(Vector<T> &&other)                 = default;
  Vector<T> &operator=(Vector &&other)      = default;

  ~Vector()                                 = default;

  inline size_t getSize() const
  {
    return _elements.size();
  }

  inline const std::vector<T> getElements() const
  {
    return _elements;
  }

  T       &operator[](const int index);
  const T &operator[](const int index) const;

  template<Numerical U>
  auto operator+(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>>;

  template<Numerical U>
  auto operator-(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>>;

  template<Numerical U>
  auto operator*(const U factor) const -> Vector<std::common_type_t<T, U>>;

  template<Numerical U>
  auto operator*(const Matrix<U> &matrix) const -> Vector<std::common_type_t<T, U>>;

  template<Numerical U>
  auto crossProduct(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>>;

  template<Numerical U>
  double dotProduct(const Vector<U> &other) const;

  template<Numerical U>
  double angleBetween(const Vector<U> &other) const;

  template<Numerical U>
  double distanceBetween(const Vector<U> &other) const;

  template<Numerical U>
  Vector<double> projectionOnto(const Vector<U> &other) const;

  double         magnitude() const;

  Vector<double> unitVector() const;

  template<Numerical U>
  bool operator==(const Vector<U> &other) const;

  operator std::string() const;
};

template<Numerical T, Numerical U>
auto operator*(const T factor, const Vector<U> &vector) -> Vector<std::common_type_t<T, U>>;

template<Numerical T, Numerical U>
auto operator*(const Matrix<T> &matrix, const Vector<U> &vector) -> Vector<std::common_type_t<T, U>>;

template<Numerical T>
std::ostream &operator<<(std::ostream &os, const Vector<T> &vector);

};    // namespace supmath

#endif    // VECTOR_H

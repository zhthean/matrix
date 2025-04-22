#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "comparison.h"
#include "exceptions/vector_exceptions.h"
#include "types.h"

namespace supmath {
template <Numerical T = double> class Vector {
private:
  std::vector<T> _elements;

  template <Numerical U>
  static auto _multiplyVectorScalar(const U factor, const Vector<T> &vector) -> Vector<std::common_type_t<T, U>> {
    using CommonType = std::common_type_t<T, U>;
    Vector<CommonType> result(vector.getSize());

    for (size_t index = 0; index < vector.getSize(); index++) {
      result[index] = vector[index] * factor;
    }

    return result;
  }

public:
  Vector() = default;
  explicit Vector(size_t size) : _elements(size, 0) {};
  explicit Vector(const std::vector<T> &vector) : _elements(vector) {};

  template <Numerical U, Numerical V>
  explicit Vector(const Vector<U> &other, V number) : _elements(other.getSize(), 0) {
    for (size_t index = 0; index < other.getSize(); index++) {
      _elements[index] = static_cast<T>(other[index]);
    }

    _elements.push_back(static_cast<T>(number));
  };

  Vector(const Vector<T> &other) = default;
  Vector<T> &operator=(const Vector &other) = default;
  Vector(Vector<T> &&other) = default;
  Vector<T> &operator=(Vector &&other) = default;

  ~Vector() = default;

  inline size_t getSize() const { return _elements.size(); }

  inline const std::vector<T> getElements() const { return _elements; }

  T &operator[](const int index) {
    const int actual_index = (index >= 0) ? index : getSize() + index;
    if (actual_index < 0 || actual_index >= getSize()) {
      throw std::out_of_range(
          "The vector has size of " + std::to_string(getSize()) + ", but trying to access element at index " +
          std::to_string(actual_index)
      );
    }

    return _elements[actual_index];
  }
  const T &operator[](const int index) const {
    const int actual_index = (index >= 0) ? index : getSize() + index;
    if (actual_index < 0 || actual_index >= getSize()) {
      throw std::out_of_range(
          "The vector has size of " + std::to_string(getSize()) + ", but trying to access element at index " +
          std::to_string(actual_index)
      );
    }

    return _elements[actual_index];
  }

  template <Numerical U> auto operator+(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>> {
    if (getSize() != other.getSize()) {
      throw VectorInconsistentDimensionException(
          "A vector of size of " + std::to_string(getSize()) + " is trying to add with a vector of size of " +
          std::to_string(other.getSize())
      );
    }

    using CommonType = std::common_type_t<T, U>;
    Vector<CommonType> result(getSize());

    for (size_t index = 0; index < getSize(); index++) {
      result[index] = _elements[index] + other[index];
    }

    return result;
  }

  template <Numerical U> auto operator-(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>> {
    if (getSize() != other.getSize()) {
      throw VectorInconsistentDimensionException(
          "A vector of size of " + std::to_string(getSize()) + " is trying to subtract with a vector of size of " +
          std::to_string(other.getSize())
      );
    }

    using CommonType = std::common_type_t<T, U>;
    Vector<CommonType> result(getSize());

    for (size_t index = 0; index < getSize(); index++) {
      result[index] = _elements[index] - other[index];
    }

    return result;
  }

  template <Numerical U> auto operator*(const U factor) const -> Vector<std::common_type_t<T, U>> {
    return _multiplyVectorScalar(factor, *this);
  }

  template <Numerical U> auto crossProduct(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>> {
    if (getSize() != 3 && other.getSize() != 3) {
      throw UnsupportVectorSizeException("Cross product methods only support for vector of size of 3");
    }

    using CommonType = std::common_type_t<T, U>;
    Vector<CommonType> result(getSize());

    result[0] = _elements[1] * other[2] - _elements[2] * other[1];
    result[1] = _elements[2] * other[0] - _elements[0] * other[2];
    result[2] = _elements[0] * other[1] - _elements[1] * other[0];

    return result;
  }

  template <Numerical U> double dotProduct(const Vector<U> &other) const {
    if (getSize() != other.getSize()) {
      throw VectorInconsistentDimensionException(
          "A vector of size of " + std::to_string(getSize()) + " is trying to dot product with a vector of size of " +
          std::to_string(other.getSize())
      );
    }

    double result = 0;
    for (size_t index = 0; index < getSize(); index++) {
      result += _elements[index] * other[index];
    }

    return result;
  }

  template <Numerical U> double angleBetween(const Vector<U> &other) const {
    auto mag = magnitude();
    auto other_mag = other.magnitude();
    if (mag == 0 || other_mag == 0) {
      throw ZeroMagnitudeException("Cannot compute angle between with a magnitude of zero");
    }

    return std::acos(dotProduct(other) / (mag * other_mag));
  }

  template <Numerical U> double distanceBetween(const Vector<U> &other) const {
    if (getSize() != other.getSize()) {
      throw VectorInconsistentDimensionException(
          "Cannot compute distance between because one vector has size of " + std::to_string(getSize()) +
          " and the other has size of " + std::to_string(other.getSize())
      );
    }

    double sum_of_diff = 0;
    for (int index = 0; index < getSize(); index++) {
      auto diff = _elements[index] - other[index];
      sum_of_diff += diff * diff;
    }

    return std::sqrt(sum_of_diff);
  }

  template <Numerical U> auto projectionOnto(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>> {
    auto other_mag = other.magnitude();
    if (other_mag == 0) {
      throw ZeroMagnitudeException("Cannot project onto a vector with zero magnitude");
    }

    return (dotProduct(other) / (other_mag * other_mag)) * other;
  }

  double magnitude() const {
    double sum_of_square = 0;
    for (auto &element : _elements) {
      sum_of_square += element * element;
    }

    return std::sqrt(sum_of_square);
  }

  Vector<double> unitVector() const {
    auto mag = magnitude();
    if (mag == 0) {
      throw ZeroMagnitudeException("Cannot compute unit vector with a magnitude of zero");
    }

    return *this * (1.0 / magnitude());
  }

  template <Numerical U> bool operator==(const Vector<U> &other) const {
    if (getSize() != other.getSize()) {
      return false;
    }

    double epsilon = std::is_same<T, float>::value ? 1e-5 : 1e-9;

    for (size_t index = 0; index < getSize(); index++) {
      if (!numeric_almost_equal(_elements[index], other[index])) {
        return false;
      }
    }

    return true;
  }

  template <Numerical U>
  friend auto operator*(const U factor, const Vector<T> &vector) -> Vector<std::common_type_t<T, U>> {
    return _multiplyVectorScalar(factor, vector);
  }

  friend std::ostream &operator<<(std::ostream &os, const Vector<T> &vector) {
    std::string message = "Vector([";

    for (size_t index = 0; index < vector.getSize(); index++) {
      message += std::to_string(vector[index]);
      if (index != vector.getSize() - 1) {
        message += ", ";
      }
    }

    message += "])";
    os << message;

    return os;
  }
};
}; // namespace supmath

#endif // VECTOR_H

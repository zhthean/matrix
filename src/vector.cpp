#include <cmath>
#include <string>

#include "exceptions/vector_exceptions.h"
#include "vector.h"

#include "comparison.h"
#include "vector_funcs.h"

using namespace supmath;

template<Numerical T>
supmath::Vector<T>::Vector(size_t size) : _elements(size, 0)
{}

template<Numerical T>
supmath::Vector<T>::Vector(const std::initializer_list<T> &list) : _elements(list)
{}

template<Numerical T>
template<Numerical U>
supmath::Vector<T>::Vector(const Vector<U> &other) : _elements(other.getSize(), 0)
{
  for ( size_t index = 0; index < other.getSize(); index++ ) {
    _elements[index] = static_cast<T>(other[index]);
  }
};

template<Numerical T>
template<Numerical U>
supmath::Vector<T>::Vector(const Vector<U> &other, double number) : _elements(other.getSize(), 0)
{
  for ( size_t index = 0; index < other.getSize(); index++ ) {
    _elements[index] = static_cast<T>(other[index]);
  }

  _elements.push_back(static_cast<T>(number));
};

template<Numerical T>
T &supmath::Vector<T>::operator[](const int index)
{
  size_t actual_index = (index >= 0) ? index : getSize() + index;
  if ( actual_index < 0 || actual_index >= getSize() ) {
    throw std::out_of_range(
        "The vector has size of " + std::to_string(getSize()) + ", but trying to access element at index "
        + std::to_string(actual_index)
    );
  }

  return _elements[actual_index];
}

template<Numerical T>
const T &supmath::Vector<T>::operator[](const int index) const
{
  size_t actual_index = (index >= 0) ? index : getSize() + index;
  if ( actual_index < 0 || actual_index >= getSize() ) {
    throw std::out_of_range(
        "The vector has size of " + std::to_string(getSize()) + ", but trying to access element at index "
        + std::to_string(actual_index)
    );
  }

  return _elements[actual_index];
}

template<Numerical T>
template<Numerical U>
auto supmath::Vector<T>::operator+(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>>
{
  if ( getSize() != other.getSize() ) {
    throw VectorInconsistentDimensionException(
        "A vector of size of " + std::to_string(getSize()) + " is trying to add with a vector of size of "
        + std::to_string(other.getSize())
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Vector<CommonType> result(getSize());

  for ( size_t index = 0; index < getSize(); index++ ) {
    result[index] = _elements[index] + other[index];
  }

  return result;
}

template<Numerical T>
template<Numerical U>
auto supmath::Vector<T>::operator-(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>>
{
  if ( getSize() != other.getSize() ) {
    throw VectorInconsistentDimensionException(
        "A vector of size of " + std::to_string(getSize()) + " is trying to subtract with a vector of size of "
        + std::to_string(other.getSize())
    );
  }

  using CommonType = std::common_type_t<T, U>;
  Vector<CommonType> result(getSize());

  for ( size_t index = 0; index < getSize(); index++ ) {
    result[index] = _elements[index] - other[index];
  }

  return result;
}

template<Numerical T>
template<Numerical U>
auto supmath::Vector<T>::operator*(const U factor) const -> Vector<std::common_type_t<T, U>>
{
  return _multiplyVectorScalar(factor, *this);
}

template<Numerical T>
template<Numerical U>
auto supmath::Vector<T>::crossProduct(const Vector<U> &other) const -> Vector<std::common_type_t<T, U>>
{
  if ( getSize() != 3 && other.getSize() != 3 ) {
    throw UnsupportVectorSizeException("Cross product is only defined for 3-dimensional vectors");
  }

  using CommonType = std::common_type_t<T, U>;
  Vector<CommonType> result(getSize());

  result[0] = _elements[1] * other[2] - _elements[2] * other[1];
  result[1] = _elements[2] * other[0] - _elements[0] * other[2];
  result[2] = _elements[0] * other[1] - _elements[1] * other[0];

  return result;
}

template<Numerical T>
template<Numerical U>
double supmath::Vector<T>::dotProduct(const Vector<U> &other) const
{
  if ( getSize() != other.getSize() ) {
    throw VectorInconsistentDimensionException(
        "A vector of size of " + std::to_string(getSize()) + " is trying to dot product with a vector of size of "
        + std::to_string(other.getSize())
    );
  }

  double result = 0;
  for ( size_t index = 0; index < getSize(); index++ ) {
    result += _elements[index] * other[index];
  }

  return result;
}

template<Numerical T>
template<Numerical U>
double supmath::Vector<T>::angleBetween(const Vector<U> &other) const
{
  auto mag       = magnitude();
  auto other_mag = other.magnitude();
  if ( mag == 0 || other_mag == 0 ) {
    throw ZeroMagnitudeException("Cannot compute angle between with a magnitude of zero");
  }

  return std::acos(dotProduct(other) / (mag * other_mag));
}

template<Numerical T>
template<Numerical U>
double supmath::Vector<T>::distanceBetween(const Vector<U> &other) const
{
  if ( getSize() != other.getSize() ) {
    throw VectorInconsistentDimensionException(
        "Cannot compute distance between because one vector has size of " + std::to_string(getSize())
        + " and the other has size of " + std::to_string(other.getSize())
    );
  }

  double sum_of_diff = 0;
  for ( int index = 0; index < getSize(); index++ ) {
    auto diff    = _elements[index] - other[index];
    sum_of_diff += diff * diff;
  }

  return std::sqrt(sum_of_diff);
}

template<Numerical T>
template<Numerical U>
Vector<double> supmath::Vector<T>::projectionOnto(const Vector<U> &other) const
{
  auto other_mag = other.magnitude();
  if ( other_mag == 0 ) {
    throw ZeroMagnitudeException("Cannot project onto a vector with zero magnitude");
  }

  return (dotProduct(other) / (other_mag * other_mag)) * other;
}

template<Numerical T>
double supmath::Vector<T>::magnitude() const
{
  double sum_of_square = 0;
  for ( auto &element : _elements ) {
    sum_of_square += element * element;
  }

  return std::sqrt(sum_of_square);
}

template<Numerical T>
Vector<double> supmath::Vector<T>::unitVector() const
{
  auto mag = magnitude();
  if ( mag == 0 ) {
    throw ZeroMagnitudeException("Cannot compute unit vector with a magnitude of zero");
  }

  return *this * (1.0 / magnitude());
}

template<Numerical T>
template<Numerical U>
bool supmath::Vector<T>::operator==(const Vector<U> &other) const
{
  if ( getSize() != other.getSize() ) {
    return false;
  }

  double epsilon = std::is_same<T, float>::value ? 1e-5 : 1e-9;

  for ( size_t index = 0; index < getSize(); index++ ) {
    if ( !numeric_almost_equal(_elements[index], other[index]) ) {
      return false;
    }
  }

  return true;
}

template<Numerical T, Numerical U>
auto supmath::operator*(const T factor, const Vector<U> &vector) -> Vector<std::common_type_t<T, U>>
{
  return _multiplyVectorScalar(factor, vector);
}

template<Numerical T>
std::ostream &supmath::operator<<(std::ostream &os, const Vector<T> &vector)
{
  std::string message = "Vector([";

  for ( size_t index = 0; index < vector.getSize(); index++ ) {
    message += std::to_string(vector[index]);
    if ( index != vector.getSize() - 1 ) {
      message += ", ";
    }
  }

  message += "])";
  os << message;

  return os;
}

// Explicit Instantiation
template class Vector<short>;
template class Vector<int>;
template class Vector<long>;
template class Vector<long long>;
template class Vector<float>;
template class Vector<double>;

template std::ostream &supmath::operator<<(std::ostream &os, const Vector<short> &matrix);
template std::ostream &supmath::operator<<(std::ostream &os, const Vector<int> &matrix);
template std::ostream &supmath::operator<<(std::ostream &os, const Vector<long> &matrix);
template std::ostream &supmath::operator<<(std::ostream &os, const Vector<long long> &matrix);
template std::ostream &supmath::operator<<(std::ostream &os, const Vector<float> &matrix);
template std::ostream &supmath::operator<<(std::ostream &os, const Vector<double> &matrix);

INSTANTIATION_VECTOR_CONSTRUCTOR(short, int)
INSTANTIATION_VECTOR_CONSTRUCTOR(short, long)
INSTANTIATION_VECTOR_CONSTRUCTOR(short, long long)
INSTANTIATION_VECTOR_CONSTRUCTOR(short, float)
INSTANTIATION_VECTOR_CONSTRUCTOR(short, double)
INSTANTIATION_VECTOR_CONSTRUCTOR(int, short)
INSTANTIATION_VECTOR_CONSTRUCTOR(int, long)
INSTANTIATION_VECTOR_CONSTRUCTOR(int, long long)
INSTANTIATION_VECTOR_CONSTRUCTOR(int, float)
INSTANTIATION_VECTOR_CONSTRUCTOR(int, double)
INSTANTIATION_VECTOR_CONSTRUCTOR(long, short)
INSTANTIATION_VECTOR_CONSTRUCTOR(long, int)
INSTANTIATION_VECTOR_CONSTRUCTOR(long, long long)
INSTANTIATION_VECTOR_CONSTRUCTOR(long, float)
INSTANTIATION_VECTOR_CONSTRUCTOR(long, double)
INSTANTIATION_VECTOR_CONSTRUCTOR(long long, short)
INSTANTIATION_VECTOR_CONSTRUCTOR(long long, int)
INSTANTIATION_VECTOR_CONSTRUCTOR(long long, long)
INSTANTIATION_VECTOR_CONSTRUCTOR(long long, float)
INSTANTIATION_VECTOR_CONSTRUCTOR(long long, double)
INSTANTIATION_VECTOR_CONSTRUCTOR(float, short)
INSTANTIATION_VECTOR_CONSTRUCTOR(float, int)
INSTANTIATION_VECTOR_CONSTRUCTOR(float, long)
INSTANTIATION_VECTOR_CONSTRUCTOR(float, long long)
INSTANTIATION_VECTOR_CONSTRUCTOR(float, double)
INSTANTIATION_VECTOR_CONSTRUCTOR(double, short)
INSTANTIATION_VECTOR_CONSTRUCTOR(double, int)
INSTANTIATION_VECTOR_CONSTRUCTOR(double, long)
INSTANTIATION_VECTOR_CONSTRUCTOR(double, long long)
INSTANTIATION_VECTOR_CONSTRUCTOR(double, float)

INSTANTIATION_VECTOR_FUNCTIONS(short, short)
INSTANTIATION_VECTOR_FUNCTIONS(short, int)
INSTANTIATION_VECTOR_FUNCTIONS(short, long)
INSTANTIATION_VECTOR_FUNCTIONS(short, long long)
INSTANTIATION_VECTOR_FUNCTIONS(short, float)
INSTANTIATION_VECTOR_FUNCTIONS(short, double)
INSTANTIATION_VECTOR_FUNCTIONS(int, short)
INSTANTIATION_VECTOR_FUNCTIONS(int, int)
INSTANTIATION_VECTOR_FUNCTIONS(int, long)
INSTANTIATION_VECTOR_FUNCTIONS(int, long long)
INSTANTIATION_VECTOR_FUNCTIONS(int, float)
INSTANTIATION_VECTOR_FUNCTIONS(int, double)
INSTANTIATION_VECTOR_FUNCTIONS(long, short)
INSTANTIATION_VECTOR_FUNCTIONS(long, int)
INSTANTIATION_VECTOR_FUNCTIONS(long, long)
INSTANTIATION_VECTOR_FUNCTIONS(long, long long)
INSTANTIATION_VECTOR_FUNCTIONS(long, float)
INSTANTIATION_VECTOR_FUNCTIONS(long, double)
INSTANTIATION_VECTOR_FUNCTIONS(long long, short)
INSTANTIATION_VECTOR_FUNCTIONS(long long, int)
INSTANTIATION_VECTOR_FUNCTIONS(long long, long)
INSTANTIATION_VECTOR_FUNCTIONS(long long, long long)
INSTANTIATION_VECTOR_FUNCTIONS(long long, float)
INSTANTIATION_VECTOR_FUNCTIONS(long long, double)
INSTANTIATION_VECTOR_FUNCTIONS(float, short)
INSTANTIATION_VECTOR_FUNCTIONS(float, int)
INSTANTIATION_VECTOR_FUNCTIONS(float, long)
INSTANTIATION_VECTOR_FUNCTIONS(float, long long)
INSTANTIATION_VECTOR_FUNCTIONS(float, float)
INSTANTIATION_VECTOR_FUNCTIONS(float, double)
INSTANTIATION_VECTOR_FUNCTIONS(double, short)
INSTANTIATION_VECTOR_FUNCTIONS(double, int)
INSTANTIATION_VECTOR_FUNCTIONS(double, long)
INSTANTIATION_VECTOR_FUNCTIONS(double, long long)
INSTANTIATION_VECTOR_FUNCTIONS(double, float)
INSTANTIATION_VECTOR_FUNCTIONS(double, double)

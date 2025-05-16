#ifndef VECTOR_FUNCS_H
#define VECTOR_FUNCS_H

#include "vector.h"

namespace supmath
{
namespace vector_funcs
{
template<Numerical T, Numerical U>
auto multiplyVectorScalar(const T factor, const Vector<U> &vector) -> Vector<std::common_type_t<T, U>>
{
  using CommonType = std::common_type_t<T, U>;
  Vector<CommonType> result(vector.getSize());

  for ( size_t index = 0; index < vector.getSize(); index++ ) {
    result[index] = vector[index] * factor;
  }

  return result;
}
}    // namespace vector_funcs
}    // namespace supmath

#endif    // VECTOR_FUNCS_H

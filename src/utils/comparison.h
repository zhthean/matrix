#ifndef COMPARISON_H
#define COMPARISON_H

#include <cmath>

namespace supmath
{
namespace comparison
{
template<typename T, typename U>
bool numeric_almost_equal(T left_num, U right_num, double epsilon = 1e-9)
{
  if constexpr ( std::is_integral_v<T> && std::is_integral_v<U> ) {
    return left_num == right_num;
  }
  else {
    return std::abs(left_num - right_num) <= epsilon;
  }
}
}    // namespace comparison
}    // namespace supmath

#endif    // COMPARISON_H

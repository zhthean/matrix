#ifndef TYPES_H
#define TYPES_H

#include <type_traits>

namespace supmath {
template <typename T>
concept Numerical =
    std::is_integral_v<T> && !std::is_same_v<T, char> && !std::is_same_v<T, bool> || std::is_floating_point_v<T>;
}

#endif // TYPES_H

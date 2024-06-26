#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>
namespace poly::utils {

template <size_t Size>
using signed_integer_t = std::conditional_t<
  Size == 8, int64_t,
  std::conditional_t<
    Size == 4, int32_t,
    std::conditional_t<Size == 2, int16_t,
                       std::conditional_t<Size == 1, int8_t, __int128_t>>>>;
template <size_t Size>
using unsigned_integer_t = std::make_unsigned_t<signed_integer_t<Size>>;

constexpr auto widen(std::signed_integral auto x) {
  return signed_integer_t<2 * sizeof(decltype(x))>(x);
}
constexpr auto widen(std::unsigned_integral auto x) {
  return unsigned_integer_t<2 * sizeof(decltype(x))>(x);
}
} // namespace poly::utils

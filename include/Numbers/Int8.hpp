#pragma once
#include "Utilities/Invariant.hpp"
#include <compare>
#include <concepts>
#include <cstdint>
#include <functional>
#include <type_traits>

namespace poly::numbers {

enum class Flag8 : unsigned char {};
// 8 bit unsigned integer, overflow UB
enum class u8 : unsigned char {};
// 8 bit signed integer, overflow UB
enum class i8 : signed char {};

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator==(u8 x, std::integral auto y) -> bool {
  return static_cast<decltype(y)>(x) == y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator==(i8 x, std::integral auto y) -> bool {
  return static_cast<decltype(y)>(x) == y;
}

template <typename OP, typename I, typename J>
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
applyop(OP op, I i, J j) {
  if constexpr (std::is_scoped_enum_v<I>) {
    if constexpr (std::is_scoped_enum_v<J>) {
      // I and J are scoped enums
      using T = std::common_type_t<std::underlying_type_t<I>,
                                   std::underlying_type_t<J>>;
      using S =
        std::conditional_t<sizeof(T) < sizeof(int), int, std::make_signed_t<T>>;
      S x = op(static_cast<S>(i), static_cast<S>(j));
      using R = std::conditional_t<
        std::same_as<T, std::underlying_type_t<I>>, I,
        std::conditional_t<std::same_as<T, std::underlying_type_t<J>>, J, T>>;
      R y = static_cast<R>(x);
      utils::invariant(y == x);
      return y;
    } else {
      // I is a scoped enum
      using T = std::common_type_t<std::underlying_type_t<I>, J>;
      using S =
        std::conditional_t<sizeof(T) < sizeof(int), int, std::make_signed_t<T>>;
      S x = op(static_cast<S>(i), static_cast<S>(j));
      using R =
        std::conditional_t<std::same_as<T, std::underlying_type_t<I>>, I, T>;
      R y = static_cast<R>(x);
      utils::invariant(y == x);
      return y;
    }
  } else {
    static_assert(std::is_scoped_enum_v<J>);
    // J is a scoped enum
    using T = std::common_type_t<I, std::underlying_type_t<J>>;
    using S =
      std::conditional_t<sizeof(T) < sizeof(int), int, std::make_signed_t<T>>;
    S x = op(static_cast<S>(i), static_cast<S>(j));
    using R =
      std::conditional_t<std::same_as<T, std::underlying_type_t<J>>, J, T>;
    R y = static_cast<R>(x);
    utils::invariant(y == x);
    return y;
  }
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator+(u8 x, std::integral auto y) {
  return applyop(std::plus<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator+(std::integral auto x, u8 y) {
  return applyop(std::plus<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator+(u8 x, u8 y) -> u8 {
  return applyop(std::plus<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator+=(u8 &x, u8 y) -> u8 {
  x = x + y;
  return x;
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator-(u8 x, std::integral auto y) {
  return applyop(std::minus<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator-(std::integral auto x, u8 y) {
  return applyop(std::minus<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator-(u8 x, u8 y) -> u8 {
  return applyop(std::minus<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator-=(u8 &x, u8 y) -> u8 {
  x = x - y;
  return x;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator*(u8 x, std::integral auto y) {
  return applyop(std::multiplies<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator*(std::integral auto x, u8 y) {
  return applyop(std::multiplies<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator*(u8 x, u8 y) -> u8 {
  return applyop(std::multiplies<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator*=(u8 &x, u8 y) -> u8 {
  x = x * y;
  return x;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator/(u8 x, std::integral auto y) {
  return applyop(std::divides<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator/(std::integral auto x, u8 y) {
  return applyop(std::divides<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator/(u8 x, u8 y) -> u8 {
  return applyop(std::divides<>{}, x, y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator/=(u8 &x, u8 y) -> u8 {
  x = x / y;
  return x;
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator+(i8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) + y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator+(std::signed_integral auto x, i8 y) {
  return x + static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator+(i8 x, i8 y) -> i8 {
  return static_cast<i8>(static_cast<signed char>(x) +
                         static_cast<signed char>(y));
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator-(i8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) - y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator-(std::signed_integral auto x, i8 y) {
  return x - static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator-(i8 x, i8 y) -> i8 {
  return static_cast<i8>(static_cast<signed char>(x) -
                         static_cast<signed char>(y));
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator*(i8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) * y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator*(std::signed_integral auto x, i8 y) {
  return x * static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator*(i8 x, i8 y) -> i8 {
  return static_cast<i8>(static_cast<signed char>(x) *
                         static_cast<signed char>(y));
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator/(i8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) / y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator/(std::signed_integral auto x, i8 y) {
  return x / static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator/(i8 x, i8 y) -> i8 {
  return static_cast<i8>(static_cast<signed char>(x) /
                         static_cast<signed char>(y));
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator&(Flag8 x, std::integral auto y) -> Flag8 {
  return static_cast<Flag8>(static_cast<decltype(y)>(x) & y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator&(std::integral auto x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(x & static_cast<decltype(x)>(y));
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator&(Flag8 x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(static_cast<unsigned char>(x) &
                            static_cast<unsigned char>(y));
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator|(Flag8 x, std::integral auto y) -> Flag8 {
  return static_cast<Flag8>(static_cast<decltype(y)>(x) | y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator|(std::integral auto x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(x | static_cast<decltype(x)>(y));
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator|(Flag8 x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(static_cast<unsigned char>(x) |
                            static_cast<unsigned char>(y));
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator^(Flag8 x, std::integral auto y) -> Flag8 {
  return static_cast<Flag8>(static_cast<decltype(y)>(x) ^ y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator^(std::integral auto x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(x ^ static_cast<decltype(x)>(y));
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator^(Flag8 x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(static_cast<unsigned char>(x) ^
                            static_cast<unsigned char>(y));
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(u8 x, std::integral auto y) -> std::strong_ordering {
  return static_cast<decltype(y)>(x) <=> y;
}
// constexpr auto operator>(UInt8 x, std::integral auto y) -> bool {
//   return static_cast<decltype(y)>(x) > y;
// }
// constexpr auto operator<(UInt8 x, std::integral auto y) -> bool {
//   return static_cast<decltype(y)>(x) < y;
// }
// constexpr auto operator>=(UInt8 x, std::integral auto y) -> bool {
//   return static_cast<decltype(y)>(x) >= y;
// }
// constexpr auto operator<=(UInt8 x, std::integral auto y) -> bool {
//   return static_cast<decltype(y)>(x) <= y;
// }

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(std::integral auto x, u8 y) -> std::strong_ordering {
  return x <=> static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(u8 x, u8 y) -> std::strong_ordering {
  return static_cast<unsigned char>(x) <=> static_cast<unsigned char>(y);
}
static_assert((static_cast<u8>(255) - uint8_t(253)) ==
              (uint8_t(255) - uint8_t(253)));

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(i8 x, std::signed_integral auto y) -> std::strong_ordering {
  return static_cast<decltype(y)>(x) <=> y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(std::signed_integral auto x, i8 y) -> std::strong_ordering {
  return x <=> static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(i8 x, i8 y) -> std::strong_ordering {
  return static_cast<signed char>(x) <=> static_cast<signed char>(y);
}

[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(Flag8 x, std::integral auto y) -> std::strong_ordering {
  return static_cast<decltype(y)>(x) <=> y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator==(Flag8 x, std::integral auto y) -> bool {
  return static_cast<decltype(y)>(x) == y;
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator==(std::integral auto x, Flag8 y) -> bool {
  return x == static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(std::integral auto x, i8 y) -> std::strong_ordering {
  return x <=> static_cast<decltype(x)>(y);
}
[[gnu::always_inline, gnu::artificial]] inline constexpr auto
operator<=>(Flag8 x, Flag8 y) -> std::strong_ordering {
  return static_cast<unsigned char>(x) <=> static_cast<unsigned char>(y);
}

} // namespace poly::numbers

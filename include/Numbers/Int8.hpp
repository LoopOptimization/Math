#pragma once
#include <compare>
#include <concepts>
#include <cstdint>

namespace poly::numbers {
enum class Flag8 : unsigned char {};
enum class UInt8 : unsigned char {};
enum class Int8 : signed char {};

constexpr auto operator+(UInt8 x, std::integral auto y) {
  return static_cast<decltype(y)>(x) + y;
}
constexpr auto operator+(std::integral auto x, UInt8 y) {
  return x + static_cast<decltype(x)>(y);
}
constexpr auto operator+(UInt8 x, UInt8 y) -> UInt8 {
  return static_cast<UInt8>(static_cast<unsigned char>(x) +
                            static_cast<unsigned char>(y));
}
constexpr auto operator+=(UInt8 &x, UInt8 y) -> UInt8 {
  x = x + y;
  return x;
}

constexpr auto operator-(UInt8 x, std::integral auto y) {
  return static_cast<decltype(y)>(x) - y;
}
static_assert((static_cast<UInt8>(255) - uint8_t(253)) ==
              (uint8_t(255) - uint8_t(253)));
constexpr auto operator-(std::integral auto x, UInt8 y) {
  return x - static_cast<decltype(x)>(y);
}
constexpr auto operator-(UInt8 x, UInt8 y) -> UInt8 {
  return static_cast<UInt8>(static_cast<unsigned char>(x) -
                            static_cast<unsigned char>(y));
}
constexpr auto operator*(UInt8 x, std::integral auto y) {
  return static_cast<decltype(y)>(x) * y;
}
constexpr auto operator*(std::integral auto x, UInt8 y) {
  return x * static_cast<decltype(x)>(y);
}
constexpr auto operator*(UInt8 x, UInt8 y) -> UInt8 {
  return static_cast<UInt8>(static_cast<unsigned char>(x) *
                            static_cast<unsigned char>(y));
}
constexpr auto operator/(UInt8 x, std::integral auto y) {
  return static_cast<decltype(y)>(x) / y;
}
constexpr auto operator/(std::integral auto x, UInt8 y) {
  return x / static_cast<decltype(x)>(y);
}
constexpr auto operator/(UInt8 x, UInt8 y) -> UInt8 {
  return static_cast<UInt8>(static_cast<unsigned char>(x) /
                            static_cast<unsigned char>(y));
}

constexpr auto operator+(Int8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) + y;
}
constexpr auto operator+(std::signed_integral auto x, Int8 y) {
  return x + static_cast<decltype(x)>(y);
}
constexpr auto operator+(Int8 x, Int8 y) -> Int8 {
  return static_cast<Int8>(static_cast<signed char>(x) +
                           static_cast<signed char>(y));
}

constexpr auto operator-(Int8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) - y;
}
constexpr auto operator-(std::signed_integral auto x, Int8 y) {
  return x - static_cast<decltype(x)>(y);
}
constexpr auto operator-(Int8 x, Int8 y) -> Int8 {
  return static_cast<Int8>(static_cast<signed char>(x) -
                           static_cast<signed char>(y));
}

constexpr auto operator*(Int8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) * y;
}
constexpr auto operator*(std::signed_integral auto x, Int8 y) {
  return x * static_cast<decltype(x)>(y);
}
constexpr auto operator*(Int8 x, Int8 y) -> Int8 {
  return static_cast<Int8>(static_cast<signed char>(x) *
                           static_cast<signed char>(y));
}

constexpr auto operator/(Int8 x, std::signed_integral auto y) {
  return static_cast<decltype(y)>(x) / y;
}
constexpr auto operator/(std::signed_integral auto x, Int8 y) {
  return x / static_cast<decltype(x)>(y);
}
constexpr auto operator/(Int8 x, Int8 y) -> Int8 {
  return static_cast<Int8>(static_cast<signed char>(x) /
                           static_cast<signed char>(y));
}

constexpr auto operator&(Flag8 x, std::integral auto y) -> Flag8 {
  return static_cast<Flag8>(static_cast<decltype(y)>(x) & y);
}
constexpr auto operator&(std::integral auto x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(x & static_cast<decltype(x)>(y));
}
constexpr auto operator&(Flag8 x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(static_cast<unsigned char>(x) &
                            static_cast<unsigned char>(y));
}
constexpr auto operator|(Flag8 x, std::integral auto y) -> Flag8 {
  return static_cast<Flag8>(static_cast<decltype(y)>(x) | y);
}
constexpr auto operator|(std::integral auto x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(x | static_cast<decltype(x)>(y));
}
constexpr auto operator|(Flag8 x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(static_cast<unsigned char>(x) |
                            static_cast<unsigned char>(y));
}
constexpr auto operator^(Flag8 x, std::integral auto y) -> Flag8 {
  return static_cast<Flag8>(static_cast<decltype(y)>(x) ^ y);
}
constexpr auto operator^(std::integral auto x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(x ^ static_cast<decltype(x)>(y));
}
constexpr auto operator^(Flag8 x, Flag8 y) -> Flag8 {
  return static_cast<Flag8>(static_cast<unsigned char>(x) ^
                            static_cast<unsigned char>(y));
}

constexpr auto operator<=>(UInt8 x,
                           std::integral auto y) -> std::strong_ordering {
  return static_cast<decltype(y)>(x) <=> y;
}
constexpr auto operator==(UInt8 x, std::integral auto y) -> bool {
  return static_cast<decltype(y)>(x) == y;
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

constexpr auto operator<=>(std::integral auto x,
                           UInt8 y) -> std::strong_ordering {
  return x <=> static_cast<decltype(x)>(y);
}
constexpr auto operator<=>(UInt8 x, UInt8 y) -> std::strong_ordering {
  return static_cast<unsigned char>(x) <=> static_cast<unsigned char>(y);
}

constexpr auto
operator<=>(Int8 x, std::signed_integral auto y) -> std::strong_ordering {
  return static_cast<decltype(y)>(x) <=> y;
}
constexpr auto operator==(Int8 x, std::integral auto y) -> bool {
  return static_cast<decltype(y)>(x) == y;
}
constexpr auto operator<=>(std::signed_integral auto x,
                           Int8 y) -> std::strong_ordering {
  return x <=> static_cast<decltype(x)>(y);
}
constexpr auto operator<=>(Int8 x, Int8 y) -> std::strong_ordering {
  return static_cast<signed char>(x) <=> static_cast<signed char>(y);
}

constexpr auto operator<=>(Flag8 x,
                           std::integral auto y) -> std::strong_ordering {
  return static_cast<decltype(y)>(x) <=> y;
}
constexpr auto operator==(Flag8 x, std::integral auto y) -> bool {
  return static_cast<decltype(y)>(x) == y;
}
constexpr auto operator==(std::integral auto x, Flag8 y) -> bool {
  return x == static_cast<decltype(x)>(y);
}
constexpr auto operator<=>(std::integral auto x,
                           Int8 y) -> std::strong_ordering {
  return x <=> static_cast<decltype(x)>(y);
}
constexpr auto operator<=>(Flag8 x, Flag8 y) -> std::strong_ordering {
  return static_cast<unsigned char>(x) <=> static_cast<unsigned char>(y);
}

} // namespace poly::numbers

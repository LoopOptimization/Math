#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#ifndef USE_MODULE
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#else

export module Widen;
import std;
#endif

#ifdef USE_MODULE
export namespace utils {
#else
namespace utils {
#endif

template <std::size_t Size>
using signed_integer_t = std::conditional_t<
  Size == 8, std::int64_t,
  std::conditional_t<
    Size == 4, std::int32_t,
    std::conditional_t<
      Size == 2, std::int16_t,
      std::conditional_t<Size == 1, std::int8_t, __int128_t>>>>;
template <std::size_t Size>
using unsigned_integer_t = std::make_unsigned_t<signed_integer_t<Size>>;

constexpr auto widen(std::signed_integral auto x) {
  return signed_integer_t<2 * sizeof(decltype(x))>(x);
}
constexpr auto widen(std::unsigned_integral auto x) {
  return unsigned_integer_t<2 * sizeof(decltype(x))>(x);
}
} // namespace utils

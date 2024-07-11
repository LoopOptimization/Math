#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#include <concepts>
#include <limits>

#ifndef USE_MODULE
#include "Utilities/Widen.cxx"
#else
export module Saturated;

import Widen;
#endif

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif
template <std::unsigned_integral T> constexpr auto add_sat(T x, T y) -> T {
  T res = x + y;
  return res | -(res < x);
}
template <std::unsigned_integral T> constexpr auto sub_sat(T x, T y) -> T {
  T res = x - y;
  return res & -(res <= x);
}

template <std::unsigned_integral T> constexpr auto mul_sat(T x, T y) -> T {
  auto prod = utils::widen(x) * utils::widen(y);
  T hi = prod >> (8 * sizeof(T));
  T lo = prod;
  return lo | -!!hi;
}

template <std::signed_integral T> constexpr auto add_sat(T x, T y) -> T {
  T res;
  if (!__builtin_add_overflow(x, y, &res)) return res;
  return x > 0 ? std::numeric_limits<T>::max() : std::numeric_limits<T>::min();
}
template <std::signed_integral T> constexpr auto sub_sat(T x, T y) -> T {
  T res;
  if (!__builtin_sub_overflow(x, y, &res)) return res;
  return x > 0 ? std::numeric_limits<T>::max() : std::numeric_limits<T>::min();
}
template <std::signed_integral T> constexpr auto mul_sat(T x, T y) -> T {
  T res;
  if (!__builtin_mul_overflow(x, y, &res)) return res;
  return ((x > 0) == (y > 0)) ? std::numeric_limits<T>::max()
                              : std::numeric_limits<T>::min();
}
} // namespace math
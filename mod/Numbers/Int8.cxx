#ifdef USE_MODULE
module;
#else
#pragma once
#endif
#include "Macros.hxx"

#ifndef USE_MODULE
#include <compare>
#include <concepts>
#include <limits>
#include <ostream>
#include <type_traits>

#include "Utilities/Invariant.cxx"
#else
export module Int8;

import Invariant;
import STL;
#endif

#ifdef USE_MODULE
export namespace numbers {
#else
namespace numbers {
#endif

template <std::integral I, int alias, bool nowrap = false> struct IntWrapper {
  enum class strong : I {};

private:
  static constexpr bool issigned = std::is_signed_v<I>;
  using T = std::conditional_t<nowrap, std::make_signed_t<I>, I>;
  TRIVIAL inline static constexpr auto create(std::integral auto x) -> strong {
    if constexpr (nowrap) {
      utils::assume(x >= std::numeric_limits<I>::min());
      utils::assume(x <= std::numeric_limits<I>::max());
    }
    return x;
  }

  TRIVIAL friend inline constexpr auto operator++(strong &x) -> strong & {
    x = static_cast<strong>(static_cast<T>(x) + T{1});
    return x;
  }
  TRIVIAL friend inline constexpr auto operator++(strong &&x)
    -> decltype(auto) {
    x = static_cast<strong>(static_cast<T>(x) + T{1});
    return x;
  }
  TRIVIAL friend inline constexpr auto operator--(strong &x) -> strong & {
    x = static_cast<strong>(static_cast<T>(x) - T{1});
    return x;
  }
  TRIVIAL friend inline constexpr auto operator--(strong &&x)
    -> decltype(auto) {
    x = static_cast<strong>(static_cast<T>(x) - T{1});
    return x;
  }
  TRIVIAL friend inline constexpr auto operator++(strong &x, int) -> strong {
    strong y = x;
    x = static_cast<strong>(static_cast<T>(x) + T{1});
    return y;
  }
  TRIVIAL friend inline constexpr auto operator--(strong &x, int) -> strong {
    strong y = x;
    x = static_cast<strong>(static_cast<T>(x) - T{1});
    return y;
  }

  TRIVIAL friend inline constexpr auto operator++(strong &&x, int) -> strong {
    strong y = x;
    x = static_cast<strong>(static_cast<T>(x) + T{1});
    return y;
  }
  TRIVIAL friend inline constexpr auto operator--(strong &&x, int) -> strong {
    strong y = x;
    x = static_cast<strong>(static_cast<T>(x) - T{1});
    return y;
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator==(strong x, J y) -> bool
  requires((issigned == std::is_signed_v<J>) ||
           (!issigned && (sizeof(J) > sizeof(I))))
  {
    if constexpr (sizeof(J) >= sizeof(I)) return static_cast<J>(x) == y;
    else return static_cast<I>(x) == static_cast<I>(y);
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator==(J y, strong x) -> bool
  requires((issigned == std::is_signed_v<J>) ||
           (!issigned && (sizeof(J) > sizeof(I))))
  {
    if constexpr (sizeof(J) >= sizeof(I)) return y == static_cast<J>(x);
    else return static_cast<I>(x) == static_cast<I>(y);
  }
  TRIVIAL friend inline constexpr auto operator==(strong x, strong y) -> bool {
    return static_cast<I>(x) == static_cast<I>(y);
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator<=>(strong x, J y)
    -> std::strong_ordering
  requires((issigned == std::is_signed_v<J>) ||
           (!issigned && (sizeof(J) > sizeof(I))))
  {
    if constexpr (sizeof(J) >= sizeof(I)) return static_cast<J>(x) <=> y;
    else return static_cast<I>(x) <=> static_cast<I>(y);
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator<=>(J y, strong x)
    -> std::strong_ordering
  requires((issigned == std::is_signed_v<J>) ||
           (!issigned && (sizeof(J) > sizeof(I))))
  {
    if constexpr (sizeof(J) >= sizeof(I)) return y <=> static_cast<J>(x);
    else return static_cast<I>(x) <=> static_cast<I>(y);
  }
  TRIVIAL friend inline constexpr auto operator<=>(strong x, strong y)
    -> std::strong_ordering {
    return static_cast<I>(x) <=> static_cast<I>(y);
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator+(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(static_cast<J>(x) + y);
    else return static_cast<strong>(static_cast<T>(x) + static_cast<T>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator+(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(y + static_cast<J>(x));
    else return static_cast<strong>(static_cast<T>(x) + static_cast<T>(y));
  }
  TRIVIAL friend inline constexpr auto operator+(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<T>(x) + static_cast<T>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator-(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(static_cast<J>(x) - y);
    else return static_cast<strong>(static_cast<T>(x) - static_cast<T>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator-(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(y - static_cast<J>(x));
    else return static_cast<strong>(static_cast<T>(y) - static_cast<T>(x));
  }
  TRIVIAL friend inline constexpr auto operator-(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<T>(x) - static_cast<T>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator*(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(static_cast<J>(x) * y);
    else return static_cast<strong>(static_cast<T>(x) * static_cast<T>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator*(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(y * static_cast<J>(x));
    else return static_cast<strong>(static_cast<T>(y) * static_cast<T>(x));
  }
  TRIVIAL friend inline constexpr auto operator*(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<T>(x) * static_cast<T>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator/(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(static_cast<J>(x) / y);
    else return static_cast<strong>(static_cast<I>(x) / static_cast<I>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator/(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(y / static_cast<J>(x));
    else return static_cast<strong>(static_cast<I>(y) / static_cast<I>(x));
  }
  TRIVIAL friend inline constexpr auto operator/(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<I>(x) / static_cast<I>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator|(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(static_cast<J>(x) | y);
    else return static_cast<strong>(static_cast<T>(x) | static_cast<T>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator|(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(y | static_cast<J>(x));
    else return static_cast<strong>(static_cast<T>(y) | static_cast<T>(x));
  }
  TRIVIAL friend inline constexpr auto operator|(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<T>(x) | static_cast<T>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator&(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return static_cast<J>(x) & y;
    else return static_cast<strong>(static_cast<T>(x) & static_cast<T>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator&(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return y & static_cast<J>(x);
    else return static_cast<strong>(static_cast<T>(y) & static_cast<T>(x));
  }
  TRIVIAL friend inline constexpr auto operator&(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<T>(x) & static_cast<T>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator^(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(static_cast<J>(x) ^ y);
    else return static_cast<strong>(static_cast<T>(x) ^ static_cast<T>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator^(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return create(y ^ static_cast<J>(x));
    else return static_cast<strong>(static_cast<T>(y) ^ static_cast<T>(x));
  }
  TRIVIAL friend inline constexpr auto operator^(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<T>(x) ^ static_cast<T>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator%(strong x, J y) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return static_cast<J>(x) % y;
    else return static_cast<strong>(static_cast<I>(x) % static_cast<I>(y));
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator%(J y, strong x) -> strong
  requires(issigned == std::is_signed_v<J>)
  {
    if constexpr (sizeof(J) > sizeof(I)) return y % static_cast<J>(x);
    else return static_cast<strong>(static_cast<I>(y) % static_cast<I>(x));
  }
  TRIVIAL friend inline constexpr auto operator%(strong x, strong y) -> strong {
    return static_cast<strong>(static_cast<I>(x) % static_cast<I>(y));
  }

  TRIVIAL friend inline constexpr auto operator<<(strong x, I y) -> strong {
    return static_cast<strong>(static_cast<I>(x) << y);
  }
  TRIVIAL friend inline constexpr auto operator<<(I y, strong x) -> strong {
    return static_cast<strong>(y << static_cast<I>(x));
  }
  TRIVIAL friend inline constexpr auto operator<<(strong x, strong y)
    -> strong {
    return static_cast<strong>(static_cast<I>(x) << static_cast<I>(y));
  }

  TRIVIAL friend inline constexpr auto operator>>(strong x, I y) -> strong {
    return static_cast<strong>(static_cast<I>(x) >> y);
  }
  TRIVIAL friend inline constexpr auto operator>>(I y, strong x) -> strong {
    return static_cast<strong>(y >> static_cast<I>(x));
  }
  TRIVIAL friend inline constexpr auto operator>>(strong x, strong y)
    -> strong {
    return static_cast<strong>(static_cast<I>(x) >> static_cast<I>(y));
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator+=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x + y;
  }

  TRIVIAL friend inline constexpr auto operator+=(strong &x, strong y)
    -> strong & {
    return x = x + y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator-=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x - y;
  }

  TRIVIAL friend inline constexpr auto operator-=(strong &x, strong y)
    -> strong & {
    return x = x - y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator*=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x * y;
  }

  TRIVIAL friend inline constexpr auto operator*=(strong &x, strong y)
    -> strong & {
    return x = x * y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator/=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x / y;
  }

  TRIVIAL friend inline constexpr auto operator/=(strong &x, strong y)
    -> strong & {
    return x = x / y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator%=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x % y;
  }

  TRIVIAL friend inline constexpr auto operator%=(strong &x, strong y)
    -> strong & {
    return x = x % y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator&=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x & y;
  }

  TRIVIAL friend inline constexpr auto operator&=(strong &x, strong y)
    -> strong & {
    return x = x & y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator|=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x | y;
  }

  TRIVIAL friend inline constexpr auto operator|=(strong &x, strong y)
    -> strong & {
    return x = x | y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator^=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x ^ y;
  }

  TRIVIAL friend inline constexpr auto operator^=(strong &x, strong y)
    -> strong & {
    return x = x ^ y;
  }

  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator<<=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x << y;
  }

  TRIVIAL friend inline constexpr auto operator<<=(strong &x, strong y)
    -> strong & {
    return x = x << y;
  }
  template <std::integral J>
  TRIVIAL friend inline constexpr auto operator>>=(strong &x, J y) -> strong &
  requires(issigned == std::is_signed_v<J>)
  {
    return x = x >> y;
  }

  TRIVIAL friend inline constexpr auto operator>>=(strong &x, strong y)
    -> strong & {
    return x = x >> y;
  }
  TRIVIAL friend inline constexpr auto operator!(strong x) -> bool {
    return !static_cast<bool>(x);
  }
  friend auto operator<<(std::ostream &os, strong x) -> std::ostream & {
    return os << static_cast<I>(x);
  }
};
static_assert(++static_cast<IntWrapper<int, 0>::strong>(3) == 4);

using i8 = IntWrapper<signed char, 0>::strong;
using u8 = IntWrapper<unsigned char, 0>::strong;
using Flag8 = IntWrapper<unsigned char, 1>::strong;

static_assert(!bool(i8{}));
static_assert(!bool(u8{}));
static_assert(bool(i8{1}));
static_assert(bool(u8{1}));

} // namespace numbers

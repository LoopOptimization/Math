#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#ifndef USE_MODULE
#include "Math/GreatestCommonDivisor.cxx"
#include "Utilities/Invariant.cxx"
#include "Utilities/TypeCompression.cxx"
#include "Utilities/Widen.cxx"
#include <concepts>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <ostream>
#include <type_traits>
#else
export module Rational;

import GCD;
import Invariant;
import STL;
import TypeCompression;
import Widen;
#endif

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif

struct Rational {
  [[no_unique_address]] int64_t numerator_{0};
  [[no_unique_address]] int64_t denominator_{1};
  // should be invariant that denominator >= 0
  constexpr Rational() = default;
  constexpr Rational(int64_t coef) : numerator_(coef) {};
  constexpr Rational(int coef) : numerator_(coef) {};
  constexpr Rational(int64_t n, int64_t d)
    : numerator_(d > 0 ? n : -n), denominator_(n ? (d > 0 ? d : -d) : 1) {}
  static constexpr auto createPositiveDenominator(int64_t n, int64_t d)
    -> Rational {
    utils::assume(d > 0);
    if (!n) Rational{0, 1};
    int64_t g = gcd(n, d);
    if (g != 1) {
      n /= g;
      d /= g;
    }
    return Rational{n, d};
  }
  static constexpr auto create(int64_t n, int64_t d) -> Rational {
    if (!n) return Rational{0, 1};
    int64_t sign = 2 * (d > 0) - 1;
    return createPositiveDenominator(n * sign, d * sign);
  }

  [[nodiscard]] constexpr auto safeAdd(Rational y) const
    -> std::optional<Rational> {
    auto [xd, yd] = divgcd(denominator_, y.denominator_);
    int64_t a, b, n, d;
    bool o1 = __builtin_smull_overflow(numerator_, yd, &a);
    bool o2 = __builtin_smull_overflow(y.numerator_, xd, &b);
    bool o3 = __builtin_smull_overflow(denominator_, yd, &d);
    bool o4 = __builtin_saddl_overflow(a, b, &n);
    if ((o1 | o2) | (o3 | o4)) return {};
    if (!n) return Rational{0, 1};
    auto [nn, nd] = divgcd(n, d);
    return Rational{nn, nd};
  }
  constexpr auto operator+(Rational y) const -> Rational {
    return *safeAdd(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  constexpr auto operator+=(Rational y) -> Rational & {
    std::optional<Rational> a = *this + y;
    invariant(a.has_value());
    *this = *a;
    return *this;
  }
  [[nodiscard]] constexpr auto safeSub(Rational y) const
    -> std::optional<Rational> {
    auto [xd, yd] = divgcd(denominator_, y.denominator_);
    int64_t a, b, n, d;
    bool o1 = __builtin_smull_overflow(numerator_, yd, &a);
    bool o2 = __builtin_smull_overflow(y.numerator_, xd, &b);
    bool o3 = __builtin_smull_overflow(denominator_, yd, &d);
    bool o4 = __builtin_ssubl_overflow(a, b, &n);
    if ((o1 | o2) | (o3 | o4)) return {};
    if (!n) return Rational{0, 1};
    auto [nn, nd] = divgcd(n, d);
    return Rational{nn, nd};
  }
  constexpr auto operator-(Rational y) const -> Rational {
    return *safeSub(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  constexpr auto operator-=(Rational y) -> Rational & {
    std::optional<Rational> a = *this - y;
    invariant(a.has_value());
    *this = *a;
    return *this;
  }
  [[nodiscard]] constexpr auto safeMul(int64_t y) const
    -> std::optional<Rational> {
    auto [xd, yn] = divgcd(denominator_, y);
    int64_t n;
    if (__builtin_mul_overflow(numerator_, yn, &n)) return {};
    return Rational{n, xd};
  }
  [[nodiscard]] constexpr auto safeMul(Rational y) const
    -> std::optional<Rational> {
    if ((numerator_ == 0) | (y.numerator_ == 0)) return Rational{0, 1};
    auto [xn, yd] = divgcd(numerator_, y.denominator_);
    auto [xd, yn] = divgcd(denominator_, y.numerator_);
    int64_t n, d;
    bool o1 = __builtin_smull_overflow(xn, yn, &n);
    bool o2 = __builtin_smull_overflow(xd, yd, &d);
    if (o1 | o2) return {};
    return Rational{n, d};
  }
  constexpr auto operator*(Rational y) const -> Rational {
    return *safeMul(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  constexpr auto operator*=(Rational y) -> Rational & {
    if ((numerator_ != 0) & (y.numerator_ != 0)) {
      auto [xn, yd] = divgcd(numerator_, y.denominator_);
      auto [xd, yn] = divgcd(denominator_, y.numerator_);
      numerator_ = xn * yn;
      denominator_ = xd * yd;
    } else {
      numerator_ = 0;
      denominator_ = 1;
    }
    return *this;
  }
  [[nodiscard]] constexpr auto inv() const -> Rational {
    if (numerator_ > 0) return Rational{denominator_, numerator_};
    invariant(denominator_ != std::numeric_limits<int64_t>::min());
    invariant(numerator_ != 0);
    return Rational{-denominator_, -numerator_};
  }
  [[nodiscard]] constexpr auto safeDiv(Rational y) const
    -> std::optional<Rational> {
    return (*this) * y.inv();
  }
  constexpr auto operator/(Rational y) const -> Rational {
    return *safeDiv(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  // *this -= a*b
  constexpr auto fnmadd(Rational a, Rational b) -> bool {
    if (std::optional<Rational> ab = a.safeMul(b)) {
      if (std::optional<Rational> c = safeSub(*ab)) {
        *this = *c;
        return false;
      }
    }
    return true;
  }
  constexpr auto div(Rational a) -> bool {
    if (std::optional<Rational> d = safeDiv(a)) {
      *this = *d;
      return false;
    }
    return true;
  }
  // Rational operator/=(Rational y) { return (*this) *= y.inv(); }

  constexpr auto operator==(Rational y) const -> bool {
    return (numerator_ == y.numerator_) & (denominator_ == y.denominator_);
  }
  constexpr auto operator!=(Rational y) const -> bool {
    return (numerator_ != y.numerator_) | (denominator_ != y.denominator_);
  }
  [[nodiscard]] constexpr auto isEqual(int64_t y) const -> bool {
    utils::invariant(denominator_ > 0);
    if (denominator_ == 1) return (numerator_ == y);
    return false;
  }
  constexpr auto operator==(int y) const -> bool { return isEqual(y); }
  constexpr auto operator==(int64_t y) const -> bool { return isEqual(y); }
  constexpr auto operator!=(int y) const -> bool { return !isEqual(y); }
  constexpr auto operator!=(int64_t y) const -> bool { return !isEqual(y); }
  constexpr auto operator<(Rational y) const -> bool {
    return (utils::widen(numerator_) * utils::widen(y.denominator_)) <
           (utils::widen(y.numerator_) * utils::widen(denominator_));
  }
  constexpr auto operator<=(Rational y) const -> bool {
    return (utils::widen(numerator_) * utils::widen(y.denominator_)) <=
           (utils::widen(y.numerator_) * utils::widen(denominator_));
  }
  constexpr auto operator>(Rational y) const -> bool {
    return (utils::widen(numerator_) * utils::widen(y.denominator_)) >
           (utils::widen(y.numerator_) * utils::widen(denominator_));
  }
  constexpr auto operator>=(Rational y) const -> bool {
    return (utils::widen(numerator_) * utils::widen(y.denominator_)) >=
           (utils::widen(y.numerator_) * utils::widen(denominator_));
  }
  constexpr auto operator>=(int y) const -> bool {
    return *this >= Rational(y);
  }
  [[nodiscard]] constexpr auto isInteger() const -> bool {
    return denominator_ == 1;
  }
  constexpr void negate() { numerator_ = -numerator_; }
  constexpr explicit operator bool() const { return numerator_ != 0; }
  constexpr explicit operator double() const {
    return double(numerator_) / double(denominator_);
  }

#ifndef NDEBUG
  [[gnu::used]] void dump() const { std::cout << *this << "\n"; }
#endif

private:
  friend constexpr auto operator+(Rational x, int64_t y) -> Rational {
    return Rational{x.numerator_ + y * x.denominator_, x.denominator_};
  }
  friend constexpr auto operator+(int64_t y, Rational x) -> Rational {
    return x + y;
  }
  friend constexpr auto operator-(Rational x, int64_t y) -> Rational {
    return Rational{x.numerator_ - y * x.denominator_, x.denominator_};
  }
  friend constexpr auto operator-(int64_t y, Rational x) -> Rational {
    return Rational{y * x.denominator_ - x.numerator_, x.denominator_};
  }
  friend constexpr auto operator*(Rational x, int64_t y) -> Rational {
    return *x.safeMul(y); // NOLINT(bugprone-unchecked-optional-access)
  }
  friend constexpr auto operator*(int64_t y, Rational x) -> Rational {
    return x * y;
  }
  friend constexpr auto operator==(int64_t x, Rational y) -> bool {
    return y == x;
  }
  friend auto operator<<(std::ostream &os, const Rational &x)
    -> std::ostream & {
    os << x.numerator_;
    if (x.denominator_ != 1) os << " // " << x.denominator_;
    return os;
  }
  friend constexpr auto gcd(Rational x, Rational y) -> std::optional<Rational> {
    return Rational{gcd(x.numerator_, y.numerator_),
                    lcm(x.denominator_, y.denominator_)};
  }
};
} // namespace math
#ifdef USE_MODULE
export {
#endif
  template <> struct std::common_type<math::Rational, int> {
    using type = math::Rational;
  };
  template <> struct std::common_type<int, math::Rational> {
    using type = math::Rational;
  };
  template <> struct std::common_type<math::Rational, int64_t> {
    using type = math::Rational;
  };
  template <> struct std::common_type<int64_t, math::Rational> {
    using type = math::Rational;
  };
#ifdef USE_MODULE
} // namespace std
#endif
static_assert(
  std::same_as<std::common_type_t<math::Rational, int>, math::Rational>);
static_assert(
  std::same_as<std::common_type_t<int, math::Rational>, math::Rational>);
static_assert(
  std::same_as<utils::decompressed_t<math::Rational>, math::Rational>);

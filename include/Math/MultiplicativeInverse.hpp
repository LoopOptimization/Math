#pragma once

#include "Numbers/Int8.hpp"
#include "Utilities/Invariant.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>
namespace poly::math {
using numbers::i8, numbers::u8;

constexpr auto widen(std::signed_integral auto x) {
  static constexpr auto sz = sizeof(decltype(x));
  if constexpr (sz == 8) return __int128_t(x);
  else if constexpr (sz == 4) return int64_t(x);
  else if constexpr (sz == 2) return int32_t(x);
  else {
    static_assert(sz == 1);
    return int16_t(x);
  }
}
constexpr auto widen(std::unsigned_integral auto x) {
  static constexpr auto sz = sizeof(decltype(x));
  if constexpr (sz == 8) return __uint128_t(x);
  else if constexpr (sz == 4) return uint64_t(x);
  else if constexpr (sz == 2) return uint32_t(x);
  else {
    static_assert(sz == 1);
    return uint16_t(x);
  }
}

#ifdef __clang__
// clang requires linking compiler-rt for 128-bit integer mul_with_overflow
// https://bugs.llvm.org/show_bug.cgi?id=16404
// Currently, instead of linking it, we're just not sanitizing those
#define CLANGNOUBSAN __attribute__((no_sanitize("undefined")))
#else
#define CLANGNOUBSAN
#endif

CLANGNOUBSAN constexpr auto _mul_high(__uint128_t a,
                                      __uint128_t b) -> __uint128_t {
  static constexpr auto shift = 16 * 4;
  __uint128_t mask = std::numeric_limits<__uint128_t>::max();
  __uint128_t a1 = a >> shift, a2 = a & mask;
  __uint128_t b1 = b >> shift, b2 = b & mask;
  __uint128_t a1b1 = a1 * b1, a1b2 = a1 * b2, a2b1 = a2 * b1, a2b2 = a2 * b2;
  __uint128_t carry =
    ((a1b2 & mask) + (a2b1 & mask) + (a2b2 >> shift)) >> shift;
  return a1b1 + (a1b2 >> shift) + (a2b1 >> shift) + carry;
}
CLANGNOUBSAN constexpr auto _mul_high(__int128_t a,
                                      __int128_t b) -> __int128_t {
  static constexpr auto shift = 16 * 8 - 1;
  auto t1 = std::bit_cast<__uint128_t>((a >> shift) & b),
       t2 = std::bit_cast<__uint128_t>((b >> shift) & a);
  return std::bit_cast<__int128_t>(
    _mul_high(std::bit_cast<__uint128_t>(a), std::bit_cast<__uint128_t>(b)) -
    t1 - t2);
}
template <std::integral T>
CLANGNOUBSAN constexpr auto _mul_high(T a, T b) -> T {
  return T((widen(a) * widen(b)) >> (8 * sizeof(T)));
}
#undef CLANGNOUBSAN

template <std::integral T> class MultiplicativeInverse {
  static constexpr bool issigned = std::is_signed_v<T>;
  using AMT = std::conditional_t<issigned, numbers::i8, bool>;
  T divisor_;
  T multiplier_;
  AMT addmul_;
  numbers::u8 shift_;

  friend constexpr auto operator/(T a, MultiplicativeInverse b) -> T {
    T x = _mul_high(a, b.multiplier_);
    if constexpr (issigned) {
      x += a * T(b.addmul_);
      return std::abs(b.divisor_) == 1 ? (a * b.divisor_)
                                       : (std::signbit(x) + (x >> T(b.shift_)));
    } else {
      x = b.addmul_ ? ((((a - x) >> 1)) + x) : x;
      return b.divisor_ == 1 ? a : x >> T(b.shift_);
    }
  }
  friend constexpr auto operator%(T a, MultiplicativeInverse b) -> T {
    return a - (a / b) * b.divisor_;
  }

public:
  constexpr auto divrem(T a) -> std::array<T, 2> {
    T d = a / (*this);
    return {d, a - d * divisor_};
  }

  constexpr MultiplicativeInverse(T d) {
    utils::invariant(d != 0);
    if constexpr (issigned) {
      using UT = std::make_unsigned_t<T>;
      UT signedmin = std::bit_cast<UT>(std::numeric_limits<T>::min());

      // Algorithm from Hacker's Delight, section 10-4
      UT ad = static_cast<UT>(std::abs(d));
      UT t = signedmin + std::signbit(d);
      UT anc = t - UT{1} - (t % ad); // absolute value of nc
      UT p = sizeof(T) * 8 - 1;
      UT q1 = signedmin / anc, r1 = signedmin % anc;
      UT q2 = signedmin / ad, r2 = signedmin % ad;
      for (;;) {
        ++p; // loop until we find a satisfactory p
        // update q1, r1 = divrem(2^p, abs(nc))
        q1 <<= 1;
        r1 <<= 1;
        if (r1 >= anc) {
          ++q1;
          r1 -= anc;
        }
        // update q2, r2 = divrem(2^p, abs(d))
        q2 <<= 1;
        r2 <<= 1;
        if (r2 >= ad) {
          ++q2;
          r2 -= ad;
        }
        UT delta = ad - r2;
        if (!(q1 < delta || (q1 == delta && r1 == 0))) break;
      }
      T m = std::bit_cast<T>(++q2);
      m = std::signbit(d) ? -m : m;
      UT s = p - sizeof(T) * 8;
      divisor_ = d;
      multiplier_ = m;
      if ((d > 0) && (m < 0)) addmul_ = i8(1);
      else if ((d < 0) && (m > 0)) addmul_ = i8(-1);
      else addmul_ = {};
      shift_ = u8(s);
    } else {
      bool add = false;
      T p = 8 * sizeof(T) - 1;
      T signedmin = T{1} << p;
      T signedmax = signedmin - T{1};
      T allones = std::numeric_limits<T>::max();
      T nc = allones - ((allones - d) % d);
      T q1 = signedmin / nc, r1 = signedmin % nc;
      T q2 = signedmax / d, r2 = signedmax % d;
      for (;;) {
        ++p;
        if (r1 >= nc - r1) {
          q1 = q1 + q1 + T{1};
          r1 = r1 + r1 - nc;
        } else {
          q1 = q1 + q1;
          r1 = r1 + r1;
        }
        if ((r2 + T{1}) >= (d - r2)) {
          add |= q2 >= signedmax;
          q2 = q2 + q2 + T{1};
          r2 = r2 + r2 + T{1} - d;
        } else {
          add |= q2 >= signedmin;
          q2 = q2 + q2;
          r2 = r2 + r2 + T{1};
        }
        T delta = d - T{1} - r2;
        if (!(p < sizeof(T) * 16 && (q1 < delta || (q1 == delta && r1 == 0))))
          break;
      }
      T m = q2 + T{1};               // resulting magic number
      T s = p - sizeof(T) * 8 - add; // resulting shift
      divisor_ = d;
      multiplier_ = m;
      addmul_ = add;
      shift_ = u8(s);
    }
  }
};

template <std::integral T>
MultiplicativeInverse(T d) -> MultiplicativeInverse<T>;

} // namespace poly::math


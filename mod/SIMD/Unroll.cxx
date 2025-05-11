#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#include "LoopMacros.hxx"
#include "Macros.hxx"
#ifndef USE_MODULE
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <type_traits>

#include "Math/AxisTypes.cxx"
#include "SIMD/Indexing.cxx"
#include "SIMD/Intrin.cxx"
#include "SIMD/Masks.cxx"
#include "SIMD/Vec.cxx"
#else
export module SIMD:Unroll;

import :Vec;
import :Intrin;
import :Index;
import :Mask;
import AxisTypes;
import std;
#endif

#ifdef USE_MODULE
export namespace simd {
#else
namespace simd {
#endif
// template <typename T, std::ptrdiff_t W, typename S>
// TRIVIAL constexpr auto vcvt(Vec<W, S> v) {
//   if constexpr (std::same_as<T, S>) return v;
//   else if constexpr (W == 1) return T(v);
//   else return __builtin_convertvector(v, Vec<W, T>);
// }
// Vector goes across cols
template <std::ptrdiff_t R, std::ptrdiff_t C, std::ptrdiff_t N, typename T>
struct Unroll {
  static constexpr std::ptrdiff_t W =
    std::ptrdiff_t(std::bit_ceil(std::size_t(N)));
  using VT = Vec<W, T>;
  static_assert(R * C > 0);
  VT data_[R * C];
  TRIVIAL constexpr auto operator[](std::ptrdiff_t i) -> VT & {
    return data_[i];
  }
  TRIVIAL constexpr auto operator[](std::ptrdiff_t r, std::ptrdiff_t c)
    -> VT & {
    return data_[r * C + c];
  }
  TRIVIAL constexpr auto operator[](std::ptrdiff_t i) const -> VT {
    return data_[i];
  }
  TRIVIAL constexpr auto operator[](std::ptrdiff_t r, std::ptrdiff_t c) const
    -> VT {
    return data_[r * C + c];
  }
  template <typename U>
  TRIVIAL constexpr operator Unroll<R, C, N, U>() const
  requires(!std::same_as<T, U>)
  {
    Unroll<R, C, N, U> x;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i)
      if constexpr (W == 1) x.data_[i] = U(data_[i]);
      else x.data_[i] = __builtin_convertvector(data_[i], Vec<W, U>);
    return x;
  }
  TRIVIAL constexpr auto operator-() {
    Unroll a;
    for (std::ptrdiff_t i = 0; i < R * C; ++i) a.data_[i] = -data_[i];
    return a;
  }
  TRIVIAL constexpr auto operator+=(const Unroll &a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] += a.data_[i];
    return *this;
  }
  TRIVIAL constexpr auto operator-=(const Unroll &a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] -= a.data_[i];
    return *this;
  }
  TRIVIAL constexpr auto operator*=(const Unroll &a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] *= a.data_[i];
    return *this;
  }
  TRIVIAL constexpr auto operator/=(const Unroll &a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] /= a.data_[i];
    return *this;
  }
  TRIVIAL constexpr auto operator+=(VT a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] += a;
    return *this;
  }
  TRIVIAL constexpr auto operator-=(VT a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] -= a;
    return *this;
  }
  TRIVIAL constexpr auto operator*=(VT a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] *= a;
    return *this;
  }
  TRIVIAL constexpr auto operator/=(VT a) -> Unroll & {
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) data_[i] /= a;
    return *this;
  }
  TRIVIAL constexpr auto operator+=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    return (*this) += vbroadcast<W, T>(a);
  }
  TRIVIAL constexpr auto operator-=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    return (*this) -= vbroadcast<W, T>(a);
  }
  TRIVIAL constexpr auto operator*=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    return (*this) *= vbroadcast<W, T>(a);
  }
  TRIVIAL constexpr auto operator/=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    return (*this) /= vbroadcast<W, T>(a);
  }

private:
  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator+(Unroll a, Unroll<R1, C1, W1, T1> b) {
    return applyop(a, b, std::plus<>{});
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator-(Unroll a, Unroll<R1, C1, W1, T1> b) {
    return applyop(a, b, std::minus<>{});
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator*(Unroll a,
                                          const Unroll<R1, C1, W1, T1> &b) {
    return applyop(a, b, std::multiplies<>{});
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator/(Unroll a,
                                          const Unroll<R1, C1, W1, T1> &b) {
    return applyop(a, b, std::divides<>{});
  }

  TRIVIAL friend constexpr auto operator+(Unroll a, VT b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a.data_[i] + b;
    return c;
  }
  TRIVIAL friend constexpr auto operator-(Unroll a, VT b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a.data_[i] - b;
    return c;
  }
  TRIVIAL friend constexpr auto operator*(Unroll a, VT b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a.data_[i] * b;
    return c;
  }
  TRIVIAL friend constexpr auto operator/(Unroll a, VT b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a.data_[i] / b;
    return c;
  }
  TRIVIAL friend constexpr auto operator+(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a + vbroadcast<W, T>(b);
  }
  TRIVIAL friend constexpr auto operator-(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a - vbroadcast<W, T>(b);
  }
  TRIVIAL friend constexpr auto operator*(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a * vbroadcast<W, T>(b);
  }
  TRIVIAL friend constexpr auto operator/(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a / vbroadcast<W, T>(b);
  }

  TRIVIAL friend constexpr auto operator+(VT a, Unroll b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a + b.data_[i];
    return c;
  }
  TRIVIAL friend constexpr auto operator-(VT a, Unroll b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a - b.data_[i];
    return c;
  }
  TRIVIAL friend constexpr auto operator*(VT a, Unroll b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a * b.data_[i];
    return c;
  }
  TRIVIAL friend constexpr auto operator/(VT a, Unroll b) -> Unroll {
    Unroll<R, C, W, T> c;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t i = 0; i < R * C; ++i) c.data_[i] = a / b.data_[i];
    return c;
  }
  TRIVIAL friend constexpr auto operator+(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) + a;
  }
  TRIVIAL friend constexpr auto operator-(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) - a;
  }
  TRIVIAL friend constexpr auto operator*(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) * a;
  }
  TRIVIAL friend constexpr auto operator/(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) / a;
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1, typename Op>
  TRIVIAL friend constexpr auto applyop(Unroll a, Unroll<R1, C1, W1, T1> b,
                                        Op op) {
    // Possibilities:
    // 1. All match
    // 2. We had separate unrolls across rows and columns, and some arrays
    // were indexed by one or two of them.
    // In the latter case, we could have arrays indexed by rows, cols, or both.
    if constexpr (!std::same_as<T, T1>) {
      using PT = std::common_type_t<T, T1>;
      return applyop(Unroll<R, C, W, PT>(a), Unroll<R1, C1, W1, PT>(b), op);
    } else if constexpr (W == W1) {
      // both were indexed by cols, and `C2`s should also match
      // or neither were, and they should still match.
      static_assert(C == C1);
      if constexpr (R == R1) {
        // Both have the same index across rows
        Unroll<R, C, W, T> c;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t i = 0; i < R * C; ++i)
          c.data_[i] = op(a.data_[i], b.data_[i]);
        return c;
      } else if constexpr (R == 1) { // R1 > 0
        // `a` was indexed across cols only
        Unroll<R1, C, W, T> z;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t r = 0; r < R1; ++r) {
          POLYMATHFULLUNROLL
          for (std::ptrdiff_t c = 0; c < C; ++c)
            z[r, c] = op(a.data_[c], b[r, c]);
        }
        return z;
      } else {
        static_assert(R1 == 1); // R > 0
        // `b` was indexed across cols only
        Unroll<R, C, W, T> z;
        if constexpr (C == 1) {
          POLYMATHFULLUNROLL
          for (std::ptrdiff_t r = 0; r < R; ++r)
            z.data_[r] = op(a.data_[r], b.vec_);
        } else {
          POLYMATHFULLUNROLL
          for (std::ptrdiff_t r = 0; r < R; ++r) {
            POLYMATHFULLUNROLL
            for (std::ptrdiff_t c = 0; c < C; ++c)
              z[r, c] = op(a[r, c], b.data_[c]);
          }
        }
        return z;
      }
    } else if constexpr (W == 1) {
      static_assert(R == 1 || C == 1);
      static constexpr std::ptrdiff_t R2 = R == 1 ? C : R;
      // `a` was indexed by row only
      Unroll<R2, C1, W1, T> z;
      static_assert(R1 == R2 || R1 == 1);
      static_assert(R2 != 1);
      if constexpr (C1 == 1) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t r = 0; r < R2; ++r)
          if constexpr (R2 == R1) z.data_[r] = op(a.data_[r], b.vec_[r]);
          else z.data_[r] = op(a.data_[r], b.vec_);
      } else {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t r = 0; r < R2; ++r) {
          POLYMATHFULLUNROLL
          for (std::ptrdiff_t c = 0; c < C1; ++c)
            if constexpr (R2 == R1) z[r, c] = op(a.data_[r], b[r, c]);
            else z[r, c] = op(a.data_[r], b.data_[c]);
        }
      }
      return z;
    } else {
      static_assert(W1 == 1);
      static_assert(R1 == 1 || C1 == 1);
      constexpr std::ptrdiff_t R2 = R1 == 1 ? C1 : R1;
      // `b` was indexed by row only
      Unroll<R2, C, W, T> z;
      static_assert(R == R2 || R == 1);
      if constexpr (R2 == 1) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          z.data_[c] = op(a.data_[c], b.vec_);
      } else if constexpr (C == 1) {
        static_assert(R != 1 && R == R2);
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t r = 0; r < R2; ++r)
          z.data_[r] = op(a.data_[r], b.data_[r]);
      } else {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t r = 0; r < R2; ++r) {
          POLYMATHFULLUNROLL
          for (std::ptrdiff_t c = 0; c < C; ++c)
            if constexpr (R == R2) z[r, c] = op(a[r, c], b.data_[r]);
            else z[r, c] = op(a.data_[c], b.data_[r]);
        }
      }
      return z;
    }
  }
};
template <std::ptrdiff_t N, typename T> struct Unroll<1, 1, N, T> {
  static constexpr std::ptrdiff_t W =
    std::ptrdiff_t(std::bit_ceil(std::size_t(N)));
  using VT = Vec<W, T>;
  VT vec_;
  TRIVIAL constexpr auto operator[](std::ptrdiff_t) -> VT & { return vec_; }
  TRIVIAL constexpr auto operator[](std::ptrdiff_t, std::ptrdiff_t) -> VT & {
    return vec_;
  }
  TRIVIAL constexpr auto operator[](std::ptrdiff_t) const -> VT { return vec_; }
  TRIVIAL constexpr auto operator[](std::ptrdiff_t, std::ptrdiff_t) const
    -> VT {
    return vec_;
  }
  TRIVIAL constexpr operator VT() const { return vec_; }
  template <typename U>
  TRIVIAL constexpr operator Unroll<1, 1, N, U>() const
  requires(!std::same_as<T, U>)
  {
    if constexpr (W == 1) return {U(vec_)};
    else return {__builtin_convertvector(vec_, Vec<W, U>)};
  }
  TRIVIAL constexpr auto operator-() { return Unroll{-vec_}; }
  TRIVIAL constexpr auto operator+=(const Unroll &a) -> Unroll & {
    vec_ += a.vec_;
    return *this;
  }
  TRIVIAL constexpr auto operator-=(const Unroll &a) -> Unroll & {
    vec_ -= a.vec_;
    return *this;
  }
  TRIVIAL constexpr auto operator*=(const Unroll &a) -> Unroll & {
    vec_ *= a.vec_;
    return *this;
  }
  TRIVIAL constexpr auto operator/=(const Unroll &a) -> Unroll & {
    vec_ /= a.vec_;
    return *this;
  }
  TRIVIAL constexpr auto operator+=(VT a) -> Unroll & {
    vec_ += a;
    return *this;
  }
  TRIVIAL constexpr auto operator-=(VT a) -> Unroll & {
    vec_ -= a;
    return *this;
  }
  TRIVIAL constexpr auto operator*=(VT a) -> Unroll & {
    vec_ *= a;
    return *this;
  }
  TRIVIAL constexpr auto operator/=(VT a) -> Unroll & {
    vec_ /= a;
    return *this;
  }
  TRIVIAL constexpr auto operator+=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    vec_ += vbroadcast<W, T>(a);
    return *this;
  }
  TRIVIAL constexpr auto operator-=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    vec_ -= vbroadcast<W, T>(a);
    return *this;
  }
  TRIVIAL constexpr auto operator*=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    vec_ *= vbroadcast<W, T>(a);
    return *this;
  }
  TRIVIAL constexpr auto operator/=(std::convertible_to<T> auto a) -> Unroll &
  requires(W != 1)
  {
    vec_ /= vbroadcast<W, T>(a);
    return *this;
  }

private:
  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator+(Unroll a, Unroll<R1, C1, W1, T1> b) {
    return applyop(a, b, std::plus<>{});
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator-(Unroll a, Unroll<R1, C1, W1, T1> b) {
    return applyop(a, b, std::minus<>{});
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator*(Unroll a,
                                          const Unroll<R1, C1, W1, T1> &b) {
    return applyop(a, b, std::multiplies<>{});
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1>
  TRIVIAL friend constexpr auto operator/(Unroll a,
                                          const Unroll<R1, C1, W1, T1> &b) {
    return applyop(a, b, std::divides<>{});
  }

  TRIVIAL friend constexpr auto operator<(Unroll a, Unroll b) {
    return cmp::lt<N, T>(a.vec_, b.vec_);
  }
  TRIVIAL friend constexpr auto operator>(Unroll a, Unroll b) {
    return cmp::gt<N, T>(a.vec_, b.vec_);
  }
  TRIVIAL friend constexpr auto operator<=(Unroll a, Unroll b) {
    return cmp::le<N, T>(a.vec_, b.vec_);
  }
  TRIVIAL friend constexpr auto operator>=(Unroll a, Unroll b) {
    return cmp::ge<N, T>(a.vec_, b.vec_);
  }
  TRIVIAL friend constexpr auto operator==(Unroll a, Unroll b) {
    return cmp::eq<N, T>(a.vec_, b.vec_);
  }
  TRIVIAL friend constexpr auto operator!=(Unroll a, Unroll b) {
    return cmp::ne<N, T>(a.vec_, b.vec_);
  }
  TRIVIAL friend constexpr auto operator<(Unroll a, VT b) {
    return cmp::lt<N, T>(a.vec_, b);
  }
  TRIVIAL friend constexpr auto operator>(Unroll a, VT b) {
    return cmp::gt<N, T>(a.vec_, b);
  }
  TRIVIAL friend constexpr auto operator<=(Unroll a, VT b) {
    return cmp::le<N, T>(a.vec_, b);
  }
  TRIVIAL friend constexpr auto operator>=(Unroll a, VT b) {
    return cmp::ge<N, T>(a.vec_, b);
  }
  TRIVIAL friend constexpr auto operator==(Unroll a, VT b) {
    return cmp::eq<N, T>(a.vec_, b);
  }
  TRIVIAL friend constexpr auto operator!=(Unroll a, VT b) {
    return cmp::ne<N, T>(a.vec_, b);
  }
  TRIVIAL friend constexpr auto operator<(VT a, Unroll b) {
    return cmp::lt<N, T>(a, b.vec_);
  }
  TRIVIAL friend constexpr auto operator>(VT a, Unroll b) {
    return cmp::gt<N, T>(a, b.vec_);
  }
  TRIVIAL friend constexpr auto operator<=(VT a, Unroll b) {
    return cmp::le<N, T>(a, b.vec_);
  }
  TRIVIAL friend constexpr auto operator>=(VT a, Unroll b) {
    return cmp::ge<N, T>(a, b.vec_);
  }
  TRIVIAL friend constexpr auto operator==(VT a, Unroll b) {
    return cmp::eq<N, T>(a, b.vec_);
  }
  TRIVIAL friend constexpr auto operator!=(VT a, Unroll b) {
    return cmp::ne<N, T>(a, b.vec_);
  }

  TRIVIAL friend constexpr auto operator+(Unroll a, VT b) -> Unroll {
    return {a.vec_ + b};
  }
  TRIVIAL friend constexpr auto operator-(Unroll a, VT b) -> Unroll {
    return {a.vec_ - b};
  }
  TRIVIAL friend constexpr auto operator*(Unroll a, VT b) -> Unroll {
    return {a.vec_ * b};
  }
  TRIVIAL friend constexpr auto operator/(Unroll a, VT b) -> Unroll {
    return {a.vec_ / b};
  }
  TRIVIAL friend constexpr auto operator+(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a + vbroadcast<W, T>(b);
  }
  TRIVIAL friend constexpr auto operator-(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a - vbroadcast<W, T>(b);
  }
  TRIVIAL friend constexpr auto operator*(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a * vbroadcast<W, T>(b);
  }
  TRIVIAL friend constexpr auto operator/(Unroll a,
                                          std::convertible_to<T> auto b)
    -> Unroll
  requires(W != 1)
  {
    return a / vbroadcast<W, T>(b);
  }

  TRIVIAL friend constexpr auto operator+(VT a, Unroll b) -> Unroll {
    return {a + b.vec_};
  }
  TRIVIAL friend constexpr auto operator-(VT a, Unroll b) -> Unroll {
    return {a - b.vec_};
  }
  TRIVIAL friend constexpr auto operator*(VT a, Unroll b) -> Unroll {
    return {a * b.vec_};
  }
  TRIVIAL friend constexpr auto operator/(VT a, Unroll b) -> Unroll {
    return {a / b.vec_};
  }
  TRIVIAL friend constexpr auto operator+(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) + a;
  }
  TRIVIAL friend constexpr auto operator-(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) - a;
  }
  TRIVIAL friend constexpr auto operator*(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) * a;
  }
  TRIVIAL friend constexpr auto operator/(T b, Unroll a) -> Unroll
  requires(W != 1)
  {
    return vbroadcast<W, T>(b) / a;
  }

  template <std::ptrdiff_t R1, std::ptrdiff_t C1, std::ptrdiff_t W1,
            typename T1, typename Op>
  TRIVIAL friend constexpr auto applyop(Unroll a, Unroll<R1, C1, W1, T1> b,
                                        Op op) {
    // Possibilities:
    // 1. All match
    // 2. We had separate unrolls across rows and columns, and some arrays
    // were indexed by one or two of them.
    // In the latter case, we could have arrays indexed by rows, cols, or both.
    if constexpr (!std::same_as<T, T1>) {
      using PT = std::common_type_t<T, T1>;
      return applyop(Unroll<1, 1, W, PT>(a), Unroll<R1, C1, W1, PT>(b), op);
    } else if constexpr (W == W1) {
      // both were indexed by cols, and `C`s should also match
      // or neither were, and they should still match.
      static_assert(C1 == 1);
      if constexpr (R1 != 1) {
        // `a` was indexed across cols only
        Unroll<R1, 1, W, T> z;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t r = 0; r < R1; ++r)
          z.data_[r] = op(a.vec_, b.data_[r]);
        return z;
      } else return Unroll<1, 1, W, T>{op(a.vec_, b.vec_)};
    } else if constexpr (W == 1) {
      // `a` was indexed by row only
      static_assert(R1 == 1);
      if constexpr (C1 != 1) {
        Unroll<1, C1, W1, T> z;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C1; ++c)
          z.data_[c] = op(a.vec_, b.data_[c]);
        return z;
      } else return Unroll<1, 1, W1, T>{op(a.vec_, b.vec_)};
    } else {
      static_assert(W1 == 1);
      static_assert(R1 == 1 || C1 == 1);
      constexpr std::ptrdiff_t R = R1 == 1 ? C1 : R1;
      // `b` was indexed by row only
      if constexpr (R != 1) {
        Unroll<R, 1, W, T> z;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t r = 0; r < R; ++r)
          z.data_[r] = op(a.vec_, b.data_[r]);
        return z;
      } else return Unroll{op(a.vec_, b.vec_)};
    }
  }
};

template <std::ptrdiff_t R, std::ptrdiff_t C, std::ptrdiff_t N, typename T,
          std::ptrdiff_t X, std::size_t NM, typename MT = mask::None<N>>
TRIVIAL constexpr auto loadunroll(const T *ptr, math::RowStride<X> rowStride,
                                  std::array<MT, NM> masks)
  -> Unroll<R, C, N, T> {
  if constexpr (R * C == 1) {
    MT msk = masks[0];
    Vec<N, T> x = load<T>(ptr, msk);
    return {x};
    // return {load(ptr, masks[0])};
  } else {
    constexpr auto W = std::ptrdiff_t(std::bit_ceil(std::size_t(N)));
    auto rs = std::ptrdiff_t(rowStride);
    Unroll<R, C, N, T> ret;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t r = 0; r < R; ++r, ptr += rs) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          ret[r, c] = load(ptr + c * W, mask::None<W>{});
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          ret[r, c] = load(ptr + c * W, masks[c]);
      } else {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          ret[r, c] = load(ptr + c * W, mask::None<W>{});
        ret[r, C - 1] = load(ptr + (C - 1) * W, masks[0]);
      }
    }
    return ret;
  }
}
template <std::ptrdiff_t R, std::ptrdiff_t C, std::ptrdiff_t N, typename T,
          std::ptrdiff_t X, std::size_t NM, typename MT = mask::None<N>>
TRIVIAL constexpr auto loadstrideunroll(const T *ptr,
                                        math::RowStride<X> rowStride,
                                        std::array<MT, NM> masks)
  -> Unroll<R, C, N, T> {
  auto s = std::int32_t(std::ptrdiff_t(rowStride));
  if constexpr (R * C == 1) return {load(ptr, masks[0], s)};
  else {
    constexpr auto W = std::ptrdiff_t(std::bit_ceil(std::size_t(N)));
    Unroll<R, C, N, T> ret;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t r = 0; r < R; ++r, ++ptr) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          ret[r, c] = load(ptr + c * W * s, mask::None<W>{}, s);
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          ret[r, c] = load(ptr + c * W * s, masks[c], s);
      } else {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          ret[r, c] = load(ptr + c * W * s, mask::None<W>{}, s);
        ret[r, C - 1] = load(ptr + (C - 1) * W * s, masks[0], s);
      }
    }
    return ret;
  }
}

// Represents a reference for a SIMD load, in particular so that we can store.
// Needs to support masking
template <std::ptrdiff_t R, std::ptrdiff_t C, std::ptrdiff_t N, typename T,
          std::ptrdiff_t X, std::ptrdiff_t NM, typename MT = mask::None<N>,
          bool Transposed = false>
struct UnrollRef {
  static constexpr std::ptrdiff_t W =
    std::ptrdiff_t(std::bit_ceil(std::size_t(N)));
  static_assert(N == W || C == 1,
                "If N != the next power of `2`, then `C` should be `1`");
  static_assert(
    NM == 0 || NM == 1 || NM == C,
    "Should have no masks, one mask for last `C`, or one mask per `C`");
  using UT = Unroll<R, C, N, T>;
  T *ptr_;
  [[no_unique_address]] math::RowStride<X> row_stride_;
  [[no_unique_address]] std::array<MT, NM> masks_;
  TRIVIAL constexpr operator UT() {
    if constexpr (!Transposed)
      return loadunroll<R, C, N, T, X, NM, MT>(ptr_, row_stride_, masks_);
    else
      return loadstrideunroll<R, C, N, T, X, NM, MT>(ptr_, row_stride_, masks_);
  }
  TRIVIAL constexpr auto operator=(UT x) -> UnrollRef &
  requires(!Transposed)
  {
    auto rs = std::ptrdiff_t(row_stride_);
    T *p = ptr_;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t r = 0; r < R; ++r, p += rs) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W, mask::None<W>{}, x[r, c]);
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c) {
          auto msk = masks_[c];
          store<T>(p + c * W, msk, x[r, c]);
          // store<T>(p + c * W, masks[c], x[r, c]);
        }
      } else { // NM == 1
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          store<T>(p + c * W, mask::None<W>{}, x[r, c]);
        store<T>(p + (C - 1) * W, masks_[0], x[r, C - 1]);
      }
    }
    return *this;
  }
  TRIVIAL constexpr auto operator=(Unroll<1, C, N, T> x) -> UnrollRef &
  requires((!Transposed) && (R != 1))
  {
    auto rs = std::ptrdiff_t(row_stride_);
    T *p = ptr_;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t r = 0; r < R; ++r, p += rs) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W, mask::None<W>{}, x[0, c]);
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c) {
          auto msk = masks_[c];
          store<T>(p + c * W, msk, x[0, c]);
        }
      } else { // NM == 1
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          store<T>(p + c * W, mask::None<W>{}, x[0, c]);
        store<T>(p + (C - 1) * W, masks_[0], x[0, C - 1]);
      }
    }
    return *this;
  }
  TRIVIAL constexpr auto operator=(Unroll<R, C, 1, T> x) -> UnrollRef &
  requires(!Transposed && (N != 1))
  {
    auto rs = std::ptrdiff_t(row_stride_);
    T *p = ptr_;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t r = 0; r < R; ++r, p += rs) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W, mask::None<W>{}, vbroadcast<W>(x[r, c]));
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W, masks_[c], vbroadcast<W>(x[r, c]));
      } else { // NM == 1
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          store<T>(p + c * W, mask::None<W>{}, x[r, c]);
        store<T>(p + (C - 1) * W, masks_[0], vbroadcast<W>(x[r, C - 1]));
      }
    }
    return *this;
  }

  TRIVIAL constexpr auto operator=(Vec<W, T> v) -> UnrollRef &
  requires(!Transposed)
  {
    auto rs = std::ptrdiff_t(row_stride_);
    T *p = ptr_;
    POLYMATHFULLUNROLL
    for (std::ptrdiff_t r = 0; r < R; ++r, p += rs) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W, mask::None<W>{}, v);
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W, masks_[c], v);
      } else { // NM == 1
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          store<T>(p + c * W, mask::None<W>{}, v);
        store<T>(p + (C - 1) * W, masks_[0], v);
      }
    }
    return *this;
  }
  TRIVIAL constexpr auto operator=(Unroll<R, C, N, T> x) -> UnrollRef &
  requires(Transposed)
  {
    auto s = std::int32_t(std::ptrdiff_t(row_stride_));
    T *p = ptr_;
    for (std::ptrdiff_t r = 0; r < R; ++r, ++p) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W * s, mask::None<W>{}, x[0, c], s);
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W * s, masks_[c], x[0, c], s);
      } else {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          store<T>(p + c * W * s, mask::None<W>{}, x[0, c], s);
        store<T>(p + (C - 1) * W * s, masks_[0], x[0, C - 1], s);
      }
    }
    return *this;
  }
  TRIVIAL constexpr auto operator=(Vec<W, T> v) -> UnrollRef &
  requires(Transposed)
  {
    auto s = std::int32_t(std::ptrdiff_t(row_stride_));
    T *p = ptr_;
    for (std::ptrdiff_t r = 0; r < R; ++r, ++p) {
      if constexpr (NM == 0) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W * s, mask::None<W>{}, v, s);
      } else if constexpr (NM == C) {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C; ++c)
          store<T>(p + c * W * s, masks_[c], v, s);
      } else {
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t c = 0; c < C - 1; ++c)
          store<T>(p + c * W * s, mask::None<W>{}, v, s);
        store<T>(p + (C - 1) * W * s, masks_[0], v, s);
      }
    }
    return *this;
  }
  TRIVIAL constexpr auto operator=(std::convertible_to<T> auto x)
    -> UnrollRef & {
    *this = Vec<W, T>{} + T(x);
    return *this;
  }
  TRIVIAL constexpr auto operator+=(const auto &x) -> UnrollRef & {
    return (*this) = Unroll<R, C, N, T>(*this) + x;
  }
  TRIVIAL constexpr auto operator-=(const auto &x) -> UnrollRef & {
    return (*this) = Unroll<R, C, N, T>(*this) - x;
  }
  TRIVIAL constexpr auto operator*=(const auto &x) -> UnrollRef & {
    return (*this) = Unroll<R, C, N, T>(*this) * x;
  }
  TRIVIAL constexpr auto operator/=(const auto &x) -> UnrollRef & {
    return (*this) = Unroll<R, C, N, T>(*this) / x;
  }
};
template <typename T, std::ptrdiff_t R, std::ptrdiff_t C, std::ptrdiff_t W,
          typename M, bool Transposed, std::ptrdiff_t X>
TRIVIAL constexpr auto ref(const T *p,
                           index::UnrollDims<R, C, W, M, Transposed, X> i)
  -> Unroll<R, C, W, T> {
  if constexpr (Transposed)
    return loadstrideunroll<R, C, W>(p, i.rs_, std::array<M, 1>{i.mask_});
  else return loadunroll<R, C, W>(p, i.rs_, std::array<M, 1>{i.mask_});
}
template <typename T, std::ptrdiff_t R, std::ptrdiff_t C, std::ptrdiff_t W,
          typename M, bool Transposed, std::ptrdiff_t X>
TRIVIAL constexpr auto ref(T *p, index::UnrollDims<R, C, W, M, Transposed, X> i)
  -> UnrollRef<R, C, W, T, X, 1, M, Transposed> {
  return {p, i.rs_, std::array<M, 1>{i.mask_}};
}

} // namespace simd

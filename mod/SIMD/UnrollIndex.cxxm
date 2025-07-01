#ifdef USE_MODULE
module;
#else
#pragma once
#endif
#include "LoopMacros.hxx"
#include "Macros.hxx"
#ifndef USE_MODULE
#include <concepts>
#include <cstddef>
#include <cstdint>

#include "SIMD/Indexing.cxx"
#include "SIMD/Intrin.cxx"
#include "SIMD/Masks.cxx"
#include "SIMD/Unroll.cxx"
#include "SIMD/Vec.cxx"
#else
export module SIMD:UnrollIndex;

import :Index;
import :Intrin;
import :Mask;
import :Unroll;
import :Vec;
import std;
#endif

#ifdef USE_MODULE
export namespace simd::index {
#else
namespace simd::index {
#endif
// Unroll rows by a factor of `R` and cols by `C`, vectorizing with width `W`
template <std::ptrdiff_t U, std::ptrdiff_t W = 1, typename M = mask::None<W>>
struct Unroll {
  std::ptrdiff_t index_;
  [[no_unique_address]] M mask_{};
  TRIVIAL explicit constexpr operator std::ptrdiff_t() const { return index_; }
  TRIVIAL explicit constexpr operator bool() const { return bool(mask_); }

  template <std::integral I>
  TRIVIAL constexpr operator ::simd::Unroll<1, U, W, I>() const {
    if constexpr (U != 1) {
      ::simd::Unroll<1, U, W, I> ret;
      Vec<W, I> r{::simd::range<W, I>()},
        ind = vbroadcast<W>(static_cast<I>(index_));
      POLYMATHFULLUNROLL
      for (std::ptrdiff_t u = 0; u < U; ++u) ret.data_[u] = r + ind + (W * u);
      return ret;
    } else return {::simd::range<W, I>() + index_};
  }
#ifdef __AVX512F__
  template <std::ptrdiff_t S>
  TRIVIAL constexpr auto sub() requires(std::same_as<M, mask::Bit<W>>) {
    static_assert((S <= W) && (U == 1));
    Unroll<1, S, mask::Bit<S>> u{index_, mask_.template sub<S>()};
    index_ += S;
    return u;
  }
#endif

private:
  TRIVIAL friend constexpr auto operator+(Unroll a, std::ptrdiff_t b)
    -> Unroll {
    return {b + a.index_, a.mask_};
  }
  TRIVIAL friend constexpr auto operator==(Unroll x, std::ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::simd::Unroll<U, 1, 1, std::int64_t> ret;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t u = 0; u < U; ++u)
          ret.data[u] = (x.index_ + u) == y;
        return ret;
      } else return ::simd::Unroll<1, 1, W, std::int64_t>{x.index_ == y};
    } else if constexpr (U > 1) {
      ::simd::Unroll<1, U, W, std::int64_t> ret;
      Vec<W, std::int64_t> v = vbroadcast<W, std::int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (std::ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, std::int64_t>() == (v - u * W);
      return ret;
    } else
      return ::simd::Unroll<1, 1, W, std::int64_t>{
        range<W, std::int64_t>() == vbroadcast<W, std::int64_t>(y - x.index_)};
  }
  TRIVIAL friend constexpr auto operator!=(Unroll x, std::ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::simd::Unroll<U, 1, 1, std::int64_t> ret;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t u = 0; u < U; ++u)
          ret.data[u] = (x.index_ + u) != y;
        return ret;
      } else return ::simd::Unroll<1, 1, W, std::int64_t>{x.index_ != y};
    } else if constexpr (U > 1) {
      ::simd::Unroll<1, U, W, std::int64_t> ret;
      Vec<W, std::int64_t> v = vbroadcast<W, std::int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (std::ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, std::int64_t>() != (v - u * W);
      return ret;
    } else
      return ::simd::Unroll<1, 1, W, std::int64_t>{
        range<W, std::int64_t>() != vbroadcast<W, std::int64_t>(y - x.index_)};
  }

  TRIVIAL friend constexpr auto operator<(Unroll x, std::ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::simd::Unroll<U, 1, 1, std::int64_t> ret;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) < y;
        return ret;
      } else return ::simd::Unroll<1, 1, W, std::int64_t>{x.index_ < y};
    } else if constexpr (U > 1) {
      ::simd::Unroll<1, U, W, std::int64_t> ret;
      Vec<W, std::int64_t> v = vbroadcast<W, std::int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (std::ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, std::int64_t>() < (v - u * W);
      return ret;
    } else
      return ::simd::Unroll<1, 1, W, std::int64_t>{
        range<W, std::int64_t>() < vbroadcast<W, std::int64_t>(y - x.index_)};
  }

  TRIVIAL friend constexpr auto operator>(Unroll x, std::ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::simd::Unroll<U, 1, 1, std::int64_t> ret;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) > y;
        return ret;
      } else return ::simd::Unroll<1, 1, W, std::int64_t>{x.index_ > y};
    } else if constexpr (U > 1) {
      ::simd::Unroll<1, U, W, std::int64_t> ret;
      Vec<W, std::int64_t> v = vbroadcast<W, std::int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (std::ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, std::int64_t>() > (v - u * W);
      return ret;
    } else
      return ::simd::Unroll<1, 1, W, std::int64_t>{
        range<W, std::int64_t>() > vbroadcast<W, std::int64_t>(y - x.index_)};
  }

  TRIVIAL friend constexpr auto operator<=(Unroll x, std::ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::simd::Unroll<U, 1, 1, std::int64_t> ret;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t u = 0; u < U; ++u)
          ret.data[u] = (x.index_ + u) <= y;
        return ret;
      } else return ::simd::Unroll<1, 1, W, std::int64_t>{x.index_ <= y};
    } else if constexpr (U > 1) {
      ::simd::Unroll<1, U, W, std::int64_t> ret;
      Vec<W, std::int64_t> v = vbroadcast<W, std::int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (std::ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, std::int64_t>() <= (v - u * W);
      return ret;
    } else
      return ::simd::Unroll<1, 1, W, std::int64_t>{
        range<W, std::int64_t>() <= vbroadcast<W, std::int64_t>(y - x.index_)};
  }

  TRIVIAL friend constexpr auto operator>=(Unroll x, std::ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::simd::Unroll<U, 1, 1, std::int64_t> ret;
        POLYMATHFULLUNROLL
        for (std::ptrdiff_t u = 0; u < U; ++u)
          ret.data[u] = (x.index_ + u) >= y;
        return ret;
      } else return ::simd::Unroll<1, 1, W, std::int64_t>{x.index_ >= y};
    } else if constexpr (U > 1) {
      ::simd::Unroll<1, U, W, std::int64_t> ret;
      Vec<W, std::int64_t> v = vbroadcast<W, std::int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (std::ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, std::int64_t>() >= (v - u * W);
      return ret;
    } else
      return ::simd::Unroll<1, 1, W, std::int64_t>{
        range<W, std::int64_t>() >= vbroadcast<W, std::int64_t>(y - x.index_)};
  }
};
template <std::ptrdiff_t U, std::ptrdiff_t W>
TRIVIAL constexpr auto unrollmask(std::ptrdiff_t L, std::ptrdiff_t i) {
  // mask applies to last iter
  // We can't check that the last iter is non-empty, because that
  // could be the loop exit condition
  auto m{mask::create<W>(i + (U - 1) * W, L)};
  return Unroll<U, W, decltype(m)>{i, m};
};
#ifdef __AVX512VL__
template <std::ptrdiff_t W>
TRIVIAL constexpr auto tailmask(std::ptrdiff_t i, std::ptrdiff_t m)
  -> Unroll<1, W, mask::Bit<W>> {
  return {i, mask::createSmallPositive<W>(m)};
}
#else
template <std::ptrdiff_t W>
TRIVIAL constexpr auto tailmask(std::ptrdiff_t i, std::ptrdiff_t m) {
  auto mask{mask::create<W>(i, i + m)};
  return Unroll<1, W, decltype(mask)>{i, mask};
}
#endif
template <std::ptrdiff_t U, std::ptrdiff_t W, typename M>
inline constexpr bool issimd<Unroll<U, W, M>> = true;
} // namespace simd::index

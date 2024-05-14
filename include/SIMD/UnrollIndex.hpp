#pragma once
#include "SIMD/Unroll.hpp"

namespace poly::simd::index {
// Unroll rows by a factor of `R` and cols by `C`, vectorizing with width `W`
template <ptrdiff_t U, ptrdiff_t W = 1, typename M = mask::None<W>>
struct Unroll {
  ptrdiff_t index;
  [[no_unique_address]] M mask{};
  explicit constexpr operator ptrdiff_t() const { return index; }
  explicit constexpr operator bool() const { return bool(mask); }
  constexpr auto operator+(ptrdiff_t b) -> Unroll { return {b + index, mask}; }

  friend constexpr auto operator==(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index + u) == y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index == y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() == (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() == vbroadcast<W, int64_t>(y - x.index)};
  }
  friend constexpr auto operator!=(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index + u) != y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index != y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() != (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() != vbroadcast<W, int64_t>(y - x.index)};
  }

  friend constexpr auto operator<(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index + u) < y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index < y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() < (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() < vbroadcast<W, int64_t>(y - x.index)};
  }

  friend constexpr auto operator>(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index + u) > y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index > y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() > (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() > vbroadcast<W, int64_t>(y - x.index)};
  }

  friend constexpr auto operator<=(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index + u) <= y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index <= y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() <= (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() <= vbroadcast<W, int64_t>(y - x.index)};
  }

  friend constexpr auto operator>=(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index + u) >= y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index >= y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() >= (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() >= vbroadcast<W, int64_t>(y - x.index)};
  }
};
template <ptrdiff_t U, ptrdiff_t W>
[[gnu::always_inline]] constexpr auto unrollmask(ptrdiff_t L, ptrdiff_t i) {
  // mask applies to last iter
  // We can't check that the last iter is non-empty, because that
  // could be the loop exit condition
  auto m{mask::create<W>(i + (U - 1) * W, L)};
  return Unroll<U, W, decltype(m)>{i, m};
};
template <ptrdiff_t U, ptrdiff_t W, typename M>
static constexpr bool issimd<Unroll<U, W, M>> = true;
} // namespace poly::simd::index


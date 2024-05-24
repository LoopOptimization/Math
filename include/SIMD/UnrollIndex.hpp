#pragma once
#include "SIMD/Indexing.hpp"
#include "SIMD/Intrin.hpp"
#include "SIMD/Masks.hpp"
#include "SIMD/Unroll.hpp"
#include "SIMD/Vec.hpp"
#include "Utilities/LoopMacros.hpp"
#include <concepts>
#include <cstddef>
#include <cstdint>

namespace poly::simd::index {
// Unroll rows by a factor of `R` and cols by `C`, vectorizing with width `W`
template <ptrdiff_t U, ptrdiff_t W = 1, typename M = mask::None<W>>
struct Unroll {
  ptrdiff_t index_;
  [[no_unique_address]] M mask_{};
  explicit constexpr operator ptrdiff_t() const { return index_; }
  explicit constexpr operator bool() const { return bool(mask_); }

  template <std::integral I>
  constexpr operator ::poly::simd::Unroll<1, U, W, I>() const {
    if constexpr (U != 1) {
      ::poly::simd::Unroll<1, U, W, I> ret;
      Vec<W, I> r{::poly::simd::range<W, I>()},
        ind = vbroadcast<W>(static_cast<I>(index_));
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u) ret.data_[u] = r + ind + (W * u);
      return ret;
    } else return {::poly::simd::range<W, I>() + index_};
  }

private:
  friend constexpr auto operator+(Unroll a, ptrdiff_t b) -> Unroll {
    return {b + a.index_, a.mask_};
  }
  friend constexpr auto operator==(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) == y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index_ == y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() == (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() == vbroadcast<W, int64_t>(y - x.index_)};
  }
  friend constexpr auto operator!=(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) != y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index_ != y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() != (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() != vbroadcast<W, int64_t>(y - x.index_)};
  }

  friend constexpr auto operator<(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) < y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index_ < y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() < (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() < vbroadcast<W, int64_t>(y - x.index_)};
  }

  friend constexpr auto operator>(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) > y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index_ > y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() > (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() > vbroadcast<W, int64_t>(y - x.index_)};
  }

  friend constexpr auto operator<=(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) <= y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index_ <= y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() <= (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() <= vbroadcast<W, int64_t>(y - x.index_)};
  }

  friend constexpr auto operator>=(Unroll x, ptrdiff_t y) {
    if constexpr (W == 1) {
      if constexpr (U > 1) {
        ::poly::simd::Unroll<U, 1, 1, int64_t> ret;
        POLYMATHFULLUNROLL
        for (ptrdiff_t u = 0; u < U; ++u) ret.data[u] = (x.index_ + u) >= y;
        return ret;
      } else return ::poly::simd::Unroll<1, 1, W, int64_t>{x.index_ >= y};
    } else if constexpr (U > 1) {
      ::poly::simd::Unroll<1, U, W, int64_t> ret;
      Vec<W, int64_t> v = vbroadcast<W, int64_t>(y - x.index_);
      POLYMATHFULLUNROLL
      for (ptrdiff_t u = 0; u < U; ++u)
        ret.data[u] = range<W, int64_t>() >= (v - u * W);
      return ret;
    } else
      return ::poly::simd::Unroll<1, 1, W, int64_t>{
        range<W, int64_t>() >= vbroadcast<W, int64_t>(y - x.index_)};
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

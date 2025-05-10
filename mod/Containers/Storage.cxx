#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#ifndef USE_MODULE
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "Math/AxisTypes.cxx"
#include "Math/MatrixDimensions.cxx"
#else
export module Storage;

import MatDim;
import AxisTypes;
import std;
#endif

#ifdef USE_MODULE
export namespace containers {
#else
namespace containers {
#endif
namespace detail {
template <class T>
concept SizeMultiple8 = (sizeof(T) % 8) == 0;

template <class S> struct DefaultCapacityType {
  using type = math::Capacity<-1, int>;
};
template <SizeMultiple8 S> struct DefaultCapacityType<S> {
  using type = math::Capacity<-1, std::ptrdiff_t>;
};
static_assert(!SizeMultiple8<std::uint32_t>);
static_assert(SizeMultiple8<std::uint64_t>);

consteval auto log2Floor(std::uint64_t x) -> std::uint64_t {
  return 63 - std::countl_zero(x);
}
consteval auto log2Ceil(std::uint64_t x) -> std::uint64_t {
  return 64 - std::countl_zero(x - 1);
}
// NOLINTNEXTLINE(misc-no-recursion)
consteval auto bisectFindSquare(std::uint64_t l, std::uint64_t h,
                                std::uint64_t N) -> std::uint64_t {
  if (l == h) return l;
  std::uint64_t m = (l + h) / 2;
  if (m * m >= N) return bisectFindSquare(l, m, N);
  return bisectFindSquare(m + 1, h, N);
}

} // namespace detail

template <class S>
using default_capacity_type_t = typename detail::DefaultCapacityType<S>::type;
static_assert(sizeof(default_capacity_type_t<std::uint32_t>) == 4);
static_assert(sizeof(default_capacity_type_t<std::uint64_t>) == 8);

template <typename T, std::ptrdiff_t N> struct Storage {
  static_assert(N > 0);
  // We can avoid `reinterpret_cast` if we have trivial/implicit lifetime types.
  static constexpr bool trivial =
    std::is_trivially_default_constructible_v<T> &&
    std::is_trivially_destructible_v<T>;
  static constexpr std::ptrdiff_t NumElt = trivial ? N : N * sizeof(T);
  using DataElt = std::conditional_t<trivial, T, char>;
  alignas(T) DataElt mem[NumElt]; // NOLINT (modernize-avoid-c-style-arrays)
  constexpr auto data() -> T * {
    if constexpr (trivial) return mem;
    else return reinterpret_cast<T *>(mem);
  }
  constexpr auto data() const -> const T * {
    if constexpr (trivial) return mem;
    else return reinterpret_cast<const T *>(mem);
  }
  constexpr Storage() {} // NOLINT (modernize-use-equals-default)
};
template <typename T> struct Storage<T, 0> {
  static constexpr auto data() -> T * { return nullptr; }
};

template <class T, class S> consteval auto PreAllocStorage() -> std::ptrdiff_t {
  static constexpr std::ptrdiff_t total_bytes = 128;
  static constexpr std::ptrdiff_t nrow = S::nrow;
  static constexpr std::ptrdiff_t nstride = S::nstride;
  // constexpr std::ptrdiff_t remainingBytes =
  //   totalBytes - sizeof(T *) - sizeof(S) -
  //   sizeof(default_capacity_type_t<S>);
  // constexpr std::ptrdiff_t N = remainingBytes / std::ptrdiff_t(sizeof(T));
  constexpr std::ptrdiff_t N = total_bytes / std::ptrdiff_t(sizeof(T));
  static_assert(N <= 128);
  if constexpr (nrow > 0 && nstride > 0) return nrow * nstride;
  else if constexpr (N <= 0) return 0;
  // else if constexpr (!math::MatrixDimension<S>) return N;
  else if constexpr (std::convertible_to<S, math::SquareDims<>>) {
    constexpr auto UN = std::uint64_t(N);
    // a fairly naive algorirthm for computing the next square `N`
    // sqrt(x) = x^(1/2) = exp2(log2(x)/2)
    constexpr std::uint64_t R = detail::log2Floor(UN) / 2;
    static_assert(R < 63);
    constexpr std::uint64_t L = std::uint64_t(1) << R;
    constexpr std::uint64_t H = std::uint64_t(1) << ((detail::log2Ceil(N) + 1) / 2);
    return std::ptrdiff_t(detail::bisectFindSquare(L, H, UN));
  } else if (nrow > 0) {
    return (N / nrow) * nrow;
  } else if (nstride > 0) {
    return (N / nstride) * nstride;
  } else return N;
}
} // namespace containers

module;

#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>

export module Storage;

import MatDim;
import AxisTypes;

template <class T>
concept SizeMultiple8 = (sizeof(T) % 8) == 0;

template <class S> struct DefaultCapacityType {
  using type = math::Capacity<-1, int>;
};
template <SizeMultiple8 S> struct DefaultCapacityType<S> {
  using type = math::Capacity<-1, std::ptrdiff_t>;
};
template <class S>
using default_capacity_type_t = typename DefaultCapacityType<S>::type;

static_assert(!SizeMultiple8<uint32_t>);
static_assert(SizeMultiple8<uint64_t>);
static_assert(sizeof(default_capacity_type_t<uint32_t>) == 4);
static_assert(sizeof(default_capacity_type_t<uint64_t>) == 8);

consteval auto log2Floor(uint64_t x) -> uint64_t {
  return 63 - std::countl_zero(x);
}
consteval auto log2Ceil(uint64_t x) -> uint64_t {
  return 64 - std::countl_zero(x - 1);
}
// NOLINTNEXTLINE(misc-no-recursion)
consteval auto bisectFindSquare(uint64_t l, uint64_t h,
                                uint64_t N) -> uint64_t {
  if (l == h) return l;
  uint64_t m = (l + h) / 2;
  if (m * m >= N) return bisectFindSquare(l, m, N);
  return bisectFindSquare(m + 1, h, N);
}

export namespace containers {
template <typename T, ptrdiff_t N> struct Storage {
  static_assert(N > 0);
  // We can avoid `reinterpret_cast` if we have trivial/implicit lifetime types.
  static constexpr bool trivial =
    std::is_trivially_default_constructible_v<T> &&
    std::is_trivially_destructible_v<T>;
  static constexpr ptrdiff_t NumElt = trivial ? N : N * sizeof(T);
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

template <class T, class S> consteval auto PreAllocStorage() -> ptrdiff_t {
  static constexpr ptrdiff_t total_bytes = 128;
  // constexpr ptrdiff_t remainingBytes =
  //   totalBytes - sizeof(T *) - sizeof(S) -
  //   sizeof(default_capacity_type_t<S>);
  // constexpr ptrdiff_t N = remainingBytes / ptrdiff_t(sizeof(T));
  constexpr ptrdiff_t N = total_bytes / ptrdiff_t(sizeof(T));
  static_assert(N <= 128);
  if constexpr (N <= 0) return 0;
  // else if constexpr (!math::MatrixDimension<S>) return N;
  else if constexpr (std::convertible_to<S, math::SquareDims<>>) {
    constexpr auto UN = uint64_t(N);
    // a fairly naive algorirthm for computing the next square `N`
    // sqrt(x) = x^(1/2) = exp2(log2(x)/2)
    constexpr uint64_t R = log2Floor(UN) / 2;
    static_assert(R < 63);
    constexpr uint64_t L = uint64_t(1) << R;
    constexpr uint64_t H = uint64_t(1) << ((log2Ceil(N) + 1) / 2);
    return ptrdiff_t(bisectFindSquare(L, H, UN));
  } else return N;
}
} // namespace containers

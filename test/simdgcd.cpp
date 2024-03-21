
#include "Containers/Tuple.hpp"
#include "Math/GreatestCommonDivisor.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include <limits>
#include <numeric>
#include <random>

template <ptrdiff_t W>
constexpr auto vecget(poly::simd::Vec<W, int64_t> v, ptrdiff_t i) -> int64_t {
  return v[i];
}
template <> constexpr auto vecget<1>(int64_t v, ptrdiff_t) -> int64_t {
  return v;
}
template <ptrdiff_t W>
constexpr void vecset(poly::simd::Vec<W, int64_t> &v, int64_t newval,
                      ptrdiff_t i) {
  v[i] = newval;
}
template <> constexpr void vecset<1>(int64_t &v, int64_t newval, ptrdiff_t) {
  v = newval;
}

template <ptrdiff_t W>
inline auto fillGCD(std::mt19937 &gen)
  -> poly::containers::Tuple<poly::simd::Vec<W, int64_t>,
                             poly::simd::Vec<W, int64_t>,
                             std::array<int64_t, size_t(W)>, int64_t> {
  std::array<int64_t, size_t(W)> z;
  poly::simd::Vec<W, int64_t> a, b;
  std::uniform_int_distribution<> distrib(std::numeric_limits<int>::min(),
                                          std::numeric_limits<int>::max());
  int64_t rg = 0;
  for (ptrdiff_t j = 0; j < W; ++j) {
    int64_t w = distrib(gen), x = distrib(gen), y = std::gcd(w, x);
    EXPECT_EQ(y, poly::math::gcd(w, x));
    z[j] = y;
    vecset<W>(a, w, j);
    vecset<W>(b, x, j);
    rg = std::gcd(rg, y);
  }
  return {a, b, z, rg};
}

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(SIMDGCDTest, BasicAssertions) {
  constexpr ptrdiff_t W = poly::simd::Width<int64_t>;
  std::random_device rd;
  std::mt19937 gen(rd());
  for (ptrdiff_t i = 0; i < 20000; ++i) {
    auto [a, b, z, rg] = fillGCD<W>(gen);
    poly::simd::Vec<W, int64_t> g = poly::math::gcd<W>(a, b);
    for (ptrdiff_t j = 0; j < W; ++j) EXPECT_EQ(vecget<W>(g, j), z[j]);
    EXPECT_EQ(rg, poly::math::gcdreduce<W>(g));
  }
}

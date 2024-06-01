#include "Math/MultiplicativeInverse.hpp"
#include <cstdint>
#include <gtest/gtest.h>

using poly::math::MultiplicativeInverse;

TEST(MultiplicativeInverse, BasicAssertions) {
  for (int32_t j = -100; j <= 100; ++j) {
    if (j == 0) continue;
    auto mij = MultiplicativeInverse(j);
    for (int32_t i = -1000; i <= 1000; ++i) {
      auto [d, r] = mij.divrem(i);
      EXPECT_EQ(i / j, d);
      EXPECT_EQ(i % j, r);
    }
  }
  for (int64_t j = -100; j <= 100; ++j) {
    if (j == 0) continue;
    auto mij = MultiplicativeInverse(j);
    for (int64_t i = -1000; i <= 1000; ++i) {
      auto [d, r] = mij.divrem(i);
      EXPECT_EQ(i / j, d);
      EXPECT_EQ(i % j, r);
    }
  }
  for (uint32_t j = 1; j <= 200; ++j) {
    auto mij = MultiplicativeInverse(j);
    for (uint32_t i = 0; i <= 2000; ++i) {
      auto [d, r] = mij.divrem(i);
      EXPECT_EQ(i / j, d);
      EXPECT_EQ(i % j, r);
    }
  }
  for (uint64_t j = 1; j <= 200; ++j) {
    auto mij = MultiplicativeInverse(j);
    for (uint64_t i = 0; i <= 2000; ++i) {
      auto [d, r] = mij.divrem(i);
      EXPECT_EQ(i / j, d);
      EXPECT_EQ(i % j, r);
    }
  }
}

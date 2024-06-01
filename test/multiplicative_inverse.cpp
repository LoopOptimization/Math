#include "Math/MultiplicativeInverse.hpp"
#include <cstdint>
#include <gtest/gtest.h>

using poly::math::MultiplicativeInverse, poly::math::cld;

TEST(MultiplicativeInverse, BasicAssertions) {
  for (int32_t j = -100; j <= 100; ++j) {
    if (j == 0) continue;
    auto mij = MultiplicativeInverse(j);
    auto mijf = MultiplicativeInverse(float(j));
    for (int32_t i = -1000; i <= 1000; ++i) {
      auto [d, r] = mij.divrem(i);
      int32_t qref = i / j, rref = i % j, cref = cld(i, j);
      EXPECT_EQ(qref, d);
      EXPECT_EQ(rref, r);
      EXPECT_EQ(cref, cld(i, mij));
      auto fi = float(i);
      auto [df, rf] = mijf.divrem(fi);
      EXPECT_EQ(qref, df);
      EXPECT_EQ(rref, rf);
      EXPECT_EQ(cref, cld(float(i), mijf));
    }
  }
  for (int64_t j = -100; j <= 100; ++j) {
    if (j == 0) continue;
    auto mij = MultiplicativeInverse(j);
    auto mijf = MultiplicativeInverse(double(j));
    for (int64_t i = -1000; i <= 1000; ++i) {
      auto [d, r] = mij.divrem(i);
      int64_t qref = i / j, rref = i % j, cref = cld(i, j);
      EXPECT_EQ(qref, d);
      EXPECT_EQ(rref, r);
      EXPECT_EQ(cref, cld(i, mij));
      auto [df, rf] = mijf.divrem(double(i));
      EXPECT_EQ(qref, df);
      EXPECT_EQ(rref, rf);
      EXPECT_EQ(cref, cld(float(i), mijf));
    }
  }
  for (uint32_t j = 1; j <= 200; ++j) {
    auto mij = MultiplicativeInverse(j);
    for (uint32_t i = 0; i <= 2000; ++i) {
      auto [d, r] = mij.divrem(i);
      uint32_t qref = i / j, rref = i % j;
      EXPECT_EQ(qref, d);
      EXPECT_EQ(rref, r);
    }
  }
  for (uint64_t j = 1; j <= 200; ++j) {
    auto mij = MultiplicativeInverse(j);
    for (uint64_t i = 0; i <= 2000; ++i) {
      auto [d, r] = mij.divrem(i);
      uint64_t qref = i / j, rref = i % j;
      EXPECT_EQ(qref, d);
      EXPECT_EQ(rref, r);
    }
  }
}

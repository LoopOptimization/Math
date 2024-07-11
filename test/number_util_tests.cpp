#include <cmath>
#include <gtest/gtest.h>
#ifndef USE_MODULE
#include "Math/Factor.cxx"
#else

import Factor;
#endif

TEST(FactorLowerBound, BasicAssertions) {

  for (int i = 0; i++ < 32;) {
    double di = i;
    for (int j = 0; j++ < i;) {
      double dj = j;
      auto [x, y] = math::lower_bound_factor(di, dj);
      EXPECT_EQ(di, x * y);
      EXPECT_LE(x, dj);
      EXPECT_EQ(std::round(x), x);
      EXPECT_EQ(std::round(y), y);
    }
  }
}
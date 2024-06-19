
#include "Math/Factor.hpp"
#include <cmath>
#include <gtest/gtest.h>

TEST(FactorLowerBound, BasicAssertions) {

  for (int i = 0; i++ < 16;) {
    double di = i;
    for (int j = 0; j++ < i;) {
      double dj = j;
      auto [x, y] = poly::math::lower_bound_factor(di, dj);
      EXPECT_EQ(di, x * y);
      EXPECT_LE(x, dj);
      EXPECT_EQ(std::round(x), x);
      EXPECT_EQ(std::round(y), y);
    }
  }
}

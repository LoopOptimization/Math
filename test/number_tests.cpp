#include "Numbers/Int8.hpp"
#include <cstdint>
#include <gtest/gtest.h>
#include <limits>

using poly::numbers::Int8, poly::numbers::UInt8, poly::numbers::Flag8;

TEST(Int8Test, BasicAssertions) {
  for (uint8_t x = 0; x < std::numeric_limits<uint8_t>::max(); ++x) {
    for (uint8_t y = 0; y < std::numeric_limits<uint8_t>::max(); ++y) {
      UInt8 ux = static_cast<UInt8>(x);
      UInt8 uy = static_cast<UInt8>(y);
      EXPECT_EQ(ux <=> uy, x <=> y);
      EXPECT_TRUE(ux + uy == static_cast<UInt8>(x + y));
      EXPECT_FALSE(ux + uy != static_cast<UInt8>(x + y));
      bool b = (ux > uy) == (x > y);
      EXPECT_TRUE(b);
      EXPECT_TRUE((ux == uy) == (x == y));
      EXPECT_TRUE((ux == y) == (x == y));
      EXPECT_TRUE((x == uy) == (x == y));
      EXPECT_TRUE((ux != uy) == (x != y));
      EXPECT_TRUE((ux != y) == (x != y));
      EXPECT_TRUE((x != uy) == (x != y));
      EXPECT_TRUE((ux > uy) == (x > y));
      EXPECT_TRUE((ux > y) == (x > y));
      EXPECT_TRUE((x > uy) == (x > y));
      EXPECT_TRUE((ux < uy) == (x < y));
      EXPECT_TRUE((ux < y) == (x < y));
      EXPECT_TRUE((x < uy) == (x < y));
      EXPECT_TRUE((ux >= uy) == (x >= y));
      EXPECT_TRUE((ux >= y) == (x >= y));
      EXPECT_TRUE((x >= uy) == (x >= y));
      EXPECT_TRUE((ux <= uy) == (x <= y));
      EXPECT_TRUE((ux <= y) == (x <= y));
      EXPECT_TRUE((x <= uy) == (x <= y));

      EXPECT_EQ(ux > uy, x > y);

      EXPECT_EQ(ux + uy, static_cast<UInt8>(x + y));
      EXPECT_EQ(ux * uy, static_cast<UInt8>(x * y));
      EXPECT_EQ(ux - uy, static_cast<UInt8>(x - y));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
      if (y) EXPECT_EQ(ux / uy, static_cast<UInt8>(x / y));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    }
  }
  for (int8_t x = std::numeric_limits<int8_t>::min();
       x < std::numeric_limits<int8_t>::max(); ++x) {
    for (int8_t y = std::numeric_limits<int8_t>::min();
         y < std::numeric_limits<int8_t>::max(); ++y) {
      Int8 ux = static_cast<Int8>(x);
      Int8 uy = static_cast<Int8>(y);

      EXPECT_EQ(ux + uy, static_cast<Int8>(x + y));
      EXPECT_EQ(ux * uy, static_cast<Int8>(x * y));
      EXPECT_EQ(ux - uy, static_cast<Int8>(x - y));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
      if (y) EXPECT_EQ(ux / uy, static_cast<Int8>(x / y));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    }
  }
}


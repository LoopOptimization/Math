#include <gtest/gtest.h>
import std;
import Tuple;

using containers::Tuple, containers::tie, containers::Add, containers::Pair;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(TupleTest, BasicAssertions) {
  double x = 2.3;
  double y = 4.5;
  long z = -5;
  unsigned long w = 15;
  Tuple t{._0 = w, ._1 = x, ._2 = y, ._3 = z};
  {
    auto [a, b, c, d] = t;
    static_assert(std::same_as<decltype(a), unsigned long>);
    static_assert(std::same_as<decltype(x), double>);
    static_assert(std::same_as<decltype(y), double>);
    static_assert(std::same_as<decltype(z), long>);
    EXPECT_EQ(a, w);
    EXPECT_EQ(b, x);
    EXPECT_EQ(c, y);
    EXPECT_EQ(d, z);
  }
  Tuple t1 = t.map([](auto arg) { return 3 * arg; });
  {
    auto [a, b, c, d] = t1;
    static_assert(std::same_as<decltype(a), unsigned long>);
    static_assert(std::same_as<decltype(x), double>);
    static_assert(std::same_as<decltype(y), double>);
    static_assert(std::same_as<decltype(z), long>);
    EXPECT_EQ(a, 3 * w);
    EXPECT_EQ(b, 3 * x);
    EXPECT_EQ(c, 3 * y);
    EXPECT_EQ(d, 3 * z);
  }
  t1.apply([](auto &arg) { arg *= 2; });
  {
    auto [a, b, c, d] = t1;
    static_assert(std::same_as<decltype(a), unsigned long>);
    static_assert(std::same_as<decltype(x), double>);
    static_assert(std::same_as<decltype(y), double>);
    static_assert(std::same_as<decltype(z), long>);
    EXPECT_EQ(a, 6 * w);
    EXPECT_EQ(b, 6 * x);
    EXPECT_EQ(c, 6 * y);
    EXPECT_EQ(d, 6 * z);
    tie(w, x, y, z) = t1;
    EXPECT_EQ(a, w);
    EXPECT_EQ(b, x);
    EXPECT_EQ(c, y);
    EXPECT_EQ(d, z);
  }

  // auto tr = tie(Add(x), y);
  double xtest = x + 4.7;
  tie(Add(x), y) = Tuple(4.7, 5.0);
  EXPECT_EQ(x, xtest);
  EXPECT_EQ(y, 5.0);
  tie(Add(x), y) = Pair(4.7, 8.0);
  EXPECT_EQ(x, xtest + 4.7);
  EXPECT_EQ(y, 8.0);
  auto t16 = Tuple(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
  EXPECT_EQ(t16.template get<0>(), 0);
  EXPECT_EQ(t16.template get<1>(), 1);
  EXPECT_EQ(t16.template get<2>(), 2);
  EXPECT_EQ(t16.template get<3>(), 3);
  EXPECT_EQ(t16.template get<4>(), 4);
  EXPECT_EQ(t16.template get<5>(), 5);
  EXPECT_EQ(t16.template get<6>(), 6);
  EXPECT_EQ(t16.template get<7>(), 7);
  EXPECT_EQ(t16.template get<8>(), 8);
  EXPECT_EQ(t16.template get<9>(), 9);
  EXPECT_EQ(t16.template get<10>(), 10);
  EXPECT_EQ(t16.template get<11>(), 11);
  EXPECT_EQ(t16.template get<12>(), 12);
  EXPECT_EQ(t16.template get<13>(), 13);
  EXPECT_EQ(t16.template get<14>(), 14);
  EXPECT_EQ(t16.template get<15>(), 15);
}

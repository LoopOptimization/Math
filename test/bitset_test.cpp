#include <gtest/gtest.h>

#ifndef USE_MODULE
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <print>

#include "Containers/BitSets.cxx"
#include "Math/ManagedArray.cxx"
#else

import BitSet;
import ManagedArray;
import std;
#endif

#define STRINGIZE_DETAIL(x) #x
#define STRINGIZE(x) STRINGIZE_DETAIL(x)
using containers::BitSet, math::Vector;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(BitSet, BasicAssertions) {
  BitSet bs(1000);
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  ASSERT_EQ(count, 0);
  bs[4] = true;
  bs[10] = true;
  bs[200] = true;
  bs[117] = true;
  bs[87] = true;
  bs[991] = true;
  ASSERT_EQ(bs.findFirstZero(), 0);
  bs[0] = true;
  ASSERT_EQ(bs.findFirstZero(), 1);
  bs.print();
  utils::print('\n');
  // EXPECT_EQ(std::ranges::begin(bs), bs.begin());
  // EXPECT_EQ(std::ranges::end(bs), bs.end());
  Vector<std::size_t> bsc{std::array{0, 4, 10, 87, 117, 200, 991}};
  std::ptrdiff_t j = 0;
  for (auto J = bs.begin(); J != decltype(bs)::end(); ++J) {
    EXPECT_EQ(*J, bsc[j++]);
    EXPECT_TRUE(bs[*J]);
    std::println("We get: {}", *J);
  }
  j = 0;
  for (auto i : bs) {
    EXPECT_EQ(i, bsc[j++]);
    EXPECT_TRUE(bs[i]);
    std::println("We get: {}", i);
  }
  EXPECT_EQ(j, bsc.size());
  EXPECT_EQ(j, bs.size());
  std::println("About to create empty!");
  BitSet empty;
  std::ptrdiff_t c = 0, d = 0;
  std::println("About to iterate empty!");
  for (auto b : empty) {
    ++c;
    d += b;
  }
  std::println("Iterated empty!");
  EXPECT_FALSE(c);
  EXPECT_FALSE(d);
  std::println("we made it to the end!");
}
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(BitSet, Insert) {
  BitSet<std::array<std::uint64_t, 2>> bs;
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  ASSERT_EQ(count, 0);
  bs.insert(1);
  bs.insert(5);
  bs.insert(6);
  bs.insert(8);
  EXPECT_EQ(bs.data_[0], 354);
  EXPECT_EQ(bs.data_[1], 0);
  bs.insert(5);
  EXPECT_EQ(bs.data_[0], 354);
}
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(BitSet, DynSize) {
  BitSet bs, bsd{BitSet<>::dense(11)};
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  ASSERT_EQ(count, 0);
  for (std::ptrdiff_t _ : bsd) ++count;
  ASSERT_EQ(count, 11);
  EXPECT_EQ(bs.data_.size(), 0);
  bs[4] = true;
  bs[10] = true;
  EXPECT_EQ(bs.data_.size(), 1);
  EXPECT_EQ(bs.data_.front(), 1040);
  for (std::ptrdiff_t i = 0; i < 11; ++i)
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
    if (!bs.contains(i)) EXPECT_TRUE(bsd.remove(i));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  EXPECT_EQ(bs, bsd);
  Vector<std::size_t> sv;
  for (auto i : bs) sv.push_back(i);
  EXPECT_EQ(sv.size(), 2);
  EXPECT_EQ(sv[0], 4);
  EXPECT_EQ(sv[1], 10);
}
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(BitSet, FixedSize) {
  std::println(STRINGIZE(__LINE__));
  BitSet<std::array<std::uint64_t, 2>> bs;
  std::println(STRINGIZE(__LINE__));
  bs[4] = true;
  std::println(STRINGIZE(__LINE__));
  bs[10] = true;
  std::println(STRINGIZE(__LINE__));
  EXPECT_EQ(bs.data_[0], 1040);
  std::println(STRINGIZE(__LINE__));
  EXPECT_EQ(bs.data_[1], 0);
  std::println(STRINGIZE(__LINE__));
  Vector<std::size_t> sv;
  std::println(STRINGIZE(__LINE__));
  for (auto i : bs) sv.push_back(i);
  std::println(STRINGIZE(__LINE__));
  EXPECT_EQ(sv.size(), 2);
  std::println(STRINGIZE(__LINE__));
  EXPECT_EQ(sv[0], 4);
  std::println(STRINGIZE(__LINE__));
  EXPECT_EQ(sv[1], 10);
  std::println(STRINGIZE(__LINE__));
}
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(BitSet, FixedSizeSmall) {
  using SB = BitSet<std::array<std::uint16_t, 1>>;
  static_assert(sizeof(SB) == 2);
  SB bs;
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  ASSERT_EQ(count, 0);
  bs[4] = true;
  bs[10] = true;
  bs[7] = true;
  bs.insert(5);
  EXPECT_EQ(bs.data_[0], 1200);
  Vector<std::size_t> sv;
  for (auto i : bs) sv.push_back(i);
  EXPECT_EQ(sv.size(), 4);
  EXPECT_EQ(sv[0], 4);
  EXPECT_EQ(sv[1], 5);
  EXPECT_EQ(sv[2], 7);
  EXPECT_EQ(sv[3], 10);
  EXPECT_EQ(SB::fromMask(1200).data_[0], 1200);
}
// NOLINTNEXTLINE()
TEST(BitSet, IterTest) {
  math::Vector<std::uint64_t, 1> data{};
  data.push_back(9223372036854775808U);
  data.push_back(1732766553700568065U);
  data.push_back(1891655728U);
  BitSet<> bs{data};
  ptrdiff_t sz = bs.size(), i = 0;
  for (std::ptrdiff_t a : bs) {
    utils::print(a, " ");
    ++i;
  }
  utils::print('\n');
  ASSERT_EQ(i, sz);
}

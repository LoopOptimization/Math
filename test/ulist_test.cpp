#include <gtest/gtest.h>
#ifndef USE_MODULE
#include "Alloc/Arena.cxx"
#include "Containers/UnrolledList.cxx"
#include <cstdint>
#else
import Arena;
import STL;
import UnrolledList;
#endif

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(UListTest, BasicAssertions) {
  alloc::OwningArena<> alloc;
  using i64 = long long;
  auto *list = alloc.create<containers::UList<i64>>();
  // for (i64 i = 0; i < 100; ++i) {
  for (i64 i = 0; i < 7; ++i) {
    list = list->push(&alloc, i);
    i64 s = i * (i + 1) / 2;
    EXPECT_EQ(list->reduce(0, [](i64 a, i64 b) { return a + b; }), s);
    i64 s2 = i * (i + 1) * (2 * i + 1) / 6;
    EXPECT_EQ(list->transform_reduce(0,
                                     [](i64 a, i64 &b) {
                                       b *= 2;
                                       return a + (b * b);
                                     }),
              s2 * 4);
    // undo the *2;
    list->forEachRev([](i64 &a) { a /= 2; });
    const auto *const_list = list;
    i64 c = 0;
    for (auto j : *const_list) c += j;
    EXPECT_EQ(c, s);
    c = 0;
    for (auto &&j : *list) c += (j += 3);
    EXPECT_EQ(c - (3 * (i + 1)), s);
    list->forEach([](i64 &a) { a -= 3; });
  }
}

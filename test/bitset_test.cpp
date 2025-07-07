import boost.ut;

import BitSet;
import CorePrint;
import ManagedArray;
import std;

#define STRINGIZE_DETAIL(x) #x
#define STRINGIZE(x) STRINGIZE_DETAIL(x)
using containers::BitSet, math::Vector;

int main() {
  using namespace boost::ut;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
"BitSet BasicAssertions"_test = [] {
  BitSet bs(1000);
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  expect(fatal(count == 0));
  bs[4] = true;
  bs[10] = true;
  bs[200] = true;
  bs[117] = true;
  bs[87] = true;
  bs[991] = true;
  expect(fatal(bs.findFirstZero() == 0));
  bs[0] = true;
  expect(fatal(bs.findFirstZero() == 1));
  bs.print();
  utils::print('\n');
  // expect(std::ranges::begin(bs), bs.begin());
  // expect(std::ranges::end(bs), bs.end());
  Vector<std::size_t> bsc{std::array{0, 4, 10, 87, 117, 200, 991}};
  std::ptrdiff_t j = 0;
  for (auto J = bs.begin(); J != decltype(bs)::end(); ++J) {
    expect(*J == bsc[j++]);
    expect(bs[*J]);
    utils::println("We get: ", *J);
  }
  j = 0;
  for (auto i : bs) {
    expect(i == bsc[j++]);
    expect(bs[i]);
    utils::println("We get: ", i);
  }
  expect(j == bsc.size());
  expect(j == bs.size());
  utils::println("About to create empty!");
  BitSet empty;
  std::ptrdiff_t c = 0, d = 0;
  utils::println("About to iterate empty!");
  for (auto b : empty) {
    ++c;
    d += b;
  }
  utils::println("Iterated empty!");
  expect(!c);
  expect(!d);
  utils::println("we made it to the end!");
};
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
"BitSet Insert"_test = [] {
  BitSet<std::array<std::uint64_t, 2>> bs;
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  expect(fatal(count == 0));
  bs.insert(1);
  bs.insert(5);
  bs.insert(6);
  bs.insert(8);
  expect(bs.data_[0] == 354);
  expect(bs.data_[1] == 0);
  bs.insert(5);
  expect(bs.data_[0] == 354);
};
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
"BitSet DynSize"_test = [] {
  BitSet bs, bsd{BitSet<>::dense(11)};
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  expect(fatal(count == 0));
  for (std::ptrdiff_t _ : bsd) ++count;
  expect(fatal(count == 11));
  expect(bs.data_.size() == 0);
  bs[4] = true;
  bs[10] = true;
  expect(bs.data_.size() == 1);
  expect(bs.data_.front() == 1040);
  for (std::ptrdiff_t i = 0; i < 11; ++i)
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
    if (!bs.contains(i)) expect(bsd.remove(i));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  expect(bs == bsd);
  Vector<std::size_t> sv;
  for (auto i : bs) sv.push_back(i);
  expect(sv.size() == 2);
  expect(sv[0] == 4);
  expect(sv[1] == 10);
};
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
"BitSet FixedSize"_test = [] {
  BitSet<std::array<std::uint64_t, 2>> bs;
  bs[4] = true;
  bs[10] = true;
  expect(bs.data_[0] == 1040);
  expect(bs.data_[1] == 0);
  Vector<std::size_t> sv;
  for (auto i : bs) sv.push_back(i);
  expect(sv.size() == 2);
  expect(sv[0] == 4);
  expect(sv[1] == 10);
};
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
"BitSet FixedSizeSmall"_test = [] {
  using SB = BitSet<std::array<std::uint16_t, 1>>;
  static_assert(sizeof(SB) == 2);
  SB bs;
  std::ptrdiff_t count = 0;
  for (std::ptrdiff_t _ : bs) ++count;
  expect(fatal(count == 0));
  bs[4] = true;
  bs[10] = true;
  bs[7] = true;
  bs.insert(5);
  expect(bs.data_[0] == 1200);
  Vector<std::size_t> sv;
  for (auto i : bs) sv.push_back(i);
  expect(sv.size() == 4);
  expect(sv[0] == 4);
  expect(sv[1] == 5);
  expect(sv[2] == 7);
  expect(sv[3] == 10);
  expect(SB::fromMask(1200).data_[0] == 1200);
};
// NOLINTNEXTLINE()
"BitSet IterTest"_test = [] {
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
  expect(fatal(i == sz));
};

  return 0;
}

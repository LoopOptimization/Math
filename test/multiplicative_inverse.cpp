import boost.ut;

import MultiplicativeInverse;
import std;

using namespace boost::ut;
using ::math::MultiplicativeInverse, ::math::cld;

void testBasicAssertions() {
  for (std::int32_t j = -100; j <= 100; ++j) {
    if (j == 0) continue;
    auto mij = MultiplicativeInverse(j);
    auto mijf = MultiplicativeInverse(float(j));
    for (std::int32_t i = -1000; i <= 1000; ++i) {
      auto [d, r] = mij.divrem(i);
      std::int32_t qref = i / j, rref = i % j, cref = cld(i, j);
      expect(qref == d);
      expect(rref == r);
      expect(cref == cld(i, mij));
      expect(mij * i == j * i);
      auto fi = float(i);
      auto [df, rf] = mijf.divrem(fi);
      expect(qref == df);
      expect(rref == rf);
      expect(cref == cld(float(i), mijf));
      expect(mijf * i == j * i);
    }
  }
  for (std::int64_t j = -100; j <= 100; ++j) {
    if (j == 0) continue;
    auto mij = MultiplicativeInverse(j);
    auto mijf = MultiplicativeInverse(double(j));
    for (std::int64_t i = -1000; i <= 1000; ++i) {
      auto [d, r] = mij.divrem(i);
      std::int64_t qref = i / j, rref = i % j, cref = cld(i, j);
      expect(qref == d);
      expect(rref == r);
      expect(cref == cld(i, mij));
      expect(mij * i == j * i);
      auto [df, rf] = mijf.divrem(double(i));
      expect(qref == df);
      expect(rref == rf);
      expect(cref == cld(float(i), mijf));
      expect(mijf * i == j * i);
    }
  }
  for (std::uint32_t j = 1; j <= 200; ++j) {
    auto mij = MultiplicativeInverse(j);
    for (std::uint32_t i = 0; i <= 2000; ++i) {
      auto [d, r] = mij.divrem(i);
      std::uint32_t qref = i / j, rref = i % j, cref = cld(i, j);
      expect(qref == d);
      expect(rref == r);
      expect(cref == cld(i, mij));
      expect(mij * i == j * i);
    }
  }
  for (std::uint64_t j = 1; j <= 200; ++j) {
    auto mij = MultiplicativeInverse(j);
    for (std::uint64_t i = 0; i <= 2000; ++i) {
      auto [d, r] = mij.divrem(i);
      std::uint64_t qref = i / j, rref = i % j, cref = cld(i, j);
      expect(qref == d);
      expect(rref == r);
      expect(cref == cld(i, mij));
      expect(mij * i == j * i);
    }
  }
#if __cpp_lib_constexpr_cmath >= 202202L
  static_assert(123456.0 / MultiplicativeInverse(5.0) == 123456 / 5);
#endif
  static_assert(123456 / MultiplicativeInverse(-5) == 123456 / -5);
  static_assert(unsigned(123456) / MultiplicativeInverse(unsigned(5)) ==
                123456 / 5);
}

int main() {
  "MultiplicativeInverse BasicAssertions"_test = [] { testBasicAssertions(); };
  return 0;
}

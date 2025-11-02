import boost.ut;

import BaseUtils;
import MultiplicativeInverse;
import StaticArray;
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

void svector() {
  using S = ::math::SVector<double, 8>;
  using M = MultiplicativeInverse<S>;
  int l2d0 = 3; // vectorization
  S u{2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 18.0};
  S c{6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0};
  M umi{u};
  // Currently using contig for load/store costs...
  // Without the need for shuffles, this should be >= discontig
  // TODO: double check for for `u` not being a power of 2
  S ufactor = u.max(bit::exp2unchecked(4));
  double lc = 0.5, sc = 1.0, ld = 4.5, sd = 9.0;

  S lcf = lc * ufactor, scf = sc * ufactor, shuf_count = u * l2d0,
    shuf_ratio = c / umi;
  auto prefer_shuf_over_gather = (lcf + shuf_count * lc) < ld * u,
       prefer_shuf_over_scatter = (scf + shuf_count * sc) < sd * u;
  S load_cost = select(prefer_shuf_over_gather, lcf * shuf_ratio, ld * c),
    stow_cost = select(prefer_shuf_over_scatter, scf * shuf_ratio, sd * c),
    comp_cost = select(prefer_shuf_over_gather, shuf_count * lc, 0.0);
  comp_cost << select(prefer_shuf_over_scatter, comp_cost + shuf_count * sc,
                      comp_cost);
  // comp_cost.print();
  // TODO: validate that these make sense?
  S load_ref{27, 24, 24, 24, 24, 24, 24, 9},
    stow_ref{54, 48, 48, 48, 48, 48, 48, 18},
    comp_ref{0, 13.5, 18, 22.5, 27, 31.5, 36, 81};
  expect(load_cost == load_ref);
  expect(stow_cost == stow_ref);
  expect(comp_cost == comp_ref);
}

int main() {
  "MultiplicativeInverseScalar"_test = [] { testBasicAssertions(); };
  "MultiplicativeInverseVector"_test = [] { svector(); };
  return 0;
}

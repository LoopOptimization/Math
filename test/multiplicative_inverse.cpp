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

void testSplit() {
  // Test SVector<double, 8> split
  using S8 = ::math::SVector<double, 8>;
  using S4 = ::math::SVector<double, 4>;
  using M8 = MultiplicativeInverse<S8>;
  // using M4 = MultiplicativeInverse<S4>;

  S8 v8{2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

  // Test StaticArray split
  auto halves = v8.split();
  expect(halves[0][0] == 2.0);
  expect(halves[0][1] == 3.0);
  expect(halves[0][2] == 4.0);
  expect(halves[0][3] == 5.0);
  expect(halves[1][0] == 6.0);
  expect(halves[1][1] == 7.0);
  expect(halves[1][2] == 8.0);
  expect(halves[1][3] == 9.0);

  // Test MultiplicativeInverse split
  M8 m8{v8};
  auto m_halves = m8.split();

  // Check that the divisors were split correctly
  expect(S4(m_halves[0]) == halves[0]);
  expect(S4(m_halves[1]) == halves[1]);

  // Check that division still works correctly with split inverses
  S4 test_val_0{10.0, 12.0, 16.0, 20.0};
  S4 test_val_1{24.0, 28.0, 32.0, 36.0};

  S4 result_0 = test_val_0 / m_halves[0];
  S4 result_1 = test_val_1 / m_halves[1];

  expect(result_0[0] == 10.0 / 2.0);
  expect(result_0[1] == 12.0 / 3.0);
  expect(result_0[2] == 16.0 / 4.0);
  expect(result_0[3] == 20.0 / 5.0);

  expect(result_1[0] == 24.0 / 6.0);
  expect(result_1[1] == 28.0 / 7.0);
  expect(result_1[2] == 32.0 / 8.0);
  expect(result_1[3] == 36.0 / 9.0);

  // Test that inverse is actually split, not recomputed
  // (this is more of a check that the implementation is efficient)
  auto inv_orig = m8.inv();
  auto inv_halves_0 = m_halves[0].inv();
  auto inv_halves_1 = m_halves[1].inv();

  expect(inv_halves_0[0] == inv_orig[0]);
  expect(inv_halves_0[1] == inv_orig[1]);
  expect(inv_halves_0[2] == inv_orig[2]);
  expect(inv_halves_0[3] == inv_orig[3]);
  expect(inv_halves_1[0] == inv_orig[4]);
  expect(inv_halves_1[1] == inv_orig[5]);
  expect(inv_halves_1[2] == inv_orig[6]);
  expect(inv_halves_1[3] == inv_orig[7]);

  // Test SVector<float, 4> split
  using F4 = ::math::SVector<float, 4>;
  using F2 = ::math::SVector<float, 2>;
  using MF4 = MultiplicativeInverse<F4>;

  F4 vf4{2.0f, 4.0f, 8.0f, 16.0f};
  auto f_halves = vf4.split();
  expect(f_halves[0][0] == 2.0f);
  expect(f_halves[0][1] == 4.0f);
  expect(f_halves[1][0] == 8.0f);
  expect(f_halves[1][1] == 16.0f);

  MF4 mf4{vf4};
  auto mf_halves = mf4.split();

  F2 test_f_0{10.0f, 20.0f};
  F2 test_f_1{40.0f, 80.0f};

  F2 result_f_0 = test_f_0 / mf_halves[0];
  F2 result_f_1 = test_f_1 / mf_halves[1];

  expect(result_f_0[0] == 10.0f / 2.0f);
  expect(result_f_0[1] == 20.0f / 4.0f);
  expect(result_f_1[0] == 40.0f / 8.0f);
  expect(result_f_1[1] == 80.0f / 16.0f);
}

int main() {
  "MultiplicativeInverseScalar"_test = [] { testBasicAssertions(); };
  "MultiplicativeInverseVector"_test = [] { svector(); };
  "MultiplicativeInverseSplit"_test = [] { testSplit(); };
  return 0;
}

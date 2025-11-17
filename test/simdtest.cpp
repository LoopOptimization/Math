import boost.ut;
import AxisTypes;
import Elementary;
import ExprTemplates;
import ManagedArray;
import SIMD;
import std;

using namespace boost::ut;

namespace {
void testBasicAssertions() {

  ::math::Vector<double> x{
    std::array{1.564299025169599, 4.328641127555183, -10.43843599926044,
               -1.650625233314754, -0.5851694806951444, 0.07422197149516746,
               -5.231238164802142, 7.298495240920298, 6.983762398719033,
               4.148462859390057, -1.5877584351154996, 0.21815818734467177,
               6.717052006977858, 1.1366111295246064, -3.2809440656442757}},
    y{std::array{4.779323575889951, 75.8411580559438, 2.9284974204673125e-5,
                 0.19192987014972496, 0.5570114511240655, 1.077045852442665,
                 0.005346900857986191, 1478.074107909007, 1078.9702561142797,
                 63.336568222798185, 0.20438323693119723, 1.2437838029152637,
                 826.3777611913293, 3.116190086521906, 0.037592750025674}},
    z{::math::length(15)};

  z << ::math::elementwise(x.view(), [](auto a) { return ::math::exp(a); });
  for (std::ptrdiff_t i = 0; i < 15; ++i) expect(approx(y[i], z[i], 1e-14));
}

void testMaskCat() {
#ifdef __AVX512VL__
  // Test mask::Bit cat
  using Bit2 = simd::mask::Bit<2>;
  using Bit4 = simd::mask::Bit<4>;

  Bit2 m1{0b11}; // both bits set
  Bit2 m2{0b01}; // only first bit set

  Bit4 m_cat = m1.cat(m2);
  expect(eq(m_cat.intmask(), 0b0111)); // m1 in lower 2 bits, m2 in upper 2 bits

  Bit2 m3{0b10}; // only second bit set
  Bit2 m4{0b00}; // no bits set

  Bit4 m_cat2 = m3.cat(m4);
  expect(
    eq(m_cat2.intmask(), 0b0010)); // m3 in lower 2 bits, m4 in upper 2 bits
#else
  // Test mask::Vector cat
  using Vec2 = simd::mask::Vector<2, sizeof(int)>;
  using Vec4 = simd::mask::Vector<4, sizeof(int)>;

  Vec2 m1{{-1, -1}}; // both elements masked
  Vec2 m2{{-1, 0}};  // first element masked

  Vec4 m_cat = m1.cat(m2);
  expect(eq(m_cat.popcount(), 3_i)); // 3 elements should be masked

  Vec2 m3{{0, -1}}; // second element masked
  Vec2 m4{{0, 0}};  // no elements masked

  Vec4 m_cat2 = m3.cat(m4);
  expect(eq(m_cat2.popcount(), 1_i)); // 1 element should be masked
#endif
}

template <std::ptrdiff_t N> void testMaskUnrollCast() {
  // Test mask::Unroll cast from Unroll<1,8,2> to Unroll<1,4,4>
  // Create a mask pattern: alternating true/false
#ifdef __AVX512VL__
  simd::mask::Unroll<1, 2 * N, 2> mask_8_2;
  for (std::ptrdiff_t i = 0; i < 2 * N; ++i) {
    mask_8_2[i] =
      simd::mask::Bit<2>{static_cast<std::uint64_t>(i % 2 ? 0b11 : 0b00)};
  }

  // Cast to Unroll<1,4,4>
  simd::mask::Unroll<1, N, 4> mask_4_4 = mask_8_2;

  // Verify the pattern is preserved
  // mask_8_2: [00, 11, 00, 11, 00, 11, 00, 11]
  // mask_4_4 after cat: [0011, 0011, 0011, 0011]
  for (std::ptrdiff_t i = 0; i < 4; ++i)
    expect(eq(mask_4_4[i].intmask(), 0b1100));
#else
  static constexpr std::size_t OldBytes = 8;
  static constexpr std::size_t NewBytes = 4;
  using Vec2 = simd::mask::Vector<2, OldBytes>;
  simd::mask::Unroll<1, 2 * N, 2, OldBytes> mask_8_2;
  for (std::ptrdiff_t i = 0; i < 2 * N; ++i)
    if (i % 2) mask_8_2[i] = Vec2{{-1, -1}};
    else mask_8_2[i] = Vec2{{0, 0}};

  // Cast to Unroll<1,4,4>
  simd::mask::Unroll<1, N, 4, NewBytes> mask_4_4 = mask_8_2;

  // Verify the pattern is preserved
  for (std::ptrdiff_t i = 0; i < 4; ++i)
    expect(eq(mask_4_4[i].popcount(), 2_i));
#endif
}

void testMaskSelectWithCast() {
  // Test select with mismatched mask and data dimensions
  // mask::Unroll<1,8,2> with simd::Unroll<1,4,4,int>
#ifdef __AVX512VL__
  simd::mask::Unroll<1, 8, 2> mask;
  // Create a mask where even indices are all-false, odd are all-true
  for (std::ptrdiff_t i = 0; i < 8; ++i) {
    mask[i] =
      simd::mask::Bit<2>{static_cast<std::uint64_t>(i % 2 ? 0b11 : 0b00)};
  }

  simd::Unroll<1, 4, 4, int> a, b;
  for (std::ptrdiff_t i = 0; i < 4; ++i) {
    a[i] = simd::Vec<4, int>{100, 101, 102, 103};
    b[i] = simd::Vec<4, int>{200, 201, 202, 203};
  }

  auto result = mask.select(a, b);

  // Verify select worked correctly
  // After cast: mask_4_4[i] = cat(mask_8_2[2*i], mask_8_2[2*i+1])
  //           = cat(0b00, 0b11) = 0b1100 for all i
  // So result should select: b[0:1], a[2:3] for each vector
  for (std::ptrdiff_t i = 0; i < 4; ++i) {
    expect(eq(result[i][0], 200_i)); // false -> b
    expect(eq(result[i][1], 201_i)); // false -> b
    expect(eq(result[i][2], 102_i)); // true -> a
    expect(eq(result[i][3], 103_i)); // true -> a
  }
#else
  using Vec2 = simd::mask::Vector<2, sizeof(int)>;
  simd::mask::Unroll<1, 8, 2, sizeof(int)> mask;
  for (std::ptrdiff_t i = 0; i < 8; ++i)
    if (i % 2) mask[i] = Vec2{{-1, -1}};
    else mask[i] = Vec2{{0, 0}};

  simd::Unroll<1, 4, 4, int> a, b;
  for (std::ptrdiff_t i = 0; i < 4; ++i) {
    a[i] = simd::Vec<4, int>{100, 101, 102, 103};
    b[i] = simd::Vec<4, int>{200, 201, 202, 203};
  }

  auto result = mask.select(a, b);

  // Verify select worked correctly
  for (std::ptrdiff_t i = 0; i < 4; ++i) {
    expect(eq(result[i][0], 200_i)); // false -> b
    expect(eq(result[i][1], 201_i)); // false -> b
    expect(eq(result[i][2], 102_i)); // true -> a
    expect(eq(result[i][3], 103_i)); // true -> a
  }
#endif
}
} // namespace

auto main() -> int {
  "ElementarySIMD BasicAssertions"_test = [] -> void { testBasicAssertions(); };
  "MaskCat"_test = [] -> void { testMaskCat(); };
  "MaskUnrollCast"_test = [] -> void {
    testMaskUnrollCast<4>();
    testMaskUnrollCast<1>();
  };
  "MaskSelectWithCast"_test = [] -> void { testMaskSelectWithCast(); };
  return 0;
}

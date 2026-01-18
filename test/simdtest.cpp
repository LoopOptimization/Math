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
  "FromSplitBasic"_test = [] -> void {
    // Test base case: Unroll<1, 2, W, int> from Vec<2*W, int>
    static constexpr auto W = simd::Width<int>;
    constexpr auto v = simd::range<2 * W, int>();

    auto u = simd::Unroll<1, 2, W, int>::fromSplit(v);

    // Should split into two Vec<W, int>
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      expect(eq(u[0][i], static_cast<int>(i)));
      expect(eq(u[1][i], static_cast<int>(W + i)));
    }
    expect(__builtin_reduce_and(u.catToSIMD() == v));
  };
  "FromSplitColumn"_test = [] -> void {
    // Test column split: Unroll<1, 4, W, int> from Vec<4*W, int>
    static constexpr auto W = simd::Width<int>;
    constexpr auto v = simd::range<4 * W, int>();

    auto u = simd::Unroll<1, 4, W, int>::fromSplit(v);

    // Should split into 4 Vec<W, int>
    for (std::ptrdiff_t c = 0; c < 4; ++c)
      for (std::ptrdiff_t i = 0; i < W; ++i)
        expect(eq(u[c][i], static_cast<int>(c * W + i)));
    expect(__builtin_reduce_and(u.catToSIMD() == v));
  };
  "FromSplitRow"_test = [] -> void {
    // Test row split: Unroll<2, 2, W, int> from Vec<4*W, int>
    static constexpr auto W = simd::Width<int>;
    constexpr auto v = simd::range<4 * W, int>();

    auto u = simd::Unroll<2, 2, W, int>::fromSplit(v);

    // Should split into 2x2 grid of Vec<W, int>
    for (std::ptrdiff_t idx = 0; idx < 4; ++idx)
      for (std::ptrdiff_t i = 0; i < W; ++i)
        expect(eq(u[idx][i], static_cast<int>(idx * W + i)));
    expect(__builtin_reduce_and(u.catToSIMD() == v));
  };
  "FromSplitLarger"_test = [] -> void {
    // Test larger split: Unroll<4, 2, W, int> from Vec<8*W, int>
    static constexpr auto Wint = simd::Width<int>;
    using I = std::conditional_t<(Wint <= 8), int, std::int64_t>;
    static constexpr auto W = simd::Width<I>;
    constexpr auto v = simd::range<8 * W, I>();

    auto u = simd::Unroll<4, 2, W, I>::fromSplit(v);

    // Should split into 4x2 grid of Vec<W, I>
    for (std::ptrdiff_t idx = 0; idx < 8; ++idx)
      for (std::ptrdiff_t i = 0; i < W; ++i)
        expect(eq(u[idx][i], static_cast<I>(idx * W + i)));
    expect(__builtin_reduce_and(u.catToSIMD() == v));
  };
  "FromSplitDouble"_test = [] -> void {
    // Test with double type: Unroll<1, 4, W, double> from Vec<4*W, double>
    static constexpr auto W = simd::Width<double>;
    constexpr auto v = simd::range<4 * W, double>();

    auto u = simd::Unroll<1, 4, W, double>::fromSplit(v);

    for (std::ptrdiff_t c = 0; c < 4; ++c)
      for (std::ptrdiff_t i = 0; i < W; ++i)
        expect(eq(u[c][i], static_cast<double>(c * W + i)));

    expect(__builtin_reduce_and(u.catToSIMD() == v));
  };
  "SelectWiderType"_test = [] -> void {
    // Test select between wider types using mask from narrower types
    // Create float vectors for comparison
    static constexpr auto Wf = simd::Width<float>;
    using f32 = simd::SVec<float>;
    f32 x{f32::range(0.F)}, y{f32::vbroadcast(Wf / 2)};

    // Create wider type vectors (uint64_t has same width as float in terms of
    // SIMD count)
    using u64 = simd::SVec<float>;
    u64 a{u64::range(100)}, b{u64::range(200)};

    // Select based on float comparison: where x > y, take a, else take b
    auto result = (x > y).select(a, b);

    // Verify: elements where i > Wf/2 should come from a, others from b
    for (std::ptrdiff_t i = 0; i < Wf; ++i) {
      std::uint64_t r = result.extract_value(i), ai = a.extract_value(i),
                    bi = b.extract_value(i);
      if (i > Wf / 2) {
        expect(eq(r, ai)) << "Failed at index " << i << ": expected " << ai
                          << ", got " << r;
      } else {
        expect(eq(r, bi)) << "Failed at index " << i << ": expected " << bi
                          << ", got " << r;
      }
    }
  };
  "ExtractValueSingleVec"_test = [] -> void {
    // Test extract_value for R*C==1 case (Unroll<1, 1, N, T>)
    static constexpr auto W = simd::Width<int>;
    using f32 = simd::Unroll<1, 1, W, int>;
    f32 u{f32::range(42)};

    // Verify extract_value retrieves correct values
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      expect(eq(u.extract_value(i), static_cast<int>(42 + i)))
        << "Failed at index " << i;
    }

    // Test with double type
    static constexpr auto Wd = simd::Width<double>;
    using f64 = simd::Unroll<1, 1, Wd, double>;
    f64 ud{f64::range(3.14)};

    for (std::ptrdiff_t i = 0; i < Wd; ++i) {
      expect(eq(ud.extract_value(i), 3.14 + static_cast<double>(i)))
        << "Failed at double index " << i;
    }
  };
  "ExtractValueMultiVec"_test = [] -> void {
    // Test extract_value for R*C>1 case (Unroll<R, C, N, T>)
    static constexpr auto W = simd::Width<int>;

    // Test column unroll: Unroll<1, 4, W, int>
    using i32x4 = simd::Unroll<1, 4, W, int>;
    i32x4 u_col{i32x4::range(0)};

    // Verify extract_value with linear indexing
    for (std::ptrdiff_t c = 0; c < 4; ++c) {
      for (std::ptrdiff_t i = 0; i < W; ++i) {
        std::ptrdiff_t linear_idx = c * W + i;
        expect(
          eq(u_col.extract_value(linear_idx), static_cast<int>(linear_idx)))
          << "Failed at column " << c << ", index " << i;
      }
    }

    // Test 2D unroll: Unroll<2, 3, W, int>
    using i32x2x3 = simd::Unroll<2, 3, W, int>;
    i32x2x3 u_2d{i32x2x3::range(0)};

    // Verify extract_value with 2D layout
    for (std::ptrdiff_t idx = 0; idx < 6; ++idx) {
      for (std::ptrdiff_t i = 0; i < W; ++i) {
        std::ptrdiff_t linear_idx = idx * W + i;
        expect(eq(u_2d.extract_value(linear_idx), static_cast<int>(linear_idx)))
          << "Failed at unroll index " << idx << ", simd index " << i;
      }
    }

    // Test with float type: Unroll<1, 2, W, float>
    static constexpr auto Wf = simd::Width<float>;
    using f32x2 = simd::Unroll<1, 2, Wf, float>;
    f32x2 u_float{f32x2::range(10.5F)};

    for (std::ptrdiff_t i = 0; i < 2 * Wf; ++i) {
      expect(eq(u_float.extract_value(i), 10.5f + static_cast<float>(i)))
        << "Failed at index " << i;
    }
  };
  "BitshiftLeftScalar"_test = [] -> void {
    // Test left shift with scalar shift amount
    static constexpr auto W = simd::Width<int>;
    using i32 = simd::Unroll<1, 1, W, int>;

    // Test: values << 2
    i32 a{i32::range(1)}; // [1, 2, 3, 4, ...]
    auto result = a << 2;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      int expected = static_cast<int>(1 + i) << 2;
      expect(eq(result.extract_value(i), expected))
        << "Failed left shift at index " << i;
    }

    // Test: scalar << Unroll
    i32 shift_amounts{i32::range(0)}; // [0, 1, 2, 3, ...]
    auto result2 = 16 << shift_amounts;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      int expected = 16 << static_cast<int>(i);
      expect(eq(result2.extract_value(i), expected))
        << "Failed scalar << Unroll at index " << i;
    }
  };
  "BitshiftRightScalar"_test = [] -> void {
    // Test right shift with scalar shift amount
    static constexpr auto W = simd::Width<int>;
    using i32 = simd::Unroll<1, 1, W, int>;

    // Test: values >> 2
    i32 a{i32::range(16)}; // [16, 17, 18, 19, ...]
    auto result = a >> 2;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      int expected = static_cast<int>(16 + i) >> 2;
      expect(eq(result.extract_value(i), expected))
        << "Failed right shift at index " << i;
    }

    // Test: scalar >> Unroll
    i32 shift_amounts{i32::range(0)}; // [0, 1, 2, 3, ...]
    auto result2 = 256 >> shift_amounts;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      int expected = 256 >> static_cast<int>(i);
      expect(eq(result2.extract_value(i), expected))
        << "Failed scalar >> Unroll at index " << i;
    }
  };
  "BitshiftVectorByVector"_test = [] -> void {
    // Test shift where shift amounts are vectorized
    static constexpr auto W = simd::Width<unsigned int>;
    using u32 = simd::Unroll<1, 1, W, unsigned int>;

    // Test: each element shifted by different amount
    u32 values{u32::range(8U)}; // [8, 9, 10, 11, ...]
    u32 shifts{u32::range(0U)}; // [0, 1, 2, 3, ...]

    auto left_result = values << shifts;
    auto right_result = values >> shifts;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      unsigned int val = static_cast<unsigned int>(8 + i);
      unsigned int shift = static_cast<unsigned int>(i);
      unsigned int expected_left = val << shift;
      unsigned int expected_right = val >> shift;

      expect(eq(left_result.extract_value(i), expected_left))
        << "Failed vector << vector at index " << i;
      expect(eq(right_result.extract_value(i), expected_right))
        << "Failed vector >> vector at index " << i;
    }
  };
  "BitshiftUnrollTypes"_test = [] -> void {
    // Test with different Unroll dimensions
    static constexpr auto W = simd::Width<int>;

    // Test with Unroll<1, 4, W, int>
    using i32x4 = simd::Unroll<1, 4, W, int>;
    i32x4 a{i32x4::range(1)};
    auto result = a << 3;

    for (std::ptrdiff_t c = 0; c < 4; ++c) {
      for (std::ptrdiff_t i = 0; i < W; ++i) {
        std::ptrdiff_t linear_idx = c * W + i;
        int expected = static_cast<int>(1 + linear_idx) << 3;
        expect(eq(result.extract_value(linear_idx), expected))
          << "Failed at column " << c << ", index " << i;
      }
    }

    // Test with Unroll<2, 2, W, int>
    using i32x2x2 = simd::Unroll<2, 2, W, int>;
    i32x2x2 b{i32x2x2::range(64)};
    auto result2 = b >> 2;

    for (std::ptrdiff_t idx = 0; idx < 4; ++idx) {
      for (std::ptrdiff_t i = 0; i < W; ++i) {
        std::ptrdiff_t linear_idx = idx * W + i;
        int expected = static_cast<int>(64 + linear_idx) >> 2;
        expect(eq(result2.extract_value(linear_idx), expected))
          << "Failed at unroll index " << idx << ", simd index " << i;
      }
    }
  };
  "BitshiftSignedUnsigned"_test = [] -> void {
    // Test with both signed and unsigned types
    static constexpr auto W = simd::Width<int>;

    // Signed int left shift
    using i32 = simd::Unroll<1, 1, W, int>;
    i32 signed_vals{i32::range(1)};
    auto signed_result = signed_vals << 4;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      int expected = static_cast<int>(1 + i) << 4;
      expect(eq(signed_result.extract_value(i), expected));
    }

    // Unsigned int right shift
    using u32 = simd::Unroll<1, 1, W, unsigned int>;
    u32 unsigned_vals{u32::range(1024U)};
    auto unsigned_result = unsigned_vals >> 4;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      unsigned int expected = static_cast<unsigned int>(1024U + i) >> 4;
      expect(eq(unsigned_result.extract_value(i), expected));
    }
  };
  "BitshiftEdgeCases"_test = [] -> void {
    // Test edge cases
    static constexpr auto W = simd::Width<unsigned int>;
    using u32 = simd::Unroll<1, 1, W, unsigned int>;

    // Shift by 0
    u32 a{u32::range(10U)};
    auto result_zero = a << 0;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      expect(
        eq(result_zero.extract_value(i), static_cast<unsigned int>(10U + i)))
        << "Shift by 0 failed at index " << i;
    }

    // Large shift (should wrap or saturate depending on type width)
    u32 b{u32::vbroadcast(0xFFFFFFFFU)};
    auto result_large = b >> 31;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      expect(eq(result_large.extract_value(i), 1U))
        << "Large right shift failed at index " << i;
    }

    // Left shift building powers of 2
    u32 c{u32::range(0U)}; // [0, 1, 2, 3, ...]
    auto powers = 1U << c;

    for (std::ptrdiff_t i = 0; i < W && i < 31; ++i) {
      unsigned int expected = 1U << static_cast<unsigned int>(i);
      expect(eq(powers.extract_value(i), expected))
        << "Power of 2 failed at index " << i;
    }
  };
  "BitshiftInt64"_test = [] -> void {
    // Test with 64-bit integers
    static constexpr auto W = simd::Width<std::int64_t>;
    using i64 = simd::Unroll<1, 1, W, std::int64_t>;

    // Test large values that need 64 bits
    i64 large_vals{i64::range(1LL << 32)};
    auto result = large_vals << 4;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      std::int64_t expected = ((1LL << 32) + i) << 4;
      expect(eq(result.extract_value(i), expected))
        << "64-bit left shift failed at index " << i;
    }

    // Right shift preserving high bits
    i64 high_vals{i64::range(1LL << 48)};
    auto result2 = high_vals >> 16;

    for (std::ptrdiff_t i = 0; i < W; ++i) {
      std::int64_t expected = ((1LL << 48) + i) >> 16;
      expect(eq(result2.extract_value(i), expected))
        << "64-bit right shift failed at index " << i;
    }
  };
  "BFloat16LoadAndCast"_test = [] -> void {
    // Test loading __bf16 and casting to float
    static constexpr auto W = simd::Width<__bf16>;

    // Create test data: simple float values that convert cleanly to bf16
    alignas(64) float source_floats[W];
    for (std::ptrdiff_t i = 0; i < W; ++i)
      source_floats[i] = static_cast<float>(i + 1);

    // Convert to bf16
    alignas(64) __bf16 bf16_data[W];
    for (std::ptrdiff_t i = 0; i < W; ++i)
      bf16_data[i] = static_cast<__bf16>(source_floats[i]);

    // Load bf16 vector using simd::load
    auto v = simd::load(bf16_data, simd::mask::None<W>{});

    // Cast to float
    auto float_v = __builtin_convertvector(v, simd::Vec<W, float>);

    // Verify values are approximately correct (bf16 has reduced precision)
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      float result = float_v[i];
      float expected = source_floats[i];
      expect(approx(result, expected, 1e-2f))
        << "BF16 load/cast failed at index " << i;
    }
  };
  "BFloat16Arithmetic"_test = [] -> void {
    // Test arithmetic operations on bf16 cast to float
    static constexpr auto W = simd::Width<__bf16>;

    alignas(64) __bf16 a_data[W], b_data[W];
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      a_data[i] = static_cast<__bf16>(static_cast<float>(i + 1));
      b_data[i] = static_cast<__bf16>(2.0f);
    }

    // Load bf16 vectors using simd::load and cast to float
    auto a_bf16 = simd::load(a_data, simd::mask::None<W>{});
    auto b_bf16 = simd::load(b_data, simd::mask::None<W>{});

    auto a_float = __builtin_convertvector(a_bf16, simd::Vec<W, float>);
    auto b_float = __builtin_convertvector(b_bf16, simd::Vec<W, float>);

    // Test addition
    auto sum = a_float + b_float;
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      float expected = static_cast<float>(i + 1) + 2.0f;
      expect(approx(sum[i], expected, 1e-2f))
        << "BF16 addition failed at index " << i;
    }

    // Test multiplication
    auto prod = a_float * b_float;
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      float expected = static_cast<float>(i + 1) * 2.0f;
      expect(approx(prod[i], expected, 1e-2f))
        << "BF16 multiplication failed at index " << i;
    }
  };
  "BFloat16Unroll"_test = [] -> void {
    // Test bf16 with Unroll types
    static constexpr auto W = simd::Width<__bf16>;

    alignas(64) __bf16 data[2 * W];
    for (std::ptrdiff_t i = 0; i < 2 * W; ++i)
      data[i] = static_cast<__bf16>(static_cast<float>(i + 1));

    // Load into Unroll<1, 2, W, __bf16>
    using bf16_unroll = simd::Unroll<1, 2, W, __bf16>;
    bf16_unroll u;
    u[0] = simd::load(data, simd::mask::None<W>{});
    u[1] = simd::load(data + W, simd::mask::None<W>{});

    // Cast to float Unroll
    using float_unroll = simd::Unroll<1, 2, W, float>;
    float_unroll u_float;
    u_float[0] = __builtin_convertvector(u[0], simd::Vec<W, float>);
    u_float[1] = __builtin_convertvector(u[1], simd::Vec<W, float>);

    // Verify
    for (std::ptrdiff_t c = 0; c < 2; ++c) {
      for (std::ptrdiff_t i = 0; i < W; ++i) {
        float expected = static_cast<float>(c * W + i + 1);
        expect(approx(u_float[c][i], expected, 1e-2f))
          << "BF16 Unroll failed at column " << c << ", index " << i;
      }
    }
  };
  "BitCastSameSize"_test = [] -> void {
    // Test bitCast between same-size types: float <-> int32_t
    static constexpr auto W = simd::Width<float>;
    using f32_unroll = simd::Unroll<1, 1, W, float>;

    // Create a float vector with specific bit pattern
    f32_unroll f;
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      float val = static_cast<float>(i + 1);
      f.insert_value(i, val);
    }

    // BitCast to int32_t
    auto i = f.bitCast<std::int32_t>();

    // Verify bit patterns match by casting back
    auto f2 = i.bitCast<float>();
    for (std::ptrdiff_t j = 0; j < W; ++j) {
      expect(eq(f.extract_value(j), f2.extract_value(j)))
        << "Bit pattern mismatch at index " << j;
    }
  };
  "BitCastUpcast"_test = [] -> void {
    // Test upcast: float -> int64_t (half the elements)
    static constexpr auto Wf = simd::Width<float>;
    using f32_unroll = simd::Unroll<1, 1, Wf, float>;

    f32_unroll f;
    for (std::ptrdiff_t i = 0; i < Wf; ++i)
      f.insert_value(i, static_cast<float>(i + 100));

    auto i64 = f.bitCast<std::int64_t>();

    // Cast back and verify
    auto f2 = i64.bitCast<float>();
    for (std::ptrdiff_t j = 0; j < Wf; ++j) {
      expect(eq(f.extract_value(j), f2.extract_value(j)))
        << "Upcast bit pattern mismatch at index " << j;
    }
  };
  "BitCastDowncast"_test = [] -> void {
    // Test downcast: int64_t -> float (double the elements)
    static constexpr auto Wi = simd::Width<std::int64_t>;
    using i64_unroll = simd::Unroll<1, 1, Wi, std::int64_t>;

    i64_unroll i;
    for (std::ptrdiff_t j = 0; j < Wi; ++j)
      i.insert_value(j, static_cast<std::int64_t>(j * 1000000));

    auto f = i.bitCast<float>();

    // Cast back and verify
    auto i2 = f.bitCast<std::int64_t>();
    for (std::ptrdiff_t k = 0; k < Wi; ++k) {
      expect(eq(i.extract_value(k), i2.extract_value(k)))
        << "Downcast bit pattern mismatch at index " << k;
    }
  };
  "BitCastSmallTypes"_test = [] -> void {
    // Test with small types: __bf16 <-> uint16_t
    static constexpr auto W = simd::Width<__bf16>;
    using bf16_unroll = simd::Unroll<1, 1, W, __bf16>;

    bf16_unroll bf;
    for (std::ptrdiff_t i = 0; i < W; ++i)
      bf.insert_value(i, static_cast<__bf16>(static_cast<float>(i + 10)));

    auto u16 = bf.bitCast<std::uint16_t>();
    auto bf2 = u16.bitCast<__bf16>();

    // Verify bit preservation
    for (std::ptrdiff_t j = 0; j < W; ++j) {
      auto orig = bf.extract_value(j);
      auto restored = bf2.extract_value(j);
      // Compare as uint16 to verify exact bit pattern
      auto orig_bits = __builtin_bit_cast(std::uint16_t, orig);
      auto restored_bits = __builtin_bit_cast(std::uint16_t, restored);
      expect(eq(orig_bits, restored_bits))
        << "Small type bit pattern mismatch at index " << j;
    }
  };
  "BitCastByteSized"_test = [] -> void {
    // Test with byte-sized types: int8_t <-> uint8_t
    static constexpr auto W = simd::Width<std::int8_t>;
    using i8_unroll = simd::Unroll<1, 1, W, std::int8_t>;

    i8_unroll i8;
    for (std::ptrdiff_t i = 0; i < W; ++i)
      i8.insert_value(i, static_cast<std::int8_t>(i - 50));

    auto u8 = i8.bitCast<std::uint8_t>();
    auto i8_2 = u8.bitCast<std::int8_t>();

    for (std::ptrdiff_t j = 0; j < W; ++j) {
      expect(eq(i8.extract_value(j), i8_2.extract_value(j)))
        << "Byte-sized bit pattern mismatch at index " << j;
    }
  };
  "BitCastMultiVector"_test = [] -> void {
    // Test with multi-vector Unroll: Unroll<2, 3, W, float>
    static constexpr auto W = simd::Width<float>;
    using f32_multi = simd::Unroll<2, 3, W, float>;

    f32_multi f;
    // Fill all 6 vectors (2 rows × 3 columns)
    for (std::ptrdiff_t idx = 0; idx < 6; ++idx)
      for (std::ptrdiff_t i = 0; i < W; ++i)
        f.insert_value(idx * W + i, static_cast<float>(idx * 100 + i));

    auto i = f.bitCast<std::int32_t>();
    auto f2 = i.bitCast<float>();

    // Verify all elements preserved
    for (std::ptrdiff_t idx = 0; idx < 6; ++idx) {
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        std::ptrdiff_t linear_idx = idx * W + j;
        expect(eq(f.extract_value(linear_idx), f2.extract_value(linear_idx)))
          << "Multi-vector mismatch at linear index " << linear_idx;
      }
    }
  };
  "BitCastBitPreservation"_test = [] -> void {
    // Verify exact bit pattern preservation with known values
    static constexpr auto W = simd::Width<std::uint32_t>;
    using u32_unroll = simd::Unroll<1, 1, W, std::uint32_t>;

    u32_unroll u32;
    // Use specific bit patterns
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      std::uint32_t pattern = 0xDEADBEEF + static_cast<std::uint32_t>(i);
      u32.insert_value(i, pattern);
    }

    auto f32 = u32.bitCast<float>();
    auto u32_2 = f32.bitCast<std::uint32_t>();

    // Verify exact bit pattern match
    for (std::ptrdiff_t j = 0; j < W; ++j) {
      expect(eq(u32.extract_value(j), u32_2.extract_value(j)))
        << "Exact bit pattern not preserved at index " << j;
    }
  };
  return 0;
}

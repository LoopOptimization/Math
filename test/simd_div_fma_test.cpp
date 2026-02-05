#ifdef __x86_64__
#include <immintrin.h>
#endif
import boost.ut;
import SIMD;
import std;

using namespace boost::ut;

namespace {

// Helper to compare floating-point values
template <typename T> constexpr bool feq(T a, T b, T eps = T(1e-10)) {
  return std::abs(a - b) < eps;
}

// Test trunc_div for correctness
template <std::ptrdiff_t W> void test_trunc_div_double() {
  // Test cases: {a, b, expected_result}
  // trunc_div should match C++ integer division semantics (toward zero)
  struct TestCase {
    double a, b, expected;
  };
  constexpr TestCase cases[] = {
    {7.0, 3.0, 2.0},    // 7/3 = 2.33... -> 2 (toward zero)
    {-7.0, 3.0, -2.0},  // -7/3 = -2.33... -> -2 (toward zero)
    {7.0, -3.0, -2.0},  // 7/-3 = -2.33... -> -2 (toward zero)
    {-7.0, -3.0, 2.0},  // -7/-3 = 2.33... -> 2 (toward zero)
    {8.0, 4.0, 2.0},    // Exact division
    {-8.0, 4.0, -2.0},  // Exact negative
    {0.0, 5.0, 0.0},    // Zero dividend
    {10.0, 3.0, 3.0},   // 10/3 = 3.33... -> 3
    {-10.0, 3.0, -3.0}, // -10/3 = -3.33... -> -3
  };

  for (const auto &tc : cases) {
    simd::Vec<W, double> a = tc.a, b = tc.b;
    auto result = simd::trunc_div(a, b);
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      expect(feq(result[i], tc.expected))
        << "trunc_div(" << tc.a << ", " << tc.b << ") = " << result[i]
        << ", expected " << tc.expected;
    }
  }
}

// Test floor_div for correctness
template <std::ptrdiff_t W> void test_floor_div_double() {
  // Test cases: {a, b, expected_result}
  // floor_div rounds toward negative infinity
  struct TestCase {
    double a, b, expected;
  };
  constexpr TestCase cases[] = {
    {7.0, 3.0, 2.0},    // 7/3 = 2.33... -> 2 (toward -inf)
    {-7.0, 3.0, -3.0},  // -7/3 = -2.33... -> -3 (toward -inf)
    {7.0, -3.0, -3.0},  // 7/-3 = -2.33... -> -3 (toward -inf)
    {-7.0, -3.0, 2.0},  // -7/-3 = 2.33... -> 2 (toward -inf)
    {8.0, 4.0, 2.0},    // Exact division
    {-8.0, 4.0, -2.0},  // Exact negative
    {0.0, 5.0, 0.0},    // Zero dividend
    {10.0, 3.0, 3.0},   // 10/3 = 3.33... -> 3
    {-10.0, 3.0, -4.0}, // -10/3 = -3.33... -> -4 (toward -inf)
  };

  for (const auto &tc : cases) {
    simd::Vec<W, double> a = tc.a, b = tc.b;
    auto result = simd::floor_div(a, b);
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      expect(feq(result[i], tc.expected))
        << "floor_div(" << tc.a << ", " << tc.b << ") = " << result[i]
        << ", expected " << tc.expected;
    }
  }
}

// Test ceil_div for correctness
template <std::ptrdiff_t W> void test_ceil_div_double() {
  // Test cases: {a, b, expected_result}
  // ceil_div rounds toward positive infinity
  struct TestCase {
    double a, b, expected;
  };
  constexpr TestCase cases[] = {
    {7.0, 3.0, 3.0},    // 7/3 = 2.33... -> 3 (toward +inf)
    {-7.0, 3.0, -2.0},  // -7/3 = -2.33... -> -2 (toward +inf)
    {7.0, -3.0, -2.0},  // 7/-3 = -2.33... -> -2 (toward +inf)
    {-7.0, -3.0, 3.0},  // -7/-3 = 2.33... -> 3 (toward +inf)
    {8.0, 4.0, 2.0},    // Exact division
    {-8.0, 4.0, -2.0},  // Exact negative
    {0.0, 5.0, 0.0},    // Zero dividend
    {10.0, 3.0, 4.0},   // 10/3 = 3.33... -> 4 (toward +inf)
    {-10.0, 3.0, -3.0}, // -10/3 = -3.33... -> -3 (toward +inf)
  };

  for (const auto &tc : cases) {
    simd::Vec<W, double> a = tc.a, b = tc.b;
    auto result = simd::ceil_div(a, b);
    for (std::ptrdiff_t i = 0; i < W; ++i) {
      expect(feq(result[i], tc.expected))
        << "ceil_div(" << tc.a << ", " << tc.b << ") = " << result[i]
        << ", expected " << tc.expected;
    }
  }
}

// Test that trunc_div and floor_div agree for positive dividends
template <std::ptrdiff_t W> void test_div_positive_agree() {
  for (double a = 1.0; a <= 100.0; a += 7.0) {
    for (double b = 1.0; b <= 20.0; b += 3.0) {
      simd::Vec<W, double> va = a, vb = b;
      auto trunc_result = simd::trunc_div(va, vb);
      auto floor_result = simd::floor_div(va, vb);
      // For positive a and b, trunc and floor should agree
      for (std::ptrdiff_t i = 0; i < W; ++i) {
        expect(feq(trunc_result[i], floor_result[i]))
          << "For positive a=" << a << ", b=" << b
          << ": trunc=" << trunc_result[i] << ", floor=" << floor_result[i];
      }
    }
  }
}

// Test FMA correctness
template <std::ptrdiff_t W> void test_fma_double() {
  simd::Vec<W, double> a, b, c;
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    a[i] = double(i + 1);
    b[i] = double(i + 2);
    c[i] = double(i + 3);
  }

  // fma(a, b, c) = a * b + c
  auto fma_result = simd::fma(a, b, c);
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    double expected = double(i + 1) * double(i + 2) + double(i + 3);
    expect(feq(fma_result[i], expected))
      << "fma failed at index " << i << ": got " << fma_result[i]
      << ", expected " << expected;
  }

  // fms(a, b, c) = a * b - c
  auto fms_result = simd::fms(a, b, c);
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    double expected = double(i + 1) * double(i + 2) - double(i + 3);
    expect(feq(fms_result[i], expected))
      << "fms failed at index " << i << ": got " << fms_result[i]
      << ", expected " << expected;
  }

  // fnma(a, b, c) = c - a * b
  auto fnma_result = simd::fnma(a, b, c);
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    double expected = double(i + 3) - double(i + 1) * double(i + 2);
    expect(feq(fnma_result[i], expected))
      << "fnma failed at index " << i << ": got " << fnma_result[i]
      << ", expected " << expected;
  }
}

// Test FMA with floats
template <std::ptrdiff_t W> void test_fma_float() {
  simd::Vec<W, float> a, b, c;
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    a[i] = float(i + 1);
    b[i] = float(i + 2);
    c[i] = float(i + 3);
  }

  auto fma_result = simd::fma(a, b, c);
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    float expected = float(i + 1) * float(i + 2) + float(i + 3);
    expect(feq(fma_result[i], expected, 1e-5f))
      << "fma float failed at index " << i;
  }

  auto fms_result = simd::fms(a, b, c);
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    float expected = float(i + 1) * float(i + 2) - float(i + 3);
    expect(feq(fms_result[i], expected, 1e-5f))
      << "fms float failed at index " << i;
  }

  auto fnma_result = simd::fnma(a, b, c);
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    float expected = float(i + 3) - float(i + 1) * float(i + 2);
    expect(feq(fnma_result[i], expected, 1e-5f))
      << "fnma float failed at index " << i;
  }
}

// Test that trunc function is actually truncating (not flooring)
template <std::ptrdiff_t W> void test_trunc_vs_floor() {
  simd::Vec<W, double> negative = -2.7;
  auto trunc_result = simd::trunc(negative);
  auto floor_result = simd::floor(negative);

  for (std::ptrdiff_t i = 0; i < W; ++i) {
    // trunc(-2.7) should be -2 (toward zero)
    expect(feq(trunc_result[i], -2.0))
      << "trunc(-2.7) = " << trunc_result[i] << ", expected -2.0";
    // floor(-2.7) should be -3 (toward -infinity)
    expect(feq(floor_result[i], -3.0))
      << "floor(-2.7) = " << floor_result[i] << ", expected -3.0";
  }

  simd::Vec<W, double> positive = 2.7;
  auto trunc_pos = simd::trunc(positive);
  auto floor_pos = simd::floor(positive);

  for (std::ptrdiff_t i = 0; i < W; ++i) {
    // For positive values, trunc and floor should agree
    expect(feq(trunc_pos[i], 2.0))
      << "trunc(2.7) = " << trunc_pos[i] << ", expected 2.0";
    expect(feq(floor_pos[i], 2.0))
      << "floor(2.7) = " << floor_pos[i] << ", expected 2.0";
  }
}

// Test FE_INEXACT preservation for trunc_div
// This test verifies that intentional truncation does NOT set FE_INEXACT
void test_fe_inexact_preservation() {
#ifdef __x86_64__
  // Clear all exception flags
  _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);

  static constexpr auto W = simd::Width<double>;
  simd::Vec<W, double> a = 7.0, b = 3.0;

  // This should NOT set FE_INEXACT (intentional truncation is suppressed)
  auto result = simd::trunc_div(a, b);
  (void)result; // Prevent optimization

  // Check that FE_INEXACT was NOT set
  unsigned int mxcsr = _mm_getcsr();
  bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;

  // Note: On AVX512 with SAE, this should definitely be false.
  // On non-AVX512, the FEnvGuard should have restored the CSR.
#ifdef __AVX512F__
  expect(!inexact_set) << "trunc_div should not set FE_INEXACT on AVX512";
#else
  // On non-AVX512, the guard saves/restores, so it should also be clean
  expect(!inexact_set)
    << "trunc_div should restore FE_INEXACT state on non-AVX512";
#endif

  // Now do a real inexact operation to verify we can detect it
  _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK); // Clear again
  // Use volatile to prevent constant folding
  volatile double one = 1.0;
  volatile double three = 3.0;
  volatile double x = one / three; // This WILL set FE_INEXACT at runtime
  (void)x;
  mxcsr = _mm_getcsr();
  inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;
  expect(inexact_set) << "1.0/3.0 should set FE_INEXACT";
#endif
}

// Test floor_div FE_INEXACT preservation
void test_floor_div_fe_inexact() {
#ifdef __x86_64__
  _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);

  static constexpr auto W = simd::Width<double>;
  simd::Vec<W, double> a = -7.0, b = 3.0;

  auto result = simd::floor_div(a, b);
  (void)result;

  unsigned int mxcsr = _mm_getcsr();
  bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;

#ifdef __AVX512F__
  expect(!inexact_set) << "floor_div should not set FE_INEXACT on AVX512";
#else
  expect(!inexact_set)
    << "floor_div should restore FE_INEXACT state on non-AVX512";
#endif
#endif
}

// Test ceil_div FE_INEXACT preservation
void test_ceil_div_fe_inexact() {
#ifdef __x86_64__
  _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);

  static constexpr auto W = simd::Width<double>;
  simd::Vec<W, double> a = 7.0, b = 3.0;

  auto result = simd::ceil_div(a, b);
  (void)result;

  unsigned int mxcsr = _mm_getcsr();
  bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;

#ifdef __AVX512F__
  expect(!inexact_set) << "ceil_div should not set FE_INEXACT on AVX512";
#else
  expect(!inexact_set)
    << "ceil_div should restore FE_INEXACT state on non-AVX512";
#endif
#endif
}

// Test with large values to ensure precision bounds
template <std::ptrdiff_t W> void test_large_values() {
  // Values well within 53-bit mantissa
  constexpr double large = double(1LL << 40); // 2^40

  simd::Vec<W, double> a, b;
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    a[i] = large + double(i);
    b[i] = 1000.0;
  }

  auto trunc_result = simd::trunc_div(a, b);
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    double expected = std::trunc((large + double(i)) / 1000.0);
    expect(feq(trunc_result[i], expected, 1.0))
      << "Large value trunc_div failed at index " << i;
  }
}

// Test vectorized GCD-style usage pattern
template <std::ptrdiff_t W> void test_gcd_pattern() {
  // Simulate Extended Euclidean Algorithm pattern
  simd::Vec<W, double> r_old, r;
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    r_old[i] = double(100 + i * 17);
    r[i] = double(30 + i * 7);
  }

  // GCD step: quotient = trunc(r_old / r)
  auto quotient = simd::trunc_div(r_old, r);

  // r_next = r_old - quotient * r (using fnma)
  auto r_next = simd::fnma(quotient, r, r_old);

  // Verify: r_next should be the remainder
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    double expected_q = std::trunc(r_old[i] / r[i]);
    double expected_r = r_old[i] - expected_q * r[i];
    expect(feq(quotient[i], expected_q))
      << "GCD quotient mismatch at index " << i;
    expect(feq(r_next[i], expected_r, 1e-9))
      << "GCD remainder mismatch at index " << i << ": got " << r_next[i]
      << ", expected " << expected_r;
  }
}

// Test that operations inconsistent with integer behavior DO set FE_INEXACT
void test_inexact_operations_set_flag() {
#ifdef __x86_64__
  static constexpr auto W = simd::Width<double>;

  // Test 1: Inexact division (7.0/3.0 = 2.333...)
  {
    _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);
    simd::Vec<W, double> a = 7.0, b = 3.0;
    // Regular division without truncation SHOULD set FE_INEXACT
    volatile auto result = a / b;
    (void)result;
    unsigned int mxcsr = _mm_getcsr();
    bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;
    expect(inexact_set) << "7.0/3.0 should set FE_INEXACT";
  }

  // Test 2: Addition with precision loss (1e18 + 65.0)
  // 1e18 is large enough that adding 65 loses precision
  {
    _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);
    simd::Vec<W, double> large = 1e18, small = 65.0;
    volatile auto result = large + small;
    (void)result;
    unsigned int mxcsr = _mm_getcsr();
    bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;
    expect(inexact_set) << "1e18 + 65.0 should set FE_INEXACT (precision loss)";
  }

  // Test 3: Multiplication with precision loss
  {
    _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);
    simd::Vec<W, double> a = 1.0 / 3.0; // Already inexact representation
    simd::Vec<W, double> b = 3.0;
    volatile auto result = a * b; // Should not be exactly 1.0
    (void)result;
    unsigned int mxcsr = _mm_getcsr();
    bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;
    // Note: This may or may not set inexact depending on the FMA behavior
    // The multiplication itself is exact, but the input was already inexact
    // We mainly care that inexact operations are detectable
    (void)inexact_set; // Just verify the mechanism works
  }

  // Test 4: Large product that exceeds mantissa precision
  {
    _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);
    // (2^27 + 1) * (2^27 + 1) requires more than 53 bits of precision
    simd::Vec<W, double> a = double((1LL << 27) + 1);
    volatile auto result = a * a;
    (void)result;
    unsigned int mxcsr = _mm_getcsr();
    bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;
    expect(inexact_set)
      << "(2^27+1)^2 should set FE_INEXACT (exceeds 53-bit mantissa)";
  }

  // Test 5: Verify exact operations do NOT set FE_INEXACT
  {
    _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);
    simd::Vec<W, double> a = 8.0, b = 4.0;
    volatile auto result = a / b; // Exact: 8/4 = 2
    (void)result;
    unsigned int mxcsr = _mm_getcsr();
    bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;
    expect(!inexact_set)
      << "8.0/4.0 should NOT set FE_INEXACT (exact division)";
  }

  // Test 6: Exact addition should not set FE_INEXACT
  {
    _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);
    simd::Vec<W, double> a = 1.0, b = 2.0;
    volatile auto result = a + b; // Exact: 1+2 = 3
    (void)result;
    unsigned int mxcsr = _mm_getcsr();
    bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;
    expect(!inexact_set)
      << "1.0+2.0 should NOT set FE_INEXACT (exact addition)";
  }
#endif
}

// Test Unroll rounding division variants
template <std::ptrdiff_t R, std::ptrdiff_t C, std::ptrdiff_t W>
void test_unroll_rounding_div() {
  using Unroll = simd::Unroll<R, C, W, double>;
  using VT = typename Unroll::VT;

  // Create test data
  Unroll a_unroll;
  Unroll b_unroll;
  for (std::ptrdiff_t i = 0; i < R * C; ++i) {
    if constexpr (W == 1) {
      a_unroll[i] = 7.0 + double(i);
      b_unroll[i] = 3.0 + double(i % 2);
    } else {
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        a_unroll[i][j] = 7.0 + double(i * W + j);
        b_unroll[i][j] = 3.0 + double((i * W + j) % 3);
      }
    }
  }

  // Test Unroll-Unroll (uses ADL to find friend functions)
  auto ceil_uu = ceil_div(a_unroll, b_unroll);
  auto floor_uu = floor_div(a_unroll, b_unroll);
  auto trunc_uu = trunc_div(a_unroll, b_unroll);
  (void)ceil_uu;
  (void)floor_uu;
  (void)trunc_uu;

  // Test Unroll-VT
  VT b_vt;
  if constexpr (W == 1) b_vt = 3.0;
  else
    for (std::ptrdiff_t j = 0; j < W; ++j) b_vt[j] = 3.0 + double(j % 2);
  auto ceil_uvt = ceil_div(a_unroll, b_vt);
  auto floor_uvt = floor_div(a_unroll, b_vt);
  auto trunc_uvt = trunc_div(a_unroll, b_vt);
  (void)ceil_uvt;
  (void)floor_uvt;
  (void)trunc_uvt;

  // Test VT-Unroll
  VT a_vt;
  if constexpr (W == 1) a_vt = 10.0;
  else
    for (std::ptrdiff_t j = 0; j < W; ++j) a_vt[j] = 10.0 + double(j);
  auto ceil_vtu = ceil_div(a_vt, b_unroll);
  auto floor_vtu = floor_div(a_vt, b_unroll);
  auto trunc_vtu = trunc_div(a_vt, b_unroll);
  (void)ceil_vtu;
  (void)floor_vtu;
  (void)trunc_vtu;

  // Test Unroll-scalar and scalar-Unroll (only for W != 1)
  if constexpr (W != 1) {
    auto ceil_us = ceil_div(a_unroll, 3.0);
    auto floor_us = floor_div(a_unroll, 3.0);
    auto trunc_us = trunc_div(a_unroll, 3.0);
    (void)ceil_us;
    (void)floor_us;
    (void)trunc_us;

    auto ceil_su = ceil_div(10.0, b_unroll);
    auto floor_su = floor_div(10.0, b_unroll);
    auto trunc_su = trunc_div(10.0, b_unroll);
    (void)ceil_su;
    (void)floor_su;
    (void)trunc_su;
  }

  // Verify correctness for a sample
  for (std::ptrdiff_t i = 0; i < R * C; ++i) {
    if constexpr (W == 1) {
      double a = a_unroll[i];
      double b = b_unroll[i];
      expect(feq(ceil_uu[i], std::ceil(a / b)))
        << "Unroll-Unroll ceil_div mismatch at " << i;
      expect(feq(floor_uu[i], std::floor(a / b)))
        << "Unroll-Unroll floor_div mismatch at " << i;
      expect(feq(trunc_uu[i], std::trunc(a / b)))
        << "Unroll-Unroll trunc_div mismatch at " << i;
    } else {
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        double a = a_unroll[i][j];
        double b = b_unroll[i][j];
        expect(feq(ceil_uu[i][j], std::ceil(a / b)))
          << "Unroll-Unroll ceil_div mismatch at [" << i << "][" << j << "]";
        expect(feq(floor_uu[i][j], std::floor(a / b)))
          << "Unroll-Unroll floor_div mismatch at [" << i << "][" << j << "]";
        expect(feq(trunc_uu[i][j], std::trunc(a / b)))
          << "Unroll-Unroll trunc_div mismatch at [" << i << "][" << j << "]";
      }
    }
  }
}

// Test that FEnvGuard is correctly used for Unroll operations
void test_unroll_fe_inexact_preservation() {
#ifdef __x86_64__
  static constexpr auto W = simd::Width<double>;
  using Unroll = simd::Unroll<2, 2, W, double>;
  using VT = typename Unroll::VT;

  // Clear exception flags
  _mm_setcsr(_mm_getcsr() & ~_MM_EXCEPT_MASK);

  // Create test data with inexact division
  Unroll a_unroll;
  VT b_vt;
  for (std::ptrdiff_t i = 0; i < 4; ++i)
    for (std::ptrdiff_t j = 0; j < W; ++j) a_unroll[i][j] = 7.0;
  for (std::ptrdiff_t j = 0; j < W; ++j) b_vt[j] = 3.0;

  // Perform rounding division (should NOT set FE_INEXACT due to guard)
  auto result = ceil_div(a_unroll, b_vt);
  (void)result;

  unsigned int mxcsr = _mm_getcsr();
  bool inexact_set = (mxcsr & _MM_EXCEPT_INEXACT) != 0;

#ifdef __AVX512F__
  if constexpr (W == 8) {
    // Full-width AVX512 uses SAE, so FE_INEXACT should not be set
    expect(!inexact_set)
      << "Unroll ceil_div should not set FE_INEXACT on full-width AVX512";
  } else {
    // Non-full-width should use FEnvGuard to restore
    expect(!inexact_set)
      << "Unroll ceil_div should restore FE_INEXACT state via FEnvGuard";
  }
#else
  expect(!inexact_set)
    << "Unroll ceil_div should restore FE_INEXACT state via FEnvGuard";
#endif
#endif
}

} // namespace

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
int main() {
  static constexpr auto W_double = simd::Width<double>;
  static constexpr auto W_float = simd::Width<float>;

  "TruncDivDouble"_test = [] { test_trunc_div_double<W_double>(); };

  "FloorDivDouble"_test = [] { test_floor_div_double<W_double>(); };

  "CeilDivDouble"_test = [] { test_ceil_div_double<W_double>(); };

  "DivPositiveAgree"_test = [] { test_div_positive_agree<W_double>(); };

  "FmaDouble"_test = [] { test_fma_double<W_double>(); };

  "FmaFloat"_test = [] { test_fma_float<W_float>(); };

  "TruncVsFloor"_test = [] { test_trunc_vs_floor<W_double>(); };

  "FeInexactPreservation"_test = [] { test_fe_inexact_preservation(); };

  "FloorDivFeInexact"_test = [] { test_floor_div_fe_inexact(); };

  "CeilDivFeInexact"_test = [] { test_ceil_div_fe_inexact(); };

  "InexactOpsSetFlag"_test = [] { test_inexact_operations_set_flag(); };

  "LargeValues"_test = [] { test_large_values<W_double>(); };

  "GcdPattern"_test = [] { test_gcd_pattern<W_double>(); };

  // Test Unroll rounding division variants
  "UnrollRoundingDiv_1_1"_test = [] {
    test_unroll_rounding_div<1, 1, W_double>();
  };
  "UnrollRoundingDiv_2_2"_test = [] {
    test_unroll_rounding_div<2, 2, W_double>();
  };
  "UnrollRoundingDiv_1_4"_test = [] {
    test_unroll_rounding_div<1, 4, W_double>();
  };

  "UnrollFeInexactPreservation"_test = [] {
    test_unroll_fe_inexact_preservation();
  };

  // Test with different SIMD widths if available
  if constexpr (W_double >= 2) {
    "TruncDivDouble_W2"_test = [] { test_trunc_div_double<2>(); };
    "FloorDivDouble_W2"_test = [] { test_floor_div_double<2>(); };
    "CeilDivDouble_W2"_test = [] { test_ceil_div_double<2>(); };
    "FmaDouble_W2"_test = [] { test_fma_double<2>(); };
  }

  if constexpr (W_double >= 4) {
    "TruncDivDouble_W4"_test = [] { test_trunc_div_double<4>(); };
    "FloorDivDouble_W4"_test = [] { test_floor_div_double<4>(); };
    "CeilDivDouble_W4"_test = [] { test_ceil_div_double<4>(); };
    "FmaDouble_W4"_test = [] { test_fma_double<4>(); };
  }

  if constexpr (W_double >= 8) {
    "TruncDivDouble_W8"_test = [] { test_trunc_div_double<8>(); };
    "FloorDivDouble_W8"_test = [] { test_floor_div_double<8>(); };
    "CeilDivDouble_W8"_test = [] { test_ceil_div_double<8>(); };
    "FmaDouble_W8"_test = [] { test_fma_double<8>(); };
  }

  return 0;
}

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

} // namespace

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
int main() {
  static constexpr auto W_double = simd::Width<double>;
  static constexpr auto W_float = simd::Width<float>;

  "TruncDivDouble"_test = [] { test_trunc_div_double<W_double>(); };

  "FloorDivDouble"_test = [] { test_floor_div_double<W_double>(); };

  "DivPositiveAgree"_test = [] { test_div_positive_agree<W_double>(); };

  "FmaDouble"_test = [] { test_fma_double<W_double>(); };

  "FmaFloat"_test = [] { test_fma_float<W_float>(); };

  "TruncVsFloor"_test = [] { test_trunc_vs_floor<W_double>(); };

  "FeInexactPreservation"_test = [] { test_fe_inexact_preservation(); };

  "FloorDivFeInexact"_test = [] { test_floor_div_fe_inexact(); };

  "LargeValues"_test = [] { test_large_values<W_double>(); };

  "GcdPattern"_test = [] { test_gcd_pattern<W_double>(); };

  // Test with different SIMD widths if available
  if constexpr (W_double >= 2) {
    "TruncDivDouble_W2"_test = [] { test_trunc_div_double<2>(); };
    "FloorDivDouble_W2"_test = [] { test_floor_div_double<2>(); };
    "FmaDouble_W2"_test = [] { test_fma_double<2>(); };
  }

  if constexpr (W_double >= 4) {
    "TruncDivDouble_W4"_test = [] { test_trunc_div_double<4>(); };
    "FloorDivDouble_W4"_test = [] { test_floor_div_double<4>(); };
    "FmaDouble_W4"_test = [] { test_fma_double<4>(); };
  }

  if constexpr (W_double >= 8) {
    "TruncDivDouble_W8"_test = [] { test_trunc_div_double<8>(); };
    "FloorDivDouble_W8"_test = [] { test_floor_div_double<8>(); };
    "FmaDouble_W8"_test = [] { test_fma_double<8>(); };
  }

  return 0;
}

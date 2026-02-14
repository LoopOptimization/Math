import Testing;

import BaseUtils;
import StaticArray;
import SIMD;
import std;

using namespace testing;

void testBasicAssertions() {

  for (int i = 0; i++ < 32;) {
    double di = i;
    for (int j = 0; j++ < i;) {
      double dj = j;
      auto [x, y] = ::math::lower_bound_factor(di, dj);
      expect(di == x * y);
      expect(x <= dj);
      expect(std::round(x) == x);
      expect(std::round(y) == y);
    }
  }
}

void testSIMDAssertions() {
  // Test with SIMD vectors of different widths
  constexpr std::ptrdiff_t W = simd::Width<double>;

  // Create test vectors for values and bounds
  for (int base_i = W; base_i <= 32; base_i += W) {
    ::math::SVector<double, W> di_vec;
    for (std::ptrdiff_t k = 0; k < W; ++k) di_vec[k] = base_i - k;

    for (int base_j = W; base_j <= base_i; base_j += W) {
      ::math::SVector<double, W> dj_vec;
      for (std::ptrdiff_t k = 0; k < W; ++k)
        dj_vec[k] = std::min(base_j - k, base_i - k - 1);

      // Get SIMD representation
      auto di_simd = di_vec.simd();
      auto dj_simd = dj_vec.simd();

      // Call lower_bound_factor with SIMD types
      auto [x_simd, y_simd] = ::math::lower_bound_factor(di_simd, dj_simd);

      // Verify results for each element
      for (std::ptrdiff_t k = 0; k < W; ++k) {
        double di = di_vec[k];
        double dj = dj_vec[k];
        double x = x_simd.vec_[k];
        double y = y_simd.vec_[k];

        if (dj > 0) { // Only test when dj is positive
          expect(std::abs(di - x * y) < 1e-10)
            << "di=" << di << " x=" << x << " y=" << y;
          expect(x <= dj + 1e-10) << "x=" << x << " dj=" << dj;
          expect(std::abs(std::round(x) - x) < 1e-10) << "x=" << x;
          expect(std::abs(std::round(y) - y) < 1e-10) << "y=" << y;
        }
      }
    }
  }
}

void testSIMDUnrollAssertions() {
  // Test with larger unrolled SIMD types
  constexpr std::ptrdiff_t W = simd::Width<double>;
  constexpr std::ptrdiff_t U = 4; // Unroll factor

  using VecType = ::math::StaticArray<double, U, W, false>;

  // Create test matrices
  for (int base_i = W; base_i <= 32; base_i += W) {
    VecType di_mat;
    for (std::ptrdiff_t r = 0; r < U; ++r)
      for (std::ptrdiff_t c = 0; c < W; ++c) di_mat[r, c] = base_i - r * W - c;

    for (int base_j = W; base_j <= base_i; base_j += W) {
      VecType dj_mat;
      for (std::ptrdiff_t r = 0; r < U; ++r)
        for (std::ptrdiff_t c = 0; c < W; ++c)
          dj_mat[r, c] = std::min(base_j - r * W - c, base_i - r * W - c - 1);

      // Get SIMD Unroll representation
      auto di_unroll = di_mat.simd();
      auto dj_unroll = dj_mat.simd();

      // Call lower_bound_factor with Unroll types
      auto [x_unroll, y_unroll] =
        ::math::lower_bound_factor(di_unroll, dj_unroll);

      // Verify results for each element
      for (std::ptrdiff_t r = 0; r < U; ++r) {
        for (std::ptrdiff_t c = 0; c < W; ++c) {
          double di = di_mat[r, c];
          double dj = dj_mat[r, c];
          double x = x_unroll[r, 0][c];
          double y = y_unroll[r, 0][c];

          if (dj > 0 && di > 0) { // Only test when values are positive
            expect(std::abs(di - x * y) < 1e-10)
              << "di=" << di << " x=" << x << " y=" << y << " at [" << r << ","
              << c << "]";
            expect(x <= dj + 1e-10)
              << "x=" << x << " dj=" << dj << " at [" << r << "," << c << "]";
            expect(std::abs(std::round(x) - x) < 1e-10)
              << "x=" << x << " at [" << r << "," << c << "]";
            expect(std::abs(std::round(y) - y) < 1e-10)
              << "y=" << y << " at [" << r << "," << c << "]";
          }
        }
      }
    }
  }
}

int main() {
  "FactorLowerBound BasicAssertions"_test = [] { testBasicAssertions(); };
  "FactorLowerBound SIMDAssertions"_test = [] { testSIMDAssertions(); };
  "FactorLowerBound SIMDUnrollAssertions"_test = [] {
    testSIMDUnrollAssertions();
  };
  return 0;
}

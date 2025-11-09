import boost.ut;
import Array;
import ExprTemplates;
import ManagedArray;
import std;

using namespace boost::ut;

auto main() -> int {
  "select with Array comparison"_test = [] {
    // Test select with mask from comparison
    ::math::Vector<int> a{std::array{1, 2, 3, 4, 5}};
    ::math::Vector<int> b{std::array{5, 4, 3, 2, 1}};
    ::math::Vector<int> result{::math::length(5)};

    // Select: where a > 2, choose a, otherwise choose b
    result << ::math::select(a > 2, a, b);

    expect(result[0] == 5_i);  // a[0] = 1 !> 2, so b[0] = 5
    expect(result[1] == 4_i);  // a[1] = 2 !> 2, so b[1] = 4
    expect(result[2] == 3_i);  // a[2] = 3 > 2, so a[2] = 3
    expect(result[3] == 4_i);  // a[3] = 4 > 2, so a[3] = 4
    expect(result[4] == 5_i);  // a[4] = 5 > 2, so a[4] = 5
  };

  "select between two arrays"_test = [] {
    ::math::Vector<int> a{std::array{10, 20, 30, 40}};
    ::math::Vector<int> b{std::array{15, 5, 35, 25}};
    ::math::Vector<int> result{::math::length(4)};

    // Select based on comparison: where a > b, choose a, otherwise choose b
    result << ::math::select(a > b, a, b);

    expect(result[0] == 15_i);  // 10 !> 15, so b[0] = 15
    expect(result[1] == 20_i);  // 20 > 5, so a[1] = 20
    expect(result[2] == 35_i);  // 30 !> 35, so b[2] = 35
    expect(result[3] == 40_i);  // 40 > 25, so a[3] = 40
  };

  "select with scalar broadcast"_test = [] {
    ::math::Vector<double> a{std::array{1.0, 2.0, 3.0, 4.0}};
    ::math::Vector<double> result{::math::length(4)};

    // Select: where a > 2.5, choose a, otherwise choose 0.0
    result << ::math::select(a > 2.5, a, 0.0);

    expect(approx(result[0], 0.0, 1e-10));
    expect(approx(result[1], 0.0, 1e-10));
    expect(approx(result[2], 3.0, 1e-10));
    expect(approx(result[3], 4.0, 1e-10));
  };

  "select with scalar on both sides"_test = [] {
    ::math::Vector<int> cond{std::array{1, 0, 1, 0, 1}};
    ::math::Vector<int> result{::math::length(5)};

    // Select: where cond != 0, choose 100, otherwise choose -100
    result << ::math::select(cond != 0, 100, -100);

    expect(result[0] == 100_i);
    expect(result[1] == -100_i);
    expect(result[2] == 100_i);
    expect(result[3] == -100_i);
    expect(result[4] == 100_i);
  };

  "select with floating point arrays"_test = [] {
    ::math::Vector<double> x{std::array{-1.5, 2.5, -3.5, 4.5, 0.0}};
    ::math::Vector<double> positive{::math::length(5)};

    // Clamp negative values to 0.0
    positive << ::math::select(x >= 0.0, x, 0.0);

    expect(approx(positive[0], 0.0, 1e-10));
    expect(approx(positive[1], 2.5, 1e-10));
    expect(approx(positive[2], 0.0, 1e-10));
    expect(approx(positive[3], 4.5, 1e-10));
    expect(approx(positive[4], 0.0, 1e-10));
  };

  "select with complex expression"_test = [] {
    ::math::Vector<int> a{std::array{1, 2, 3, 4, 5}};
    ::math::Vector<int> b{std::array{10, 20, 30, 40, 50}};
    ::math::Vector<int> c{std::array{5, 15, 25, 35, 45}};
    ::math::Vector<int> result{::math::length(5)};

    // Select based on condition involving expressions
    result << ::math::select((a + b) > c, a * 2, b / 2);

    // a[0]+b[0]=11 > c[0]=5, so a[0]*2=2
    expect(result[0] == 2_i);
    // a[1]+b[1]=22 > c[1]=15, so a[1]*2=4
    expect(result[1] == 4_i);
    // a[2]+b[2]=33 > c[2]=25, so a[2]*2=6
    expect(result[2] == 6_i);
    // a[3]+b[3]=44 > c[3]=35, so a[3]*2=8
    expect(result[3] == 8_i);
    // a[4]+b[4]=55 > c[4]=45, so a[4]*2=10
    expect(result[4] == 10_i);
  };

  "select with equality comparison"_test = [] {
    ::math::Vector<int> a{std::array{1, 2, 3, 2, 1}};
    ::math::Vector<int> result{::math::length(5)};

    // Find all elements equal to 2
    result << ::math::select(a == 2, 999, 0);

    expect(result[0] == 0_i);
    expect(result[1] == 999_i);
    expect(result[2] == 0_i);
    expect(result[3] == 999_i);
    expect(result[4] == 0_i);
  };

  "select with long arrays for SIMD"_test = [] {
    // Create longer arrays to ensure SIMD paths are exercised
    constexpr std::ptrdiff_t N = 64;
    ::math::Vector<double> a{::math::length(N)};
    ::math::Vector<double> b{::math::length(N)};
    ::math::Vector<double> result{::math::length(N)};

    // Initialize with pattern
    for (std::ptrdiff_t i = 0; i < N; ++i) {
      a[i] = static_cast<double>(i);
      b[i] = static_cast<double>(N - i);
    }

    // Select based on midpoint
    result << ::math::select(a < (N / 2), a, b);

    for (std::ptrdiff_t i = 0; i < N; ++i) {
      if (i < N / 2) {
        expect(approx(result[i], static_cast<double>(i), 1e-10));
      } else {
        expect(approx(result[i], static_cast<double>(N - i), 1e-10));
      }
    }
  };

  "select with nested operations"_test = [] {
    ::math::Vector<int> a{std::array{1, 2, 3, 4, 5}};
    ::math::Vector<int> b{std::array{2, 4, 6, 8, 10}};
    ::math::Vector<int> result{::math::length(5)};

    // Nested select: abs-like operation using select
    // result = (a > 0) ? a : -a, then select based on comparison with b
    auto abs_a = ::math::select(a > 0, a, -a);
    result << ::math::select(abs_a < b, abs_a, b);

    expect(result[0] == 1_i);
    expect(result[1] == 2_i);
    expect(result[2] == 3_i);
    expect(result[3] == 4_i);
    expect(result[4] == 5_i);
  };

  "select with inequality operators"_test = [] {
    ::math::Vector<double> temps{std::array{-5.0, 0.0, 15.0, 25.0, 35.0}};
    ::math::Vector<double> result{::math::length(5)};

    // Classify temperatures: >= 20 is hot (1.0), otherwise cold (0.0)
    result << ::math::select(temps >= 20.0, 1.0, 0.0);

    expect(approx(result[0], 0.0, 1e-10));
    expect(approx(result[1], 0.0, 1e-10));
    expect(approx(result[2], 0.0, 1e-10));
    expect(approx(result[3], 1.0, 1e-10));
    expect(approx(result[4], 1.0, 1e-10));
  };

  "select with less than comparison"_test = [] {
    ::math::Vector<int> values{std::array{-10, -5, 0, 5, 10}};
    ::math::Vector<int> result{::math::length(5)};

    // Select based on < 0 threshold
    result << ::math::select(values < 0, -1, 1);

    expect(result[0] == -1_i);
    expect(result[1] == -1_i);
    expect(result[2] == 1_i);
    expect(result[3] == 1_i);
    expect(result[4] == 1_i);
  };

  "select with less than or equal comparison"_test = [] {
    ::math::Vector<int> values{std::array{1, 5, 10, 15, 20}};
    ::math::Vector<int> result{::math::length(5)};

    // Select based on <= 10 threshold
    result << ::math::select(values <= 10, values, 10);

    expect(result[0] == 1_i);
    expect(result[1] == 5_i);
    expect(result[2] == 10_i);
    expect(result[3] == 10_i);   // Clamped to 10
    expect(result[4] == 10_i);   // Clamped to 10
  };

  return 0;
}

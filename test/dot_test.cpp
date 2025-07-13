import boost.ut;

import ArrayConcepts;
import ManagedArray;
import Reductions;
import StaticArray;
import std;

using namespace ::math;
using namespace boost::ut;

int main() {
  "DotProduct BasicAssertions"_test = [] {
    // Test dynamic vectors constructed from arrays
    Vector<double> a(std::array<double, 3>{1.0, 2.0, 3.0});
    Vector<double> b(std::array<double, 3>{4.0, 5.0, 6.0});
    expect(dot(a, b) == 32.0_d);

    // Test static arrays
    SVector<double, 3> sa{1.0, 2.0, 3.0};
    SVector<double, 3> sb{4.0, 5.0, 6.0};
    expect(dot(sa, sb) == 32.0_d);

    // Test mixed types
    expect(dot(a, sa) == 14.0_d);
    expect(dot(sa, b) == 32.0_d);
  };

  "DotProduct TransposedViews"_test = [] {
    // Test with matrix views
    StaticArray<double, 2, 3> A{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    auto a_row = A[0, _]; // Row vector [1, 2, 3]
    auto a_col = A[_, 0]; // Column vector [1, 4]

    Vector<double> b(std::array<double, 3>{1.0, 2.0, 3.0});
    Vector<double> c(std::array<double, 2>{1.0, 4.0});

    expect(dot(a_row, b) == 14.0_d); // 1*1 + 2*2 + 3*3 = 14
    expect(dot(a_col, c) == 17.0_d); // 1*1 + 4*4 = 17
  };

  "DotProduct Flattening"_test = [] {
    // Test flattening behavior with compatible layouts
    StaticArray<int, 2, 2> a{1, 2, 3, 4};
    StaticArray<int, 2, 2> b{2, 3, 4, 5};
    static_assert(StaticArray<int, 2, 2>::can_flatten);
    static_assert(!are_both_transposed<decltype(a), decltype(b)>());
    static_assert(can_flatten_together<decltype(a), decltype(b)>());
    static_assert(CanFlatten<decltype(a)>);
    expect(dot(a.flatview(), b.flatview()) == 40); // 1*2 + 2*3 + 3*4 + 4*5 = 40
    expect(dot(a, b) == 40);                       // 1*2 + 2*3 + 3*4 + 4*5 = 40

    // Test with different element types
    SVector<float, 3> af{1.0f, 2.0f, 3.0f};
    SVector<float, 3> bf{4.0f, 5.0f, 6.0f};
    expect(dot(af, bf) == 32.0f);
  };

  "DotProduct EdgeCases"_test = [] {
    // Test single element vectors
    SVector<double, 1> single_a{5.0};
    SVector<double, 1> single_b{3.0};
    expect(dot(single_a, single_b) == 15.0_d);

    // Test zero vectors
    SVector<double, 3> zero_a{0.0, 0.0, 0.0};
    SVector<double, 3> zero_b{1.0, 2.0, 3.0};
    expect(dot(zero_a, zero_b) == 0.0_d);
  };

  "DotProduct VariableSizes"_test = [] {
    Vector<double> a(length(64));
    Vector<double> b(length(64));

    // Fill with simple patterns
    for (std::ptrdiff_t i = 0; i < 64; ++i) {
      a[i] = static_cast<double>(i + 1); // [1, 2, 3, ...]
      b[i] = static_cast<double>(i + 1); // [1, 2, 3, ...]
    }
    // Test various sizes to cover different SIMD paths
    for (std::ptrdiff_t size = 1; size <= 64; ++size) {

      // Expected result: 1^2 + 2^2 + 3^2 + ... + size^2 =
      // size*(size+1)*(2*size+1)/6
      double expected =
        static_cast<double>(size * (size + 1) * (2 * size + 1)) / 6.0;
      expect(dot(a[_(size)], b[_(size)]) == expected);
    }
  };
}

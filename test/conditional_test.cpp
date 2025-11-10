import boost.ut;
import Array;
import ArrayParse;
import ExprTemplates;
import ManagedArray;
import std;

using namespace ::math;
using utils::operator""_mat;
using namespace boost::ut;

int main() {
  "Conditional with SIMD mask test"_test = [] {
    // Test conditional operation with vector masks
    // This tests the code path where W > 1 (SIMD width > 1)
    // and ensures we're accessing c.mask_ instead of c.vec_

    auto a = "[10 20 30 40 50 60 70 80]"_mat;
    auto b = "[5 15 25 35 45 55 65 75]"_mat;

    // Create a mask: elements where a > b
    auto mask = a > b;

    // Test conditional addition: if (mask) then a + b else a
    Vector<std::int64_t> result1{conditional(std::plus<>{}, mask, a, b)};
    auto expected1 = "[15 35 55 75 95 115 135 155]"_mat;
    expect(result1 == expected1);

    // Test conditional subtraction: if (mask) then a - b else a
    Vector<std::int64_t> result2{conditional(std::minus<>{}, mask, a, b)};
    auto expected2 = "[5 5 5 5 5 5 5 5]"_mat;
    expect(result2 == expected2);

    // Test conditional multiplication: if (mask) then a * b else a
    Vector<std::int64_t> result3{conditional(std::multiplies<>{}, mask, a, b)};
    auto expected3 = "[50 300 750 1400 2250 3300 4550 6000]"_mat;
    expect(result3 == expected3);
  };

  "Conditional with different mask patterns"_test = [] {
    // Test with different comparison patterns
    auto x = "[1 2 3 4 5 6 7 8 9 10]"_mat;
    auto y = "[5 5 5 5 5 5 5 5 5 5]"_mat;

    // Mask: x < y (first 4 elements true)
    auto mask_lt = x < y;
    Vector<std::int64_t> result_lt{conditional(std::plus<>{}, mask_lt, x, y)};
    auto expected_lt = "[6 7 8 9 5 6 7 8 9 10]"_mat;
    expect(result_lt == expected_lt);

    // Mask: x >= y (last 6 elements true)
    auto mask_ge = x >= y;
    Vector<std::int64_t> result_ge{conditional(std::plus<>{}, mask_ge, x, y)};
    auto expected_ge = "[1 2 3 4 10 11 12 13 14 15]"_mat;
    expect(result_ge == expected_ge);
  };

  "Conditional with scalar operands"_test = [] {
    // Test mixing vector mask with scalar operands
    auto v = "[2 4 6 8 10 12 14 16]"_mat;
    auto threshold = "[5 5 5 5 5 5 5 5]"_mat;
    auto mask = v > threshold;

    // Test with scalar second operand
    Vector<std::int64_t> result{conditional(std::plus<>{}, mask, v, 100)};
    auto expected = "[2 4 106 108 110 112 114 116]"_mat;
    expect(result == expected);
  };

  "Conditional with expression templates"_test = [] {
    // Test conditional with expression template inputs
    auto a = "[10 20 30 40]"_mat;
    auto b = "[5 15 25 35]"_mat;

    // Mask from expression: (a * 2) > (b * 3)
    auto mask = (a * 2) > (b * 3);

    Vector<std::int64_t> result{conditional(std::plus<>{}, mask, a, b)};
    // a * 2 = [20 40 60 80]
    // b * 3 = [15 45 75 105]
    // mask = [true false false false]
    // result[0] = 10 + 5 = 15 (mask true)
    // result[1] = 20 (mask false, keep a)
    // result[2] = 30 (mask false, keep a)
    // result[3] = 40 (mask false, keep a)
    auto expected = "[15 20 30 40]"_mat;
    expect(result == expected);
  };

  "Conditional with matrix operations"_test = [] {
    // Test with matrices to ensure 2D indexing works
    auto A = "[10 20 30; 40 50 60]"_mat;
    auto B = "[15 15 15; 45 45 45]"_mat;

    auto mask = A > B;

    IntMatrix<> result{conditional(std::plus<>{}, mask, A, B)};
    // Mask: [false true true; false true true]
    // Result: [10 35 45; 40 95 105]
    auto expected = "[10 35 45; 40 95 105]"_mat;
    expect(result == expected);
  };

  return 0;
}

import boost.ut;

import Arena;
import Comparisons;
import ManagedArray;
import StaticArray;
import std;

using namespace ::math;
using namespace boost::ut;

int main() {
  "find_if vector tests"_test = [] {
    // Test with StaticArray (vector)
    StaticArray<std::int64_t, 1, 8> vec{1, 2, 3, 4, 5, 6, 7, 8};

    // Find first element > 3
    auto idx = find_first(vec, [](auto x) { return x > 3; });
    expect(idx == 3_i); // Should find 4 at index 3

    // Find first element == 7
    idx = find_first(vec, [](auto x) { return x == 7; });
    expect(idx == 6_i); // Should find 7 at index 6

    // Find first element > 10 (doesn't exist)
    idx = find_first(vec, [](auto x) { return x > 10; });
    expect(idx == -1_i); // Should return -1

    // Find first element < 0 (doesn't exist)
    idx = find_first(vec, [](auto x) { return x < 0; });
    expect(idx == -1_i); // Should return -1

    // Find first element == 1
    idx = find_first(vec, [](auto x) { return x == 1; });
    expect(idx == 0_i); // Should find 1 at index 0
  };

  "find_if matrix tests"_test = [] {
    // Test with StaticArray (matrix)
    StaticArray<std::int64_t, 3, 4> mat{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    // Find first element > 7
    auto idx = find_first(mat, [](auto x) { return x > 7; });
    expect(idx == 7_i); // Should find 8 at linear index 7 (row 1, col 3)

    idx = find_first(mat.t(), [](auto x) { return x > 7; });
    expect(idx == 2_i); // Should find 8 at linear index 7 (row 1, col 3)

    // Find first element == 10
    idx = find_first(mat, [](auto x) { return x == 10; });
    expect(idx == 9_i); // Should find 10 at linear index 9 (row 2, col 1)

    // Find first element > 15 (doesn't exist)
    idx = find_first(mat, [](auto x) { return x > 15; });
    expect(idx == -1_i); // Should return -1

    // Find first element == 1
    idx = find_first(mat, [](auto x) { return x == 1; });
    expect(idx == 0_i); // Should find 1 at linear index 0
  };

  "find_if dynamic array tests"_test = [] {
    // Test with ManagedArray (dynamic size)
    ManagedArray<std::int64_t, Length<>> vec(length(10));
    for (std::ptrdiff_t i = 0; i < 10; ++i)
      vec[i] = i * 2; // [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

    // Find first odd element (doesn't exist)
    // auto idx = find_if(vec, [](auto x) { return x % 2 == 1; });
    // expect(idx == -1_i); // Should return -1

    // Find first element >= 10
    auto idx = find_first(vec, [](auto x) { return x >= 10; });
    expect(idx == 5_i); // Should find 10 at index 5

    // Find first element == 0
    idx = find_first(vec, [](auto x) { return x == 0; });
    expect(idx == 0_i); // Should find 0 at index 0
  };

  "find_if SIMD vector tests"_test = [] {
    // Test with larger array to trigger SIMD paths
    SVector<std::int64_t, 32> large_vec;
    for (std::ptrdiff_t i = 0; i < 32; ++i)
      large_vec[i] = i + 1; // [1, 2, 3, ..., 32]

    // Find first element > 20
    auto idx = find_first(large_vec, [](auto x) { return x > 20; });
    expect(idx == 20_i); // Should find 21 at index 20

    // Find first element == 15
    idx = find_first(large_vec, [](auto x) { return x == 15; });
    expect(idx == 14_i); // Should find 15 at index 14

    // Find first element > 50 (doesn't exist)
    idx = find_first(large_vec, [](auto x) { return x > 50; });
    expect(idx == -1_i); // Should return -1
  };

  "find_if negative values tests"_test = [] {
    // Test with negative values
    SVector<std::int64_t, 8> vec{-4, -3, -2, -1, 0, 1, 2, 3};

    // Find first positive element
    auto idx = find_first(vec, [](auto x) { return x > 0; });
    expect(idx == 5_i); // Should find 1 at index 5

    // Find first negative element
    idx = find_first(vec, [](auto x) { return x < 0; });
    expect(idx == 0_i); // Should find -4 at index 0

    // Find first element == 0
    idx = find_first(vec, [](auto x) { return x == 0; });
    expect(idx == 4_i); // Should find 0 at index 4
  };
}

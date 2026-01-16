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

    expect(vec.findFirstEqIdx(1) == 0_i);
    expect(vec.findFirstEqIdx(2) == 1_i);
    expect(vec.findFirstEqIdx(3) == 2_i);
    expect(vec.findFirstEqIdx(4) == 3_i);
    expect(vec.findFirstEqIdx(5) == 4_i);
    expect(vec.findFirstEqIdx(6) == 5_i);
    expect(vec.findFirstEqIdx(7) == 6_i);
    expect(vec.findFirstEqIdx(8) == 7_i);
    expect(vec.findFirstEqIdx(9) == -1_i);
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

  "find_if uint8_t small vector tests"_test = [] {
    // Test with uint8_t small vector
    SVector<std::uint8_t, 16> vec{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

    // Find first element > 10
    auto idx = find_first(vec, [](auto x) { return x > 10; });
    expect(idx == 10_i); // Should find 11 at index 10

    // Find first element == 5
    idx = find_first(vec, [](auto x) { return x == 5; });
    expect(idx == 4_i); // Should find 5 at index 4

    // Find first element > 20 (doesn't exist)
    idx = find_first(vec, [](auto x) { return x > 20; });
    expect(idx == -1_i); // Should return -1

    // Test findFirstEqIdx
    expect(vec.findFirstEqIdx(1) == 0_i);
    expect(vec.findFirstEqIdx(8) == 7_i);
    expect(vec.findFirstEqIdx(16) == 15_i);
    expect(vec.findFirstEqIdx(20) == -1_i);
  };

  "find_if uint8_t medium vector tests"_test = [] {
    // Test with 32-element uint8_t vector (AVX256 path)
    SVector<std::uint8_t, 32> vec;
    for (std::ptrdiff_t i = 0; i < 32; ++i)
      vec[i] = static_cast<std::uint8_t>(i * 3); // [0, 3, 6, 9, ..., 93]

    // Find first element >= 50
    auto idx = find_first(vec, [](auto x) { return x >= 50; });
    expect(idx == 17_i); // Should find 51 at index 17

    // Find first element == 0
    idx = find_first(vec, [](auto x) { return x == 0; });
    expect(idx == 0_i); // Should find 0 at index 0

    // Find first element > 100 (doesn't exist)
    idx = find_first(vec, [](auto x) { return x > 100; });
    expect(idx == -1_i); // Should return -1
  };

  "find_if uint8_t large vector tests"_test = [] {
    // Test with 64-element uint8_t vector (AVX512 path)
    SVector<std::uint8_t, 64> vec;
    for (std::ptrdiff_t i = 0; i < 64; ++i)
      vec[i] = static_cast<std::uint8_t>(i + 100); // [100, 101, ..., 163]

    // Find first element > 150
    auto idx = find_first(vec, [](auto x) { return x > 150; });
    expect(idx == 51_i); // Should find 151 at index 51

    // Find first element == 120
    idx = find_first(vec, [](auto x) { return x == 120; });
    expect(idx == 20_i); // Should find 120 at index 20

    // Find first element < 100 (doesn't exist)
    idx = find_first(vec, [](auto x) { return x < 100; });
    expect(idx == -1_i); // Should return -1

    // Test with all matching
    SVector<std::uint8_t, 64> all_same;
    for (std::ptrdiff_t i = 0; i < 64; ++i)
      all_same[i] = 42;

    idx = find_first(all_same, [](auto x) { return x == 42; });
    expect(idx == 0_i); // Should find at first position
  };

  "find_if uint8_t very large vector tests"_test = [] {
    // Test with very large vector to ensure multi-chunk processing
    SVector<std::uint8_t, 200> vec;
    for (std::ptrdiff_t i = 0; i < 200; ++i)
      vec[i] = static_cast<std::uint8_t>(i % 256);

    // Find element at end
    auto idx = find_first(vec, [](auto x) { return x == 199; });
    expect(idx == 199_i);

    // Find element in middle
    idx = find_first(vec, [](auto x) { return x == 100; });
    expect(idx == 100_i);

    // Find non-existent (would be > 255)
    idx = find_first(vec, [](auto x) { return x == 255; });
    expect(idx == -1_i);
  };
}

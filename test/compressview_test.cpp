import boost.ut;

import Arena;
import Array;
import AxisTypes;
import ManagedArray;
import SIMD;
import std;

using namespace ::math;
using namespace boost::ut;

int main() {
  "CompressView Basic Operations"_test = [] {
    // Test basic compress construction and indexing
    auto vec = ManagedArray<std::int64_t, Length<>>{length(4)};
    for (std::ptrdiff_t i = 0; i < 4; ++i) vec[i] = i + 1; // vec = [1, 2, 3, 4]

    // Create compress with mask 0x35 = 0b00110101 (selects elements 0, 2, 4, 5)
    std::uint8_t mask = 0x35;
    auto filtered = vec[Compress(mask)];

    // Test compress store: v[Compress(0x35)] << [10, 20, 30, 40]
    auto source = ManagedArray<std::int64_t, Length<>>{length(8)};
    for (std::ptrdiff_t i = 0; i < 8; ++i)
      source[i] = 10 * (i + 1); // vec = [10, 20, 30, 40, 50, 60, 70, 80]

    filtered << source;

    // Check that compress store worked correctly
    // Original positions: 0, 2, 4, 5 should now contain [10, 20, 30, 40]
    expect(vec[0] == 10_i); // Position 0: 10
    expect(vec[1] == 30_i); // Position 1: unchanged
    expect(vec[2] == 50_i); // Position 2: 20
    expect(vec[3] == 60_i); // Position 3: unchanged
  };

  "CompressView Expression Usage"_test = [] {
    // Test basic compressed indexing
    auto x = ManagedArray<std::int64_t, Length<>>{length(8)};

    // Initialize vectors
    for (std::ptrdiff_t i = 0; i < 8; ++i)
      x[i] = i + 1; // x = [1, 2, 3, 4, 5, 6, 7, 8]

    // Compress mask 0x29 = 0b00101001 (selects positions 0, 3, 5)
    std::uint8_t mask = 0x29;

    // Create compressed view - should compress x[0], x[3], x[5] into contiguous
    // storage
    auto result = ManagedArray<std::int64_t, Length<>>{length(3)};
    auto compressed_result = result[Compress(mask)];

    compressed_result << x; // Should store x[0]=1, x[3]=4, x[5]=6 into result

    // Check compressed storage worked
    expect(result[0] == 1_i); // x[0]
    expect(result[1] == 4_i); // x[3]
    expect(result[2] == 6_i); // x[5]
  };

  "CompressView SIMD Width Handling"_test = [] {
    // Test with different vector sizes to verify SIMD and scalar paths
    static constexpr std::ptrdiff_t N = 9;

    auto vec = ManagedArray<std::int32_t, Length<>>{length(N)};
    for (std::ptrdiff_t i = 0; i < N; ++i)
      vec[i] = i * 10; // [0, 10, 20, 30, ...]

    // Use a pattern that tests SIMD alignment
    std::uint16_t mask = 0x555d; // 0101010101011101
    auto filtered = vec[Compress(mask)];

    auto new_values = ManagedArray<std::int32_t, Length<>>{length(16)};
    for (std::ptrdiff_t i = 0; i < 16; ++i)
      new_values[i] = (i + 1) * 100; // [100, 200, 300, ...]

    filtered << new_values;

    // Check alternating positions got updated
    expect(vec[0] == 100_i);
    expect(vec[1] == 300_i);
    expect(vec[2] == 400_i);
    expect(vec[3] == 500_i);
    expect(vec[4] == 700_i);
    expect(vec[5] == 900_i);
    expect(vec[6] == 1100_i);
    expect(vec[7] == 1300_i);
    expect(vec[8] == 1500_i);
  };
}

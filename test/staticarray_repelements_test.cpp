import boost.ut;
import StaticArray;
import SIMD;
import std;

using namespace boost::ut;

namespace {

void testRepElements4WithDouble() {
  // Test with double (typically Width=2 on 128-bit SIMD)
  // Use 2 elements to match natural SIMD width
  ::math::StaticArray<double, 1, 2> arr;
  arr[0] = 1.5;
  arr[1] = 2.5;

  auto result = arr.repElements4();

  // Result should be [1.5, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5]
  static_assert(decltype(result)::numCol() == 8);
  expect(eq(result[0], 1.5));
  expect(eq(result[1], 1.5));
  expect(eq(result[2], 1.5));
  expect(eq(result[3], 1.5));
  expect(eq(result[4], 2.5));
  expect(eq(result[5], 2.5));
  expect(eq(result[6], 2.5));
  expect(eq(result[7], 2.5));
}

void testRepElements4WithFloat() {
  // Test with float (typically Width=4 on 128-bit SIMD)
  // Use 4 elements to match natural SIMD width
  ::math::StaticArray<float, 1, 4> arr;
  arr[0] = 1.0f;
  arr[1] = 2.0f;
  arr[2] = 3.0f;
  arr[3] = 4.0f;

  auto result = arr.repElements4();

  // Result should be [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
  static_assert(decltype(result)::numCol() == 16);
  for (std::ptrdiff_t i = 0; i < 4; ++i) expect(eq(result[i], 1.0f));
  for (std::ptrdiff_t i = 4; i < 8; ++i) expect(eq(result[i], 2.0f));
  for (std::ptrdiff_t i = 8; i < 12; ++i) expect(eq(result[i], 3.0f));
  for (std::ptrdiff_t i = 12; i < 16; ++i) expect(eq(result[i], 4.0f));
}

void testRepElements4MultiRow() {
  // Test with multiple rows
  // Using double with 2 columns (matches Width=2)
  ::math::StaticArray<double, 2, 2> arr;
  arr[0, 0] = 1.0;
  arr[0, 1] = 2.0;
  arr[1, 0] = 3.0;
  arr[1, 1] = 4.0;

  auto result = arr.repElements4();

  // Result should have 2 rows and 8 columns
  static_assert(decltype(result)::numRow() == 2);
  static_assert(decltype(result)::numCol() == 8);

  // First row: [1, 1, 1, 1, 2, 2, 2, 2]
  for (std::ptrdiff_t i = 0; i < 4; ++i) expect(eq(result[0, i], 1.0));
  for (std::ptrdiff_t i = 4; i < 8; ++i) expect(eq(result[0, i], 2.0));

  // Second row: [3, 3, 3, 3, 4, 4, 4, 4]
  for (std::ptrdiff_t i = 0; i < 4; ++i) expect(eq(result[1, i], 3.0));
  for (std::ptrdiff_t i = 4; i < 8; ++i) expect(eq(result[1, i], 4.0));
}

void testRepElements4WithInt() {
  // Test with int (typically Width=4 on 128-bit SIMD)
  // Use 4 elements to match natural SIMD width
  ::math::StaticArray<int, 1, 4> arr;
  arr[0] = 10;
  arr[1] = 20;
  arr[2] = 30;
  arr[3] = 40;

  auto result = arr.repElements4();

  // Result should be [10, 10, 10, 10, 20, 20, 20, 20, 30, 30, 30, 30, 40, 40,
  // 40, 40]
  static_assert(decltype(result)::numCol() == 16);
  for (std::ptrdiff_t i = 0; i < 4; ++i) expect(eq(result[i], 10_i));
  for (std::ptrdiff_t i = 4; i < 8; ++i) expect(eq(result[i], 20_i));
  for (std::ptrdiff_t i = 8; i < 12; ++i) expect(eq(result[i], 30_i));
  for (std::ptrdiff_t i = 12; i < 16; ++i) expect(eq(result[i], 40_i));
}

} // namespace

auto main() -> int {
  "RepElements4 Double"_test = [] -> void { testRepElements4WithDouble(); };
  "RepElements4 Float"_test = [] -> void { testRepElements4WithFloat(); };
  "RepElements4 Int"_test = [] -> void { testRepElements4WithInt(); };
  "RepElements4 MultiRow"_test = [] -> void { testRepElements4MultiRow(); };
  return 0;
}

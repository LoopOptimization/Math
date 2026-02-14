import Testing;

import Sort;
import StaticArray;
import std;

using ::math::SVector;
using namespace testing;

auto main() -> int {
  "sort_fp64"_test = [] -> void {
    using SV16 = SVector<double, 16>;
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV16 x{};
      for (int j = 0; j < 16; ++j) x[j] = dist(rng);
      SV16 sorted = x;
      std::ranges::sort(sorted);
      expect(utils::sort(x) == sorted);
    }
  };
  "sort_fp32"_test = [] -> void {
    using SV16 = SVector<float, 16>;
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV16 x{};
      for (int j = 0; j < 16; ++j) x[j] = dist(rng);
      SV16 sorted = x;
      std::ranges::sort(sorted);
      expect(utils::sort(x) == sorted);
    }
  };
  "sort_fp64_32"_test = [] -> void {
    using SV32 = SVector<double, 32>;
    using SV16 = SVector<double, 16>;
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV32 x{};
      for (int j = 0; j < 32; ++j) x[j] = dist(rng);
      SV32 sorted = x;
      std::ranges::sort(sorted);
      expect(utils::sort(x) == sorted);

      // Test topHalf: sorted topHalf should equal first 16 of fully sorted
      SV16 top = utils::topHalf(x);
      SV16 sorted_top = utils::sort(top);
      SV16 expected_top{};
      for (int j = 0; j < 16; ++j) expected_top[j] = sorted[j];
      expect(sorted_top == expected_top);
    }
  };
  "sort_fp32_32"_test = [] -> void {
    using SV32 = SVector<float, 32>;
    using SV16 = SVector<float, 16>;
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV32 x{};
      for (int j = 0; j < 32; ++j) x[j] = dist(rng);
      SV32 sorted = x;
      std::ranges::sort(sorted);
      expect(utils::sort(x) == sorted);

      // Test topHalf: sorted topHalf should equal first 16 of fully sorted
      SV16 top = utils::topHalf(x);
      SV16 sorted_top = utils::sort(top);
      SV16 expected_top{};
      for (int j = 0; j < 16; ++j) expected_top[j] = sorted[j];
      expect(sorted_top == expected_top);
    }
  };

  "sortperm_fp64_16"_test = [] -> void {
    using SV16 = SVector<double, 16>;
    using SV16u16 = SVector<std::uint16_t, 16>;
    using SV16u64 = SVector<std::uint64_t, 16>;
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist{0.123, 1024.8};
    std::uniform_int_distribution<std::uint16_t> idist16{0, 65535};
    std::uniform_int_distribution<std::uint64_t> idist64{
      0, std::numeric_limits<std::uint64_t>::max()};

    for (int i = 0; i < 100; ++i) {
      SV16 keys{}, values{};
      SV16u16 ivalues16{};
      SV16u64 ivalues64{};
      for (int j = 0; j < 16; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
        ivalues16[j] = idist16(rng);
        ivalues64[j] = idist64(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<16>::make(keys);

      // Verify sorted_keys is correct
      SV16 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV16 permuted_values = perm(values);
      SV16u16 permuted_ivalues16 = perm(ivalues16);
      SV16u64 permuted_ivalues64 = perm(ivalues64);

      std::array<std::tuple<double, double, std::uint16_t, std::uint64_t>, 16>
        tuples;
      for (int j = 0; j < 16; ++j)
        tuples[j] = {keys[j], values[j], ivalues16[j], ivalues64[j]};
      std::ranges::sort(tuples, [](auto &a, auto &b) {
        return std::get<0>(a) < std::get<0>(b);
      });

      for (int j = 0; j < 16; ++j) {
        expect(permuted_values[j] == std::get<1>(tuples[j]));
        expect(permuted_ivalues16[j] == std::get<2>(tuples[j]));
        expect(permuted_ivalues64[j] == std::get<3>(tuples[j]));
      }
    }
  };

  "sortperm_fp32_16"_test = [] -> void {
    using SV16 = SVector<float, 16>;
    using SV16u16 = SVector<std::uint16_t, 16>;
    using SV16u64 = SVector<std::uint64_t, 16>;
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> dist{0.123f, 1024.8f};
    std::uniform_int_distribution<std::uint16_t> idist16{0, 65535};
    std::uniform_int_distribution<std::uint64_t> idist64{
      0, std::numeric_limits<std::uint64_t>::max()};

    for (int i = 0; i < 100; ++i) {
      SV16 keys{}, values{};
      SV16u16 ivalues16{};
      SV16u64 ivalues64{};
      for (int j = 0; j < 16; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
        ivalues16[j] = idist16(rng);
        ivalues64[j] = idist64(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<16>::make(keys);

      // Verify sorted_keys is correct
      SV16 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV16 permuted_values = perm(values);
      SV16u16 permuted_ivalues16 = perm(ivalues16);
      SV16u64 permuted_ivalues64 = perm(ivalues64);

      std::array<std::tuple<float, float, std::uint16_t, std::uint64_t>, 16>
        tuples;
      for (int j = 0; j < 16; ++j)
        tuples[j] = {keys[j], values[j], ivalues16[j], ivalues64[j]};
      std::ranges::sort(tuples, [](auto &a, auto &b) {
        return std::get<0>(a) < std::get<0>(b);
      });

      for (int j = 0; j < 16; ++j) {
        expect(permuted_values[j] == std::get<1>(tuples[j]));
        expect(permuted_ivalues16[j] == std::get<2>(tuples[j]));
        expect(permuted_ivalues64[j] == std::get<3>(tuples[j]));
      }
    }
  };

  "sortperm_fp64_32"_test = [] -> void {
    using SV32 = SVector<double, 32>;
    using SV16 = SVector<double, 16>;
    using SV32u16 = SVector<std::uint16_t, 32>;
    using SV16u16 = SVector<std::uint16_t, 16>;
    using SV32u64 = SVector<std::uint64_t, 32>;
    using SV16u64 = SVector<std::uint64_t, 16>;
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist{0.123, 1024.8};
    std::uniform_int_distribution<std::uint16_t> idist16{0, 65535};
    std::uniform_int_distribution<std::uint64_t> idist64{
      0, std::numeric_limits<std::uint64_t>::max()};

    for (int i = 0; i < 100; ++i) {
      SV32 keys{}, values{};
      SV32u16 ivalues16{};
      SV32u64 ivalues64{};
      for (int j = 0; j < 32; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
        ivalues16[j] = idist16(rng);
        ivalues64[j] = idist64(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<32>::make(keys);

      // Verify sorted_keys is correct
      SV32 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV32 permuted_values = perm(values);
      SV32u16 permuted_ivalues16 = perm(ivalues16);
      SV32u64 permuted_ivalues64 = perm(ivalues64);

      std::array<std::tuple<double, double, std::uint16_t, std::uint64_t>, 32>
        tuples;
      for (int j = 0; j < 32; ++j)
        tuples[j] = {keys[j], values[j], ivalues16[j], ivalues64[j]};
      std::ranges::sort(tuples, [](auto &a, auto &b) {
        return std::get<0>(a) < std::get<0>(b);
      });

      for (int j = 0; j < 32; ++j) {
        expect(permuted_values[j] == std::get<1>(tuples[j]));
        expect(permuted_ivalues16[j] == std::get<2>(tuples[j]));
        expect(permuted_ivalues64[j] == std::get<3>(tuples[j]));
      }

      // Test TopHalfPerm
      auto [thperm, top_keys] = utils::TopHalfPerm::make(keys);

      // Use SortPerm to sort the top keys and get the permutation
      auto [top_perm, sorted_top_keys] = utils::SortPerm<16>::make(top_keys);

      // Verify sorted top_keys equal first 16 of fully sorted keys
      SV16 expected_top_keys{};
      std::memcpy(&expected_top_keys, &expected_keys, sizeof(SV16));
      expect(sorted_top_keys == expected_top_keys);

      // Verify applying thperm to keys gives top_keys
      expect(thperm(keys) == top_keys);

      // Verify applying thperm to values then sorting matches key-value sort
      SV16 sorted_top_values = top_perm(thperm(values));
      SV16u16 sorted_top_ivalues16 = top_perm(thperm(ivalues16));
      SV16u64 sorted_top_ivalues64 = top_perm(thperm(ivalues64));
      for (int j = 0; j < 16; ++j) {
        expect(sorted_top_values[j] == std::get<1>(tuples[j]));
        expect(sorted_top_ivalues16[j] == std::get<2>(tuples[j]));
        expect(sorted_top_ivalues64[j] == std::get<3>(tuples[j]));
      }
    }
  };

  "sortperm_fp32_32"_test = [] -> void {
    using SV32 = SVector<float, 32>;
    using SV16 = SVector<float, 16>;
    using SV32u16 = SVector<std::uint16_t, 32>;
    using SV16u16 = SVector<std::uint16_t, 16>;
    using SV32u64 = SVector<std::uint64_t, 32>;
    using SV16u64 = SVector<std::uint64_t, 16>;
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> dist{0.123f, 1024.8f};
    std::uniform_int_distribution<std::uint16_t> idist16{0, 65535};
    std::uniform_int_distribution<std::uint64_t> idist64{
      0, std::numeric_limits<std::uint64_t>::max()};

    for (int i = 0; i < 100; ++i) {
      SV32 keys{}, values{};
      SV32u16 ivalues16{};
      SV32u64 ivalues64{};
      for (int j = 0; j < 32; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
        ivalues16[j] = idist16(rng);
        ivalues64[j] = idist64(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<32>::make(keys);

      // Verify sorted_keys is correct
      SV32 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV32 permuted_values = perm(values);
      SV32u16 permuted_ivalues16 = perm(ivalues16);
      SV32u64 permuted_ivalues64 = perm(ivalues64);

      std::array<std::tuple<float, float, std::uint16_t, std::uint64_t>, 32>
        tuples;
      for (int j = 0; j < 32; ++j)
        tuples[j] = {keys[j], values[j], ivalues16[j], ivalues64[j]};
      std::ranges::sort(tuples, [](auto &a, auto &b) {
        return std::get<0>(a) < std::get<0>(b);
      });

      for (int j = 0; j < 32; ++j) {
        expect(permuted_values[j] == std::get<1>(tuples[j]));
        expect(permuted_ivalues16[j] == std::get<2>(tuples[j]));
        expect(permuted_ivalues64[j] == std::get<3>(tuples[j]));
      }

      // Test TopHalfPerm
      auto [thperm, top_keys] = utils::TopHalfPerm::make(keys);

      // Use SortPerm to sort the top keys and get the permutation
      auto [top_perm, sorted_top_keys] = utils::SortPerm<16>::make(top_keys);

      // Verify sorted top_keys equal first 16 of fully sorted keys
      SV16 expected_top_keys{};
      std::memcpy(&expected_top_keys, &expected_keys, sizeof(SV16));
      expect(sorted_top_keys == expected_top_keys);

      // Verify applying thperm to keys gives top_keys
      expect(thperm(keys) == top_keys);

      // Verify applying thperm to values then sorting matches key-value sort
      SV16 sorted_top_values = top_perm(thperm(values));
      SV16u16 sorted_top_ivalues16 = top_perm(thperm(ivalues16));
      SV16u64 sorted_top_ivalues64 = top_perm(thperm(ivalues64));
      for (int j = 0; j < 16; ++j) {
        expect(sorted_top_values[j] == std::get<1>(tuples[j]));
        expect(sorted_top_ivalues16[j] == std::get<2>(tuples[j]));
        expect(sorted_top_ivalues64[j] == std::get<3>(tuples[j]));
      }
    }
  };

  return 0;
}

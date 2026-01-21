import boost.ut;

import Sort;
import StaticArray;
import std;

using ::math::SVector;
using namespace boost::ut;

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
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV16 keys{}, values{};
      for (int j = 0; j < 16; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<double, 16>::make(keys);

      // Verify sorted_keys is correct
      SV16 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV16 permuted_values = perm(values);

      std::array<std::pair<double, double>, 16> pairs;
      for (int j = 0; j < 16; ++j) pairs[j] = {keys[j], values[j]};
      std::ranges::sort(pairs,
                        [](auto &a, auto &b) { return a.first < b.first; });

      for (int j = 0; j < 16; ++j)
        expect(permuted_values[j] == pairs[j].second);
    }
  };

  "sortperm_fp32_16"_test = [] -> void {
    using SV16 = SVector<float, 16>;
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> dist{0.123f, 1024.8f};

    for (int i = 0; i < 100; ++i) {
      SV16 keys{}, values{};
      for (int j = 0; j < 16; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<float, 16>::make(keys);

      // Verify sorted_keys is correct
      SV16 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV16 permuted_values = perm(values);

      std::array<std::pair<float, float>, 16> pairs;
      for (int j = 0; j < 16; ++j) pairs[j] = {keys[j], values[j]};
      std::ranges::sort(pairs,
                        [](auto &a, auto &b) { return a.first < b.first; });

      for (int j = 0; j < 16; ++j)
        expect(permuted_values[j] == pairs[j].second);
    }
  };

  "sortperm_fp64_32"_test = [] -> void {
    using SV32 = SVector<double, 32>;
    using SV16 = SVector<double, 16>;
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV32 keys{}, values{};
      for (int j = 0; j < 32; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<double, 32>::make(keys);

      // Verify sorted_keys is correct
      SV32 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV32 permuted_values = perm(values);

      std::array<std::pair<double, double>, 32> pairs;
      for (int j = 0; j < 32; ++j) pairs[j] = {keys[j], values[j]};
      std::ranges::sort(pairs,
                        [](auto &a, auto &b) { return a.first < b.first; });

      for (int j = 0; j < 32; ++j)
        expect(permuted_values[j] == pairs[j].second);

      // Test TopHalfPerm
      auto [thperm, top_keys] = utils::TopHalfPerm<double>::make(keys);

      // Use SortPerm to sort the top keys and get the permutation
      auto [top_perm, sorted_top_keys] =
        utils::SortPerm<double, 16>::make(top_keys);

      // Verify sorted top_keys equal first 16 of fully sorted keys
      SV16 expected_top_keys{};
      std::memcpy(&expected_top_keys, &expected_keys, sizeof(SV16));
      expect(sorted_top_keys == expected_top_keys);

      // Verify applying thperm to keys gives top_keys
      expect(thperm(keys) == top_keys);

      // Verify applying thperm to values then sorting matches key-value sort
      SV16 sorted_top_values = top_perm(thperm(values));
      for (int j = 0; j < 16; ++j)
        expect(sorted_top_values[j] == pairs[j].second);
    }
  };

  "sortperm_fp32_32"_test = [] -> void {
    using SV32 = SVector<float, 32>;
    using SV16 = SVector<float, 16>;
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> dist{0.123f, 1024.8f};

    for (int i = 0; i < 100; ++i) {
      SV32 keys{}, values{};
      for (int j = 0; j < 32; ++j) {
        keys[j] = dist(rng);
        values[j] = dist(rng);
      }

      auto [perm, sorted_keys] = utils::SortPerm<float, 32>::make(keys);

      // Verify sorted_keys is correct
      SV32 expected_keys = keys;
      std::ranges::sort(expected_keys);
      expect(sorted_keys == expected_keys);

      // Verify applying perm to keys gives sorted_keys
      expect(perm(keys) == sorted_keys);

      // Verify applying perm to values matches key-value sort
      SV32 permuted_values = perm(values);

      std::array<std::pair<float, float>, 32> pairs;
      for (int j = 0; j < 32; ++j) pairs[j] = {keys[j], values[j]};
      std::ranges::sort(pairs,
                        [](auto &a, auto &b) { return a.first < b.first; });

      for (int j = 0; j < 32; ++j)
        expect(permuted_values[j] == pairs[j].second);

      // Test TopHalfPerm
      auto [thperm, top_keys] = utils::TopHalfPerm<float>::make(keys);

      // Use SortPerm to sort the top keys and get the permutation
      auto [top_perm, sorted_top_keys] =
        utils::SortPerm<float, 16>::make(top_keys);

      // Verify sorted top_keys equal first 16 of fully sorted keys
      SV16 expected_top_keys{};
      std::memcpy(&expected_top_keys, &expected_keys, sizeof(SV16));
      expect(sorted_top_keys == expected_top_keys);

      // Verify applying thperm to keys gives top_keys
      expect(thperm(keys) == top_keys);

      // Verify applying thperm to values then sorting matches key-value sort
      SV16 sorted_top_values = top_perm(thperm(values));
      for (int j = 0; j < 16; ++j)
        expect(sorted_top_values[j] == pairs[j].second);
    }
  };

  return 0;
}

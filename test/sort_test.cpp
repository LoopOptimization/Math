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
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV32 x{};
      for (int j = 0; j < 32; ++j) x[j] = dist(rng);
      SV32 sorted = x;
      std::ranges::sort(sorted);
      expect(utils::sort(x) == sorted);
    }
  };
  "sort_fp32_32"_test = [] -> void {
    using SV32 = SVector<float, 32>;
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> dist{0.123, 1024.8};

    for (int i = 0; i < 100; ++i) {
      SV32 x{};
      for (int j = 0; j < 32; ++j) x[j] = dist(rng);
      SV32 sorted = x;
      std::ranges::sort(sorted);
      expect(utils::sort(x) == sorted);
    }
  };
  return 0;
}

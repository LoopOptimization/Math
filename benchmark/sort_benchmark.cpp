
import Nanobench;
import Sort;
import StaticArray;
import std;

namespace {

// Helper to fill an SVector with random values
template <typename T, std::ptrdiff_t N>
auto fill_random(std::mt19937_64 &rng) -> math::SVector<T, N> {
  math::SVector<T, N> v;
  std::uniform_real_distribution<T> dist(T(-1000), T(1000));
  for (std::ptrdiff_t i = 0; i < N; ++i) v[i] = dist(rng);
  return v;
}

// Benchmark utils::sort vs std::ranges::sort
template <typename T, std::ptrdiff_t N>
void BM_sort_simd_impl(Bench &bench, const std::string &type_name) {
  std::mt19937_64 rng{42};
  static constexpr std::size_t X = 32768 / sizeof(math::SVector<T, N>);
  std::array<math::SVector<T, N>, X> xs;
  for (std::size_t i = 0; i < X; ++i) xs[i] = fill_random<T, N>(rng);

  bench.run("utils::sort<" + type_name + "," + std::to_string(N) + ">", [&] {
    for (std::size_t i = 0; i < X; ++i) {
      auto sorted = utils::sort(xs[i]);
      doNotOptimizeAway(sorted);
    }
  });
}

template <typename T, std::ptrdiff_t N>
void BM_sort_ranges_impl(Bench &bench, const std::string &type_name) {
  std::mt19937_64 rng{42};
  static constexpr std::size_t X = 32768 / sizeof(math::SVector<T, N>);
  std::array<math::SVector<T, N>, X> xs;
  for (std::size_t i = 0; i < X; ++i) xs[i] = fill_random<T, N>(rng);

  bench.run("std::ranges::sort<" + type_name + "," + std::to_string(N) + ">",
            [&] {
              for (std::size_t i = 0; i < X; ++i) {
                std::array<T, N> arr;
                std::memcpy(&arr, &xs[i], sizeof(arr));
                std::ranges::sort(arr);
                doNotOptimizeAway(arr);
              }
            });
}

// Benchmark SortPerm + apply vs sorting pairs
template <typename T, std::ptrdiff_t N>
void BM_sortperm_impl(Bench &bench, const std::string &type_name) {
  std::mt19937_64 rng{42};
  static constexpr std::size_t X = 32768 / (2 * sizeof(math::SVector<T, N>));
  std::array<math::SVector<T, N>, X> keys, values;
  for (std::size_t i = 0; i < X; ++i) {
    keys[i] = fill_random<T, N>(rng);
    values[i] = fill_random<T, N>(rng);
  }

  bench.run("SortPerm<" + type_name + "," + std::to_string(N) + ">", [&] {
    for (std::size_t i = 0; i < X; ++i) {
      auto [perm, sorted_keys] = utils::SortPerm<T, N>::make(keys[i]);
      auto sorted_values = perm(values[i]);
      doNotOptimizeAway(sorted_keys);
      doNotOptimizeAway(sorted_values);
    }
  });
}

template <typename T, std::ptrdiff_t N>
void BM_sort_pairs_impl(Bench &bench, const std::string &type_name) {
  std::mt19937_64 rng{42};
  static constexpr std::size_t X = 32768 / (2 * sizeof(math::SVector<T, N>));
  std::array<math::SVector<T, N>, X> keys, values;
  for (std::size_t i = 0; i < X; ++i) {
    keys[i] = fill_random<T, N>(rng);
    values[i] = fill_random<T, N>(rng);
  }

  bench.run("sort_pairs<" + type_name + "," + std::to_string(N) + ">", [&] {
    for (std::size_t j = 0; j < X; ++j) {
      std::array<std::pair<T, T>, std::size_t(N)> pairs;
      for (std::ptrdiff_t i = 0; i < N; ++i)
        pairs[std::size_t(i)] = {keys[j][i], values[j][i]};

      std::ranges::sort(
        pairs, [](const auto &a, const auto &b) { return a.first < b.first; });
      doNotOptimizeAway(pairs);
    }
  });
}

} // namespace

// Export benchmark functions for benchmark_main.cpp
void BM_sort_float16(Bench &bench) {
  BM_sort_simd_impl<float, 16>(bench, "float");
  BM_sort_ranges_impl<float, 16>(bench, "float");
}

void BM_sort_double16(Bench &bench) {
  BM_sort_simd_impl<double, 16>(bench, "double");
  BM_sort_ranges_impl<double, 16>(bench, "double");
}

void BM_sort_float32(Bench &bench) {
  BM_sort_simd_impl<float, 32>(bench, "float");
  BM_sort_ranges_impl<float, 32>(bench, "float");
}

void BM_sort_double32(Bench &bench) {
  BM_sort_simd_impl<double, 32>(bench, "double");
  BM_sort_ranges_impl<double, 32>(bench, "double");
}

void BM_sortperm_float16(Bench &bench) {
  BM_sortperm_impl<float, 16>(bench, "float");
  BM_sort_pairs_impl<float, 16>(bench, "float");
}

void BM_sortperm_double16(Bench &bench) {
  BM_sortperm_impl<double, 16>(bench, "double");
  BM_sort_pairs_impl<double, 16>(bench, "double");
}

void BM_sortperm_float32(Bench &bench) {
  BM_sortperm_impl<float, 32>(bench, "float");
  BM_sort_pairs_impl<float, 32>(bench, "float");
}

void BM_sortperm_double32(Bench &bench) {
  BM_sortperm_impl<double, 32>(bench, "double");
  BM_sort_pairs_impl<double, 32>(bench, "double");
}

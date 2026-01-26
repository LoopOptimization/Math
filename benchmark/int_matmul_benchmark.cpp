
import Nanobench;
import ManagedArray;
import SIMD;
import StaticArray;
import std;

namespace {

// Helper to fill an SMatrix with random values
template <std::ptrdiff_t N>
auto fill_random_matrix(std::mt19937_64 &rng)
  -> math::SMatrix<std::int64_t, N, N> {
  math::SMatrix<std::int64_t, N, N> m;
  std::uniform_int_distribution<std::int64_t> dist(-100, 100);
  for (std::ptrdiff_t i = 0; i < N; ++i)
    for (std::ptrdiff_t j = 0; j < N; ++j) m[i, j] = dist(rng);
  return m;
}

// Helper to fill a DenseMatrix with random values
void fill_random_dense_matrix(math::DenseMatrix<std::int64_t> &m,
                              std::mt19937_64 &rng) {
  std::uniform_int_distribution<std::int64_t> dist(-100, 100);
  for (std::ptrdiff_t i = 0; i < m.numRow(); ++i)
    for (std::ptrdiff_t j = 0; j < m.numCol(); ++j) m[i, j] = dist(rng);
}

// Benchmark C = A * B for static square matrices of size N
template <std::ptrdiff_t N>
void BM_matmul_impl(Bench &bench, const std::string &size_str) {
  std::mt19937_64 rng{42};
  static constexpr std::size_t X =
    32768 / sizeof(math::SMatrix<std::int64_t, N, N>);
  static constexpr std::size_t num_mats = X > 2 ? X : 2;
  std::array<math::SMatrix<std::int64_t, N, N>, num_mats> as, bs;
  for (std::size_t i = 0; i < num_mats; ++i) {
    as[i] = fill_random_matrix<N>(rng);
    bs[i] = fill_random_matrix<N>(rng);
  }

  bench.run("static_matmul<int64," + size_str + "x" + size_str + ">", [&] {
    for (std::size_t i = 0; i < num_mats; ++i) {
      math::SMatrix<std::int64_t, N, N> c;
      c << as[i] * bs[i];
      doNotOptimizeAway(c);
    }
  });
}

// Benchmark D = A * B + C for static square matrices of size N
template <std::ptrdiff_t N>
void BM_matmuladd_impl(Bench &bench, const std::string &size_str) {
  std::mt19937_64 rng{42};
  static constexpr std::size_t X =
    32768 / sizeof(math::SMatrix<std::int64_t, N, N>);
  static constexpr std::size_t num_mats = X > 2 ? X : 2;
  std::array<math::SMatrix<std::int64_t, N, N>, num_mats> as, bs, cs;
  for (std::size_t i = 0; i < num_mats; ++i) {
    as[i] = fill_random_matrix<N>(rng);
    bs[i] = fill_random_matrix<N>(rng);
    cs[i] = fill_random_matrix<N>(rng);
  }

  bench.run("static_matmuladd<int64," + size_str + "x" + size_str + ">", [&] {
    for (std::size_t i = 0; i < num_mats; ++i) {
      math::SMatrix<std::int64_t, N, N> d;
      d << as[i] * bs[i] + cs[i];
      doNotOptimizeAway(d);
    }
  });
}

// Benchmark C = A * B for dynamic DenseMatrix of size N
void BM_dense_matmul_impl(Bench &bench, std::ptrdiff_t N) {
  std::mt19937_64 rng{42};
  std::size_t mat_size = std::size_t(N * N) * sizeof(std::int64_t);
  std::size_t num_mats = 32768 / mat_size;
  if (num_mats < 2) num_mats = 2;

  math::DenseDims dim{math::row(N), math::col(N)};
  std::vector<math::DenseMatrix<std::int64_t>> as(num_mats), bs(num_mats);
  for (std::size_t i = 0; i < num_mats; ++i) {
    as[i].resize(dim);
    bs[i].resize(dim);
    fill_random_dense_matrix(as[i], rng);
    fill_random_dense_matrix(bs[i], rng);
  }

  math::DenseMatrix<std::int64_t> c{dim};
  bench.run("dense_matmul<int64," + std::to_string(N) + "x" +
              std::to_string(N) + ">",
            [&] {
              for (std::size_t i = 0; i < num_mats; ++i) {
                c << as[i] * bs[i];
                doNotOptimizeAway(c);
              }
            });
}

// Benchmark D = A * B + C for dynamic DenseMatrix of size N
void BM_dense_matmuladd_impl(Bench &bench, std::ptrdiff_t N) {
  std::mt19937_64 rng{42};
  std::size_t mat_size = std::size_t(N * N) * sizeof(std::int64_t);
  std::size_t num_mats = 32768 / mat_size;
  if (num_mats < 2) num_mats = 2;

  math::DenseDims dim{math::row(N), math::col(N)};
  std::vector<math::DenseMatrix<std::int64_t>> as(num_mats), bs(num_mats),
    cs(num_mats);
  for (std::size_t i = 0; i < num_mats; ++i) {
    as[i].resize(dim);
    bs[i].resize(dim);
    cs[i].resize(dim);
    fill_random_dense_matrix(as[i], rng);
    fill_random_dense_matrix(bs[i], rng);
    fill_random_dense_matrix(cs[i], rng);
  }

  math::DenseMatrix<std::int64_t> d{dim};
  bench.run("dense_matmuladd<int64," + std::to_string(N) + "x" +
              std::to_string(N) + ">",
            [&] {
              for (std::size_t i = 0; i < num_mats; ++i) {
                d << as[i] * bs[i] + cs[i];
                doNotOptimizeAway(d);
              }
            });
}

} // namespace

// Export benchmark functions for benchmark_main.cpp - static matrices
void BM_int_matmul_2(Bench &bench) {
  BM_matmul_impl<2>(bench, "2");
  BM_matmuladd_impl<2>(bench, "2");
}

void BM_int_matmul_3(Bench &bench) {
  BM_matmul_impl<3>(bench, "3");
  BM_matmuladd_impl<3>(bench, "3");
}

void BM_int_matmul_4(Bench &bench) {
  BM_matmul_impl<4>(bench, "4");
  BM_matmuladd_impl<4>(bench, "4");
}

void BM_int_matmul_5(Bench &bench) {
  BM_matmul_impl<5>(bench, "5");
  BM_matmuladd_impl<5>(bench, "5");
}

void BM_int_matmul_6(Bench &bench) {
  BM_matmul_impl<6>(bench, "6");
  BM_matmuladd_impl<6>(bench, "6");
}

void BM_int_matmul_7(Bench &bench) {
  BM_matmul_impl<7>(bench, "7");
  BM_matmuladd_impl<7>(bench, "7");
}

void BM_int_matmul_8(Bench &bench) {
  BM_matmul_impl<8>(bench, "8");
  BM_matmuladd_impl<8>(bench, "8");
}

void BM_int_matmul_9(Bench &bench) {
  BM_matmul_impl<9>(bench, "9");
  BM_matmuladd_impl<9>(bench, "9");
}

void BM_int_matmul_10(Bench &bench) {
  BM_matmul_impl<10>(bench, "10");
  BM_matmuladd_impl<10>(bench, "10");
}

void BM_int_matmul_11(Bench &bench) {
  BM_matmul_impl<11>(bench, "11");
  BM_matmuladd_impl<11>(bench, "11");
}

void BM_int_matmul_12(Bench &bench) {
  BM_matmul_impl<12>(bench, "12");
  BM_matmuladd_impl<12>(bench, "12");
}

void BM_int_matmul_13(Bench &bench) {
  BM_matmul_impl<13>(bench, "13");
  BM_matmuladd_impl<13>(bench, "13");
}

void BM_int_matmul_14(Bench &bench) {
  BM_matmul_impl<14>(bench, "14");
  BM_matmuladd_impl<14>(bench, "14");
}

void BM_int_matmul_15(Bench &bench) {
  BM_matmul_impl<15>(bench, "15");
  BM_matmuladd_impl<15>(bench, "15");
}

void BM_int_matmul_16(Bench &bench) {
  BM_matmul_impl<16>(bench, "16");
  BM_matmuladd_impl<16>(bench, "16");
}

// Export benchmark functions for benchmark_main.cpp - dense matrices
void BM_dense_int_matmul(Bench &bench, std::ptrdiff_t size) {
  BM_dense_matmul_impl(bench, size);
  BM_dense_matmuladd_impl(bench, size);
}

// SIMD-optimized matmul kernel using wideload/widestore
void simd_matmul_kernel(std::int64_t *__restrict c_data,
                        const std::int64_t *__restrict a_data,
                        const std::int64_t *__restrict b_data,
                        std::ptrdiff_t N) {
  using namespace simd;
  constexpr std::ptrdiff_t L = 16;
  using V = SVec<std::int64_t, L>;

  // Process 4 output rows at a time
  std::ptrdiff_t i = 0;
  for (; i + 4 <= N; i += 4) {
    V c0{}, c1{}, c2{}, c3{};
    const std::int64_t *a0 = a_data + i * N;
    const std::int64_t *a1 = a_data + (i + 1) * N;
    const std::int64_t *a2 = a_data + (i + 2) * N;
    const std::int64_t *a3 = a_data + (i + 3) * N;

    for (std::ptrdiff_t k = 0; k < N; ++k) {
      auto b_row = wideload<L>(b_data + k * N, N);
      c0 = c0 + V::vbroadcast(a0[k]) * b_row;
      c1 = c1 + V::vbroadcast(a1[k]) * b_row;
      c2 = c2 + V::vbroadcast(a2[k]) * b_row;
      c3 = c3 + V::vbroadcast(a3[k]) * b_row;
    }

    widestore<L>(c_data + i * N, c0, N);
    widestore<L>(c_data + (i + 1) * N, c1, N);
    widestore<L>(c_data + (i + 2) * N, c2, N);
    widestore<L>(c_data + (i + 3) * N, c3, N);
  }

  // Handle remaining 0-3 rows
  switch (N - i) {
  case 3: {
    V c0{}, c1{}, c2{};
    const std::int64_t *a0 = a_data + i * N;
    const std::int64_t *a1 = a_data + (i + 1) * N;
    const std::int64_t *a2 = a_data + (i + 2) * N;
    for (std::ptrdiff_t k = 0; k < N; ++k) {
      auto b_row = wideload<L>(b_data + k * N, N);
      c0 = c0 + V::vbroadcast(a0[k]) * b_row;
      c1 = c1 + V::vbroadcast(a1[k]) * b_row;
      c2 = c2 + V::vbroadcast(a2[k]) * b_row;
    }
    widestore<L>(c_data + i * N, c0, N);
    widestore<L>(c_data + (i + 1) * N, c1, N);
    widestore<L>(c_data + (i + 2) * N, c2, N);
    break;
  }
  case 2: {
    V c0{}, c1{};
    const std::int64_t *a0 = a_data + i * N;
    const std::int64_t *a1 = a_data + (i + 1) * N;
    for (std::ptrdiff_t k = 0; k < N; ++k) {
      auto b_row = wideload<L>(b_data + k * N, N);
      c0 = c0 + V::vbroadcast(a0[k]) * b_row;
      c1 = c1 + V::vbroadcast(a1[k]) * b_row;
    }
    widestore<L>(c_data + i * N, c0, N);
    widestore<L>(c_data + (i + 1) * N, c1, N);
    break;
  }
  case 1: {
    V c0{};
    const std::int64_t *a0 = a_data + i * N;
    for (std::ptrdiff_t k = 0; k < N; ++k) {
      auto b_row = wideload<L>(b_data + k * N, N);
      c0 = c0 + V::vbroadcast(a0[k]) * b_row;
    }
    widestore<L>(c_data + i * N, c0, N);
    break;
  }
  default: break;
  }
}

// Benchmark for SIMD manual matmul
void BM_simd_manual_matmul_impl(Bench &bench, std::ptrdiff_t N) {
  std::mt19937_64 rng{42};
  std::size_t mat_size = std::size_t(N * N) * sizeof(std::int64_t);
  std::size_t num_mats = 32768 / mat_size;
  if (num_mats < 2) num_mats = 2;

  math::DenseDims dim{math::row(N), math::col(N)};
  std::vector<math::DenseMatrix<std::int64_t>> as(num_mats), bs(num_mats);
  for (std::size_t i = 0; i < num_mats; ++i) {
    as[i].resize(dim);
    bs[i].resize(dim);
    fill_random_dense_matrix(as[i], rng);
    fill_random_dense_matrix(bs[i], rng);
  }

  math::DenseMatrix<std::int64_t> c{dim};
  bench.run("simd_manual_matmul<int64," + std::to_string(N) + "x" +
              std::to_string(N) + ">",
            [&] {
              for (std::size_t i = 0; i < num_mats; ++i) {
                simd_matmul_kernel(c.data(), as[i].data(), bs[i].data(), N);
                doNotOptimizeAway(c);
              }
            });
}

// Export benchmark function for benchmark_main.cpp
void BM_simd_manual_matmul(Bench &bench, std::ptrdiff_t size) {
  BM_simd_manual_matmul_impl(bench, size);
}

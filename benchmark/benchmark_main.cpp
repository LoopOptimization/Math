import Nanobench;
import std;

// Dual benchmark functions
template <std::ptrdiff_t M, std::ptrdiff_t N> void BM_dualprod(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_manual(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_simdarray(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_manual_tuple(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_simdarray_tuple(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualdivsum(Bench &bench, std::ptrdiff_t len);

// Simplex benchmark functions
void BM_Simplex0(Bench &bench);
void BM_Simplex1(Bench &bench);

// Integer matrix benchmark functions
void BM_solve_system(Bench& bench, std::ptrdiff_t size);

// Tensor benchmark functions
void BM_dual8x2dApI(Bench& bench, std::ptrdiff_t size);
void BM_dual8x2BmApI(Bench& bench, std::ptrdiff_t size);
void BM_dual8x2BmApI_manual(Bench& bench, std::ptrdiff_t size);
void BM_dual7x2dApI(Bench& bench, std::ptrdiff_t size);
void BM_dual7x2BmApI(Bench& bench, std::ptrdiff_t size);
void BM_dual7x2BmApI_manual(Bench& bench, std::ptrdiff_t size);
void BM_dual6x2dApI(Bench& bench, std::ptrdiff_t size);
void BM_dual6x2BmApI(Bench& bench, std::ptrdiff_t size);
void BM_dual6x2BmApI_manual(Bench& bench, std::ptrdiff_t size);

// Expm benchmark functions
void BM_expm_dual1(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual3(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual4(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual5(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual6(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual7(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual8(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual1x2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual2x2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual3x2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual4x2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual5x2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual6x2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual7x2(Bench& bench, std::ptrdiff_t size);
void BM_expm_dual8x2(Bench& bench, std::ptrdiff_t size);

int main() {
  Bench bench;

  // Dual product benchmarks
  bench.title("Dual Product Benchmarks");
  bench.run("BM_dualprod<6,2>", [&] { BM_dualprod<6, 2>(bench); });
  bench.run("BM_dualprod<7,2>", [&] { BM_dualprod<7, 2>(bench); });
  bench.run("BM_dualprod<8,2>", [&] { BM_dualprod<8, 2>(bench); });

  bench.run("BM_dualprod_manual<6,2>",
            [&] { BM_dualprod_manual<6, 2>(bench); });
  bench.run("BM_dualprod_manual<7,2>",
            [&] { BM_dualprod_manual<7, 2>(bench); });
  bench.run("BM_dualprod_manual<8,2>",
            [&] { BM_dualprod_manual<8, 2>(bench); });

  bench.run("BM_dualprod_simdarray<6,2>",
            [&] { BM_dualprod_simdarray<6, 2>(bench); });
  bench.run("BM_dualprod_simdarray<7,2>",
            [&] { BM_dualprod_simdarray<7, 2>(bench); });
  bench.run("BM_dualprod_simdarray<8,2>",
            [&] { BM_dualprod_simdarray<8, 2>(bench); });

  bench.run("BM_dualprod_manual_tuple<6,2>",
            [&] { BM_dualprod_manual_tuple<6, 2>(bench); });
  bench.run("BM_dualprod_manual_tuple<7,2>",
            [&] { BM_dualprod_manual_tuple<7, 2>(bench); });
  bench.run("BM_dualprod_manual_tuple<8,2>",
            [&] { BM_dualprod_manual_tuple<8, 2>(bench); });

  bench.run("BM_dualprod_simdarray_tuple<6,2>",
            [&] { BM_dualprod_simdarray_tuple<6, 2>(bench); });
  bench.run("BM_dualprod_simdarray_tuple<7,2>",
            [&] { BM_dualprod_simdarray_tuple<7, 2>(bench); });
  bench.run("BM_dualprod_simdarray_tuple<8,2>",
            [&] { BM_dualprod_simdarray_tuple<8, 2>(bench); });

  // Dual division sum benchmarks with different sizes
  bench.title("Dual Division Sum Benchmarks");
  for (std::ptrdiff_t len = 1; len <= 1024; len *= 2) {
    bench.run("BM_dualdivsum<7,0>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<7, 0>(bench, len); });
    bench.run("BM_dualdivsum<8,0>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<8, 0>(bench, len); });
    bench.run("BM_dualdivsum<7,2>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<7, 2>(bench, len); });
    bench.run("BM_dualdivsum<8,2>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<8, 2>(bench, len); });
    bench.run("BM_dualdivsum<7,4>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<7, 4>(bench, len); });
    bench.run("BM_dualdivsum<8,4>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<8, 4>(bench, len); });
  }

  // Simplex benchmarks
  bench.title("Simplex Benchmarks");
  bench.run("BM_Simplex0", [&] { BM_Simplex0(bench); });
  bench.run("BM_Simplex1", [&] { BM_Simplex1(bench); });

  // Integer matrix benchmarks
  bench.title("Integer Matrix Benchmarks");
  for (std::ptrdiff_t size = 2; size <= 10; ++size) {
    bench.run("BM_solve_system_size=" + std::to_string(size), [&] { BM_solve_system(bench, size); });
  }

  // Tensor benchmarks
  bench.title("Tensor Benchmarks");
  for (std::ptrdiff_t size = 2; size <= 10; ++size) {
    bench.run("BM_dual8x2dApI_size=" + std::to_string(size), [&] { BM_dual8x2dApI(bench, size); });
    bench.run("BM_dual8x2BmApI_size=" + std::to_string(size), [&] { BM_dual8x2BmApI(bench, size); });
    bench.run("BM_dual8x2BmApI_manual_size=" + std::to_string(size), [&] { BM_dual8x2BmApI_manual(bench, size); });
    bench.run("BM_dual7x2dApI_size=" + std::to_string(size), [&] { BM_dual7x2dApI(bench, size); });
    bench.run("BM_dual7x2BmApI_size=" + std::to_string(size), [&] { BM_dual7x2BmApI(bench, size); });
    bench.run("BM_dual7x2BmApI_manual_size=" + std::to_string(size), [&] { BM_dual7x2BmApI_manual(bench, size); });
    bench.run("BM_dual6x2dApI_size=" + std::to_string(size), [&] { BM_dual6x2dApI(bench, size); });
    bench.run("BM_dual6x2BmApI_size=" + std::to_string(size), [&] { BM_dual6x2BmApI(bench, size); });
    bench.run("BM_dual6x2BmApI_manual_size=" + std::to_string(size), [&] { BM_dual6x2BmApI_manual(bench, size); });
  }

  // Expm benchmarks
  bench.title("Expm Benchmarks");
  for (std::ptrdiff_t size = 2; size <= 10; ++size) {
    bench.run("BM_expm_dual1_size=" + std::to_string(size), [&] { BM_expm_dual1(bench, size); });
    bench.run("BM_expm_dual2_size=" + std::to_string(size), [&] { BM_expm_dual2(bench, size); });
    bench.run("BM_expm_dual3_size=" + std::to_string(size), [&] { BM_expm_dual3(bench, size); });
    bench.run("BM_expm_dual4_size=" + std::to_string(size), [&] { BM_expm_dual4(bench, size); });
    bench.run("BM_expm_dual5_size=" + std::to_string(size), [&] { BM_expm_dual5(bench, size); });
    bench.run("BM_expm_dual6_size=" + std::to_string(size), [&] { BM_expm_dual6(bench, size); });
    bench.run("BM_expm_dual7_size=" + std::to_string(size), [&] { BM_expm_dual7(bench, size); });
    bench.run("BM_expm_dual8_size=" + std::to_string(size), [&] { BM_expm_dual8(bench, size); });
    bench.run("BM_expm_dual1x2_size=" + std::to_string(size), [&] { BM_expm_dual1x2(bench, size); });
    bench.run("BM_expm_dual2x2_size=" + std::to_string(size), [&] { BM_expm_dual2x2(bench, size); });
    bench.run("BM_expm_dual3x2_size=" + std::to_string(size), [&] { BM_expm_dual3x2(bench, size); });
    bench.run("BM_expm_dual4x2_size=" + std::to_string(size), [&] { BM_expm_dual4x2(bench, size); });
    bench.run("BM_expm_dual5x2_size=" + std::to_string(size), [&] { BM_expm_dual5x2(bench, size); });
    bench.run("BM_expm_dual6x2_size=" + std::to_string(size), [&] { BM_expm_dual6x2(bench, size); });
    bench.run("BM_expm_dual7x2_size=" + std::to_string(size), [&] { BM_expm_dual7x2(bench, size); });
    bench.run("BM_expm_dual8x2_size=" + std::to_string(size), [&] { BM_expm_dual8x2(bench, size); });
  }

  return 0;
}

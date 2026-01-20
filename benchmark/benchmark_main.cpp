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
void BM_solve_system(Bench &bench, std::ptrdiff_t size);

// Tensor benchmark functions
void BM_dual8x2dApI(Bench &bench, std::ptrdiff_t size);
void BM_dual8x2BmApI(Bench &bench, std::ptrdiff_t size);
void BM_dual8x2BmApI_manual(Bench &bench, std::ptrdiff_t size);
void BM_dual7x2dApI(Bench &bench, std::ptrdiff_t size);
void BM_dual7x2BmApI(Bench &bench, std::ptrdiff_t size);
void BM_dual7x2BmApI_manual(Bench &bench, std::ptrdiff_t size);
void BM_dual6x2dApI(Bench &bench, std::ptrdiff_t size);
void BM_dual6x2BmApI(Bench &bench, std::ptrdiff_t size);
void BM_dual6x2BmApI_manual(Bench &bench, std::ptrdiff_t size);

// Expm benchmark functions
void BM_expm_dual1(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual3(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual4(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual5(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual6(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual7(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual8(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual1x2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual2x2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual3x2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual4x2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual5x2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual6x2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual7x2(Bench &bench, std::ptrdiff_t size);
void BM_expm_dual8x2(Bench &bench, std::ptrdiff_t size);

// Sort benchmark functions
void BM_sort_float16(Bench &bench);
void BM_sort_double16(Bench &bench);
void BM_sort_float32(Bench &bench);
void BM_sort_double32(Bench &bench);
void BM_sortperm_float16(Bench &bench);
void BM_sortperm_double16(Bench &bench);
void BM_sortperm_float32(Bench &bench);
void BM_sortperm_double32(Bench &bench);

int main() {
  Bench bench;
  // bench.minEpochTime();

  // // Dual product benchmarks
  // bench.title("Dual Product Benchmarks");
  // BM_dualprod<6, 2>(bench);
  // BM_dualprod<7, 2>(bench);
  // BM_dualprod<8, 2>(bench);

  // BM_dualprod_manual<6, 2>(bench);
  // BM_dualprod_manual<7, 2>(bench);
  // BM_dualprod_manual<8, 2>(bench);

  // BM_dualprod_simdarray<6, 2>(bench);
  // BM_dualprod_simdarray<7, 2>(bench);
  // BM_dualprod_simdarray<8, 2>(bench);

  // BM_dualprod_manual_tuple<6, 2>(bench);
  // BM_dualprod_manual_tuple<7, 2>(bench);
  // BM_dualprod_manual_tuple<8, 2>(bench);

  // BM_dualprod_simdarray_tuple<6, 2>(bench);
  // BM_dualprod_simdarray_tuple<7, 2>(bench);
  // BM_dualprod_simdarray_tuple<8, 2>(bench);

  // // Dual division sum benchmarks with different sizes
  // bench.title("Dual Division Sum Benchmarks");
  // for (std::ptrdiff_t len = 1; len <= 1024; len *= 2) {
  //   BM_dualdivsum<7, 0>(bench, len);
  //   BM_dualdivsum<8, 0>(bench, len);
  //   BM_dualdivsum<7, 2>(bench, len);
  //   BM_dualdivsum<8, 2>(bench, len);
  //   BM_dualdivsum<7, 4>(bench, len);
  //   BM_dualdivsum<8, 4>(bench, len);
  // }

  // // Simplex benchmarks
  // bench.title("Simplex Benchmarks");
  // BM_Simplex0(bench);
  // BM_Simplex1(bench);

  // // Integer matrix benchmarks
  // bench.title("Integer Matrix Benchmarks");
  // for (std::ptrdiff_t size = 2; size <= 10; ++size)
  //   BM_solve_system(bench, size);

  // // Tensor benchmarks
  // bench.title("Tensor Benchmarks");
  // for (std::ptrdiff_t size = 2; size <= 10; ++size) {
  //   BM_dual8x2dApI(bench, size);
  //   BM_dual8x2BmApI(bench, size);
  //   BM_dual8x2BmApI_manual(bench, size);
  //   BM_dual7x2dApI(bench, size);
  //   BM_dual7x2BmApI(bench, size);
  //   BM_dual7x2BmApI_manual(bench, size);
  //   BM_dual6x2dApI(bench, size);
  //   BM_dual6x2BmApI(bench, size);
  //   BM_dual6x2BmApI_manual(bench, size);
  // }

  // // Expm benchmarks
  // bench.title("Expm Benchmarks");
  // for (std::ptrdiff_t size = 2; size <= 10; ++size) {
  //   BM_expm_dual1(bench, size);
  //   BM_expm_dual2(bench, size);
  //   BM_expm_dual3(bench, size);
  //   BM_expm_dual4(bench, size);
  //   BM_expm_dual5(bench, size);
  //   BM_expm_dual6(bench, size);
  //   BM_expm_dual7(bench, size);
  //   BM_expm_dual8(bench, size);
  //   BM_expm_dual1x2(bench, size);
  //   BM_expm_dual2x2(bench, size);
  //   BM_expm_dual3x2(bench, size);
  //   BM_expm_dual4x2(bench, size);
  //   BM_expm_dual5x2(bench, size);
  //   BM_expm_dual6x2(bench, size);
  //   BM_expm_dual7x2(bench, size);
  //   BM_expm_dual8x2(bench, size);
  // }

  // Sort benchmarks
  bench.title("Sort Benchmarks");
  BM_sort_float16(bench);
  BM_sort_double16(bench);
  BM_sort_float32(bench);
  BM_sort_double32(bench);

  // SortPerm benchmarks
  bench.title("SortPerm Benchmarks");
  BM_sortperm_float16(bench);
  BM_sortperm_double16(bench);
  BM_sortperm_float32(bench);
  BM_sortperm_double32(bench);

  return 0;
}

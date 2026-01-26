import Nanobench;
import std;

// Simplex benchmark functions
void BM_Simplex0(Bench &bench);
void BM_Simplex1(Bench &bench);

// Integer matrix benchmark functions
void BM_solve_system(Bench &bench, std::ptrdiff_t size);

// Integer matrix multiplication benchmark functions
void BM_int_matmul_2(Bench &bench);
void BM_int_matmul_3(Bench &bench);
void BM_int_matmul_4(Bench &bench);
void BM_int_matmul_5(Bench &bench);
void BM_int_matmul_6(Bench &bench);
void BM_int_matmul_7(Bench &bench);
void BM_int_matmul_8(Bench &bench);
void BM_int_matmul_9(Bench &bench);
void BM_int_matmul_10(Bench &bench);
void BM_int_matmul_11(Bench &bench);
void BM_int_matmul_12(Bench &bench);
void BM_int_matmul_13(Bench &bench);
void BM_int_matmul_14(Bench &bench);
void BM_int_matmul_15(Bench &bench);
void BM_int_matmul_16(Bench &bench);

// Dense integer matrix multiplication benchmark functions
void BM_dense_int_matmul(Bench &bench, std::ptrdiff_t size);

// SIMD manual integer matrix multiplication benchmark functions
void BM_simd_manual_matmul(Bench &bench, std::ptrdiff_t size);

// Sort benchmark functions
void BM_sort_float16(Bench &bench);
void BM_sort_double16(Bench &bench);
void BM_sort_float32(Bench &bench);
void BM_sort_double32(Bench &bench);
void BM_sortperm_float16(Bench &bench);
void BM_sortperm_double16(Bench &bench);
void BM_sortperm_float32(Bench &bench);
void BM_sortperm_double32(Bench &bench);

auto main() -> int {
  Bench bench;

  // Simplex benchmarks
  bench.title("Simplex Benchmarks");
  BM_Simplex0(bench);
  BM_Simplex1(bench);

  // Integer matrix benchmarks
  bench.title("Integer Matrix Benchmarks");
  for (std::ptrdiff_t size = 2; size <= 10; ++size)
    BM_solve_system(bench, size);

  // Integer matrix multiplication benchmarks
  bench.title("Integer Matrix Multiplication Benchmarks");
  BM_int_matmul_2(bench);
  BM_int_matmul_3(bench);
  BM_int_matmul_4(bench);
  BM_int_matmul_5(bench);
  BM_int_matmul_6(bench);
  BM_int_matmul_7(bench);
  BM_int_matmul_8(bench);
  BM_int_matmul_9(bench);
  BM_int_matmul_10(bench);
  BM_int_matmul_11(bench);
  BM_int_matmul_12(bench);
  BM_int_matmul_13(bench);
  BM_int_matmul_14(bench);
  BM_int_matmul_15(bench);
  BM_int_matmul_16(bench);

  // Dense integer matrix multiplication benchmarks
  bench.title("Dense Integer Matrix Multiplication Benchmarks");
  for (std::ptrdiff_t size = 2; size <= 16; ++size)
    BM_dense_int_matmul(bench, size);

  // Manual SIMD integer matrix multiplication benchmarks
  bench.title("Manual SIMD Integer Matrix Multiplication Benchmarks");
  for (std::ptrdiff_t size = 2; size <= 16; ++size)
    BM_simd_manual_matmul(bench, size);

  // Sort benchmarks
  bench.title("Sort Benchmarks").warmup(1024).epochIterations(2048);
  BM_sort_float16(bench);
  BM_sort_double16(bench);
  BM_sort_float32(bench);
  BM_sort_double32(bench);

  // SortPerm benchmarks
  bench.title("SortPerm Benchmarks").warmup(1024).epochIterations(2048);
  BM_sortperm_float16(bench);
  BM_sortperm_double16(bench);
  BM_sortperm_float32(bench);
  BM_sortperm_double32(bench);

  return 0;
}

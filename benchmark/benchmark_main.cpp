import Nanobench;
import std;

// Simplex benchmark functions
void BM_Simplex0(Bench &bench);
void BM_Simplex1(Bench &bench);

// Integer matrix benchmark functions
void BM_solve_system(Bench &bench, std::ptrdiff_t size);

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

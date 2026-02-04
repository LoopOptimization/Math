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

// GCD benchmark functions
void init_gcd_data();
void BM_gcd_scalar_int64(Bench &bench);
void BM_gcd_svec_int64(Bench &bench);
void BM_gcd_svec_double(Bench &bench);
void BM_gcdreduce_svec_int64(Bench &bench);
void BM_dgcdx_scalar_int64(Bench &bench);
void BM_dgcdx_vec_int64(Bench &bench);
void BM_dgcdx_vec_double(Bench &bench);

auto main(int argc, char *argv[]) -> int {
  Bench bench;
  std::string_view filter = (argc == 2) ? argv[1] : "";
  auto should_run = [&](std::string_view title) {
    return filter.empty() || filter == title;
  };

  // Simplex benchmarks
  if (should_run("Simplex Benchmarks")) {
    bench.title("Simplex Benchmarks");
    BM_Simplex0(bench);
    BM_Simplex1(bench);
  }

  // Integer matrix benchmarks
  if (should_run("Integer Matrix Benchmarks")) {
    bench.title("Integer Matrix Benchmarks");
    for (std::ptrdiff_t size = 2; size <= 10; ++size)
      BM_solve_system(bench, size);
  }

  // Sort benchmarks
  if (should_run("Sort Benchmarks")) {
    bench.title("Sort Benchmarks").warmup(1024).epochIterations(2048);
    BM_sort_float16(bench);
    BM_sort_double16(bench);
    BM_sort_float32(bench);
    BM_sort_double32(bench);
  }

  // SortPerm benchmarks
  if (should_run("SortPerm Benchmarks")) {
    bench.title("SortPerm Benchmarks").warmup(1024).epochIterations(2048);
    BM_sortperm_float16(bench);
    BM_sortperm_double16(bench);
    BM_sortperm_float32(bench);
    BM_sortperm_double32(bench);
  }

  // GCD benchmarks
  if (should_run("GCD Benchmarks")) {
    init_gcd_data();
    bench.title("GCD Benchmarks").warmup(1024).epochIterations(2048);
    BM_gcd_scalar_int64(bench);
    BM_gcd_svec_int64(bench);
    BM_gcd_svec_double(bench);
    BM_gcdreduce_svec_int64(bench);
    BM_dgcdx_scalar_int64(bench);
    BM_dgcdx_vec_int64(bench);
    BM_dgcdx_vec_double(bench);
  }

  return 0;
}

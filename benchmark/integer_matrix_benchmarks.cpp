import ManagedArray;
import ArrayParse;
import Nanobench;
import NormalForm;
import std;

static void BM_solve_system(benchmark::State &state) {
  using math::_;
  std::mt19937_64 rng0;
  std::uniform_int_distribution<> distrib(-10, 10);
  std::ptrdiff_t d = state.range(0), num_iter = 100;
  math::DenseDims dim{math::DenseDims<>{math::row(num_iter * d), math::col(d)}};
  math::DenseMatrix<std::int64_t> A{dim}, B{dim}, C{dim}, D{dim};
  for (auto &&x : C) x = distrib(rng0);
  for (auto &&x : D) x = distrib(rng0);
  for (auto b : state) {
    A << C;
    B << D;
    for (std::ptrdiff_t n = 0; n < num_iter; ++n)
      math::NormalForm::solveSystem(A[_(n * d, (n + 1) * d), _],
                                    B[_(n * d, (n + 1) * d), _]);
  }
}
BENCHMARK(BM_solve_system)->DenseRange(2, 10, 1);

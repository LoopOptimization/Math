import Nanobench;
import std;
import ManagedArray;
import ArrayParse;
import Nanobench;
import NormalForm;
import std;

void BM_solve_system(Bench &bench, std::ptrdiff_t size) {
  using math::_;
  std::mt19937_64 rng0;
  std::uniform_int_distribution<> distrib(-10, 10);
  std::ptrdiff_t d = size, num_iter = 100;
  math::DenseDims dim{math::DenseDims<>{math::row(num_iter * d), math::col(d)}};
  math::DenseMatrix<std::int64_t> A{dim}, B{dim}, C{dim}, D{dim};
  for (auto &&x : C) x = distrib(rng0);
  for (auto &&x : D) x = distrib(rng0);
  bench.run("BM_solve_system_size=" + std::to_string(size), [&] {
    A << C;
    B << D;
    for (std::ptrdiff_t n = 0; n < num_iter; ++n)
      math::NormalForm::solveSystem(A[_(n * d, (n + 1) * d), _],
                                    B[_(n * d, (n + 1) * d), _]);
  });
}

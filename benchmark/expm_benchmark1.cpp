import ExpMat;
import Nanobench;
import std;

void BM_expm_dual1(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<double, 1>;
  SquareMatrix<D> A{SquareDims{math::row(size)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  bench.run("BM_expm_dual1_size=" + std::to_string(size), [&] { expbench(A); });
}

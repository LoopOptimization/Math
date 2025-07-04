import ExpMat;
import Nanobench;
import std;

static void BM_expm_dual5x2(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 5>, 2>;
  SquareMatrix<D> A{SquareDims{math::row(size)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  bench.run([&] { expbench(A); });
}

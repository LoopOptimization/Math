import ExpMat;
import Nanobench;
import std;

void BM_expm(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  SquareMatrix<double> A{SquareDims{math::row(size)}};
  for (auto &&a : A) a = URand<double>{}(rng0);
  bench.run("BM_expm_size=" + std::to_string(size), [&] { expbench(A); });
}

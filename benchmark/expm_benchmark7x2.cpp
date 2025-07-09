import ExpMat;
import Nanobench;
import std;

void BM_expm_dual7x2(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 7>, 2>;
  SquareMatrix<D> A{SquareDims{math::row(size)}};
#ifndef NDEBUG
  static_assert(
    std::same_as<math::ElementwiseBinaryOp<math::Array<D, SquareDims<>, true>,
                                           math::Array<D, SquareDims<>, true>,
                                           std::plus<>>,
                 decltype(A + A)>);
  static_assert(
    std::same_as<
      double, math::scalarize_via_cast_t<math::Array<D, SquareDims<>, true>>>);
  static_assert(
    std::same_as<double, math::scalarize_via_cast_t<decltype(A + A)>>);
  static_assert(
    std::same_as<math::ElementwiseBinaryOp<math::Array<D, SquareDims<>, true>,
                                           double, std::multiplies<>>,
                 decltype(A * 2.3)>);
  static_assert(
    std::same_as<double, math::scalarize_via_cast_t<decltype(A * 2.3)>>);
  static_assert(
    std::same_as<double, math::scalarize_via_cast_t<decltype(2.3 * A)>>);
  static_assert(
    math::ScalarizeViaCastTo<math::scalarize_via_cast_t<decltype(view(A))>,
                             decltype(2.3 * A)>());
#endif
  for (auto &&a : A) a = URand<D>{}(rng0);
  bench.run("BM_expm_dual7x2_size=" + std::to_string(size),
            [&] { expbench(A); });
}

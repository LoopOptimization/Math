
import Array;
import ExprTemplateUtils;
import ManagedArray;
import Nanobench;
import RandDual;
import StaticArray;
import std;
import Tuple;

using math::Dual, math::SquareDims, math::SquareMatrix, math::MutArray,
  math::Array, math::URand, containers::tie, containers::Tuple;

namespace {
#ifndef NDEBUG
template <typename A, typename... As, typename B, typename... Bs>
constexpr void tuplecheck(Tuple<A, As...> &, const Tuple<B, Bs...> &) {
  using C = math::scalarize_via_cast_t<
    std::remove_cvref_t<decltype(std::declval<A>().view())>>;
  static_assert(
    !std::same_as<C, void> &&
    detail::ScalarizeViaCastTo<C, As..., decltype(std::declval<B>().view()),
                               Bs...>());
}
#endif

template <typename T, typename S>
[[gnu::noinline]] void eltmul(MutArray<T, S> A, double t) {
  A *= t;
}

template <typename T, typename S>
[[gnu::noinline]] void eltadd(MutArray<T, S> C, Array<T, S> A, Array<T, S> B) {
  static_assert(
    std::same_as<double, math::scalarize_via_cast_t<MutArray<T, S>>>);
  static_assert(
    std::same_as<double, math::scalarize_via_cast_t<decltype(A + B)>>);
  C << A + B;
}
template <typename T, typename S>
[[gnu::noinline]] void eltaddsub(MutArray<T, S> C, MutArray<T, S> D,
                                 Array<T, S> A, Array<T, S> B) {
#ifndef NDEBUG
  {
    auto lval{tie(C, D)};
    tuplecheck(lval, Tuple(A + B, A - B));
  }
#endif
  tie(C, D) << Tuple(A + B, A - B);
}

template <std::ptrdiff_t N>
void BM_dualNdoublemul(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<double, N>;
  SquareMatrix<D> A{SquareDims{math::row(size)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  double t = 1.00000000001;
  bench.run([&] { eltmul(A, t); });
}

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualMxNdoublemul(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, M>, N>;
  SquareMatrix<D> A{SquareDims{math::row(size)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  double t = 1.00000000001;
  bench.run([&] { eltmul(A, t); });
}

template <std::ptrdiff_t N>
void BM_dualNadd(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using T = Dual<double, N>;
  SquareDims<> dim{math::row(size)};
  SquareMatrix<T> A{dim}, B{dim}, C{dim};
  for (auto &&a : A) a = URand<T>{}(rng0);
  for (auto &&b : B) b = URand<T>{}(rng0);
  bench.run([&] { eltadd(C, A, B); });
}
template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualMxNadd(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using T = Dual<Dual<double, M>, N>;
  SquareDims<> dim{math::row(size)};
  SquareMatrix<T> A{dim}, B{dim}, C{dim};
  for (auto &&a : A) a = URand<T>{}(rng0);
  for (auto &&b : B) b = URand<T>{}(rng0);
  bench.run([&] { eltadd(C, A, B); });
}
template <std::ptrdiff_t N>
void BM_dualNaddsub(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using T = Dual<double, N>;
  SquareDims<> dim{math::row(size)};
  SquareMatrix<T> A{dim}, B{dim}, C{dim}, D{dim};
  for (auto &&a : A) a = URand<T>{}(rng0);
  for (auto &&b : B) b = URand<T>{}(rng0);
  bench.run([&] { eltaddsub(C, D, A, B); });
}
template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualMxNaddsub(Bench &bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using T = Dual<Dual<double, M>, N>;
  SquareDims<> dim{math::row(size)};
  SquareMatrix<T> A{dim}, B{dim}, C{dim}, D{dim};
  for (auto &&a : A) a = URand<T>{}(rng0);
  for (auto &&b : B) b = URand<T>{}(rng0);
  bench.run([&] { eltaddsub(C, D, A, B); });
}
} // namespace

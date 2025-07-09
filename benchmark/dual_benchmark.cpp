
import Array;
import BaseUtils;
import ManagedArray;
import Nanobench;
import RandDual;
import SIMD;
import StaticArray;
import std;
import Tuple;

using math::Dual, math::SquareMatrix, math::URand;
namespace {

[[gnu::noinline]] void prod(auto &c, const auto &a, const auto &b) {
  c = a * b;
}

template <typename T, std::ptrdiff_t N, bool SIMDArray = false>
struct ManualDual {
  T value;
  math::SVector<T, N> partials;
  auto grad() -> math::SVector<T, N> & { return partials; }
};
template <std::floating_point T, std::ptrdiff_t N, bool SIMDArray>
struct ManualDual<ManualDual<T, N, SIMDArray>, 2, false> {
  using V = ManualDual<T, N, SIMDArray>;
  V value;
  containers::Tuple<V, V> partials{V{}, V{}};
  struct Gradient {
    containers::Tuple<V, V> &partials;
    auto operator[](std::ptrdiff_t i) -> V & {
      utils::invariant(i == 0 || i == 1);
      if (i == 0) return partials._0;
      return partials._1;
    }
  };
  constexpr auto grad() -> Gradient { return {partials}; }
  [[nodiscard]] constexpr auto grad() const -> std::array<V, 2> {
    return {partials._0, partials._1};
  }
};

template <std::floating_point T, std::ptrdiff_t N>
struct ManualDual<T, N, false> {
  using P = simd::Vec<std::ptrdiff_t(std::bit_ceil(std::size_t(N))), T>;
  T value;
  P partials;
  auto grad() -> P & { return partials; }
};
template <std::floating_point T, std::ptrdiff_t N>
struct ManualDual<T, N, true> {
  using P = math::StaticArray<T, 1, N, false>;
  T value;
  P partials;
  auto grad() -> P & { return partials; }
};
template <typename T, std::ptrdiff_t N, bool B>
[[gnu::always_inline]] constexpr auto operator*(ManualDual<T, N, B> a,
                                                ManualDual<T, N, B> b)
  -> ManualDual<T, N, B> {
  if constexpr ((!B) && (!std::floating_point<T>) && (N == 2))
    return {a.value * b.value,
            {a.value * b.grad()[0] + b.value + a.grad()[0],
             a.value * b.grad()[1] + b.value * a.grad()[1]}};
  else return {a.value * b.value, a.value * b.partials + b.value * a.partials};
}
template <typename T, std::ptrdiff_t N, bool B>
[[gnu::always_inline]] constexpr auto operator*(ManualDual<T, N, B> a, T b)
  -> ManualDual<T, N, B> {
  return {a.value * b, b * a.partials};
}
template <typename T, std::ptrdiff_t N, bool B>
[[gnu::always_inline]] constexpr auto operator*(T a, ManualDual<T, N, B> b)
  -> ManualDual<T, N, B> {
  return {b.value * a, a * b.partials};
}
template <typename T, std::ptrdiff_t N, bool B>
[[gnu::always_inline]] constexpr auto operator+(ManualDual<T, N, B> a,
                                                ManualDual<T, N, B> b)
  -> ManualDual<T, N, B> {
  return {a.value + b.value, a.partials + b.partials};
}
template <typename T, std::ptrdiff_t N, bool B>
[[gnu::always_inline]] constexpr auto operator+(ManualDual<T, N, B> a, T b)
  -> ManualDual<T, N, B> {
  return {a.value + b, a.partials};
}
template <typename T, std::ptrdiff_t N, bool B>
[[gnu::always_inline]] constexpr auto operator+(T a, ManualDual<T, N, B> b)
  -> ManualDual<T, N, B> {
  return {b.value + a, b.partials};
}

// template <typename T, std::ptrdiff_t M, std::ptrdiff_t N>
// [[gnu::noinline]] void prod_manual(ManualDual<ManualDual<T, M>, N> &c,
//                                    const ManualDual<ManualDual<T, M>, N> &a,
//                                    const ManualDual<ManualDual<T, M>, N> &b)
//                                    {
//   // return {val * other.val, val * other.partials + other.val * partials};
//   c.value= a.value* b.value;
//   c.partials = a.value* b.partials+ b.value* a.partials;
// }

template <std::ptrdiff_t M, std::ptrdiff_t N, bool SIMDArray, bool Outer>
auto setup_manual() {
  using D = ManualDual<ManualDual<double, M, SIMDArray>, N, Outer>;
  std::mt19937_64 rng0;
  D a{}, b{}, c{};
  a.value.value = URand<double>{}(rng0);
  b.value.value = URand<double>{}(rng0);
  for (std::ptrdiff_t j = 0; j < M; ++j) {
    a.value.partials[j] = URand<double>{}(rng0);
    b.value.partials[j] = URand<double>{}(rng0);
  }
  for (std::ptrdiff_t i = 0; i < N; ++i) {
    a.grad()[i].value = URand<double>{}(rng0);
    b.grad()[i].value = URand<double>{}(rng0);
    for (std::ptrdiff_t j = 0; j < M; ++j) {
      a.grad()[i].partials[j] = URand<double>{}(rng0);
      b.grad()[i].partials[j] = URand<double>{}(rng0);
    }
  }
  return std::array<D, 3>{a, b, c};
}
} // namespace
template <std::ptrdiff_t M, std::ptrdiff_t N> void BM_dualprod(Bench &bench) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, M>, N>;
  D a = URand<D>{}(rng0), b = URand<D>{}(rng0), c;
  bench.run("BM_dualprod<" + std::to_string(M) + "," + std::to_string(N) + ">",
            [&] {
              prod(c, a, b);
              doNotOptimizeAway(c);
            });
}

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_manual(Bench &bench) {
  auto [a, b, c] = setup_manual<M, N, false, true>();
  bench.run("BM_dualprod_manual<" + std::to_string(M) + "," +
              std::to_string(N) + ">",
            [&] {
              prod(c, a, b);
              doNotOptimizeAway(c);
            });
}
template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_simdarray(Bench &bench) {
  auto [a, b, c] = setup_manual<M, N, true, true>();
  bench.run("BM_dualprod_simdarray<" + std::to_string(M) + "," +
              std::to_string(N) + ">",
            [&] {
              prod(c, a, b);
              doNotOptimizeAway(c);
            });
}
template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_manual_tuple(Bench &bench) {
  auto [a, b, c] = setup_manual<M, N, false, false>();
  bench.run("BM_dualprod_manual_tuple<" + std::to_string(M) + "," +
              std::to_string(N) + ">",
            [&] {
              prod(c, a, b);
              doNotOptimizeAway(c);
            });
}
template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_simdarray_tuple(Bench &bench) {
  auto [a, b, c] = setup_manual<M, N, true, false>();
  bench.run("BM_dualprod_simdarray_tuple<" + std::to_string(M) + "," +
              std::to_string(N) + ">",
            [&] {
              prod(c, a, b);
              doNotOptimizeAway(c);
            });
}
template void BM_dualprod<6, 2>(Bench &bench);
template void BM_dualprod<7, 2>(Bench &bench);
template void BM_dualprod<8, 2>(Bench &bench);

template void BM_dualprod_manual<6, 2>(Bench &bench);
template void BM_dualprod_manual<7, 2>(Bench &bench);
template void BM_dualprod_manual<8, 2>(Bench &bench);

template void BM_dualprod_simdarray<6, 2>(Bench &bench);
template void BM_dualprod_simdarray<7, 2>(Bench &bench);
template void BM_dualprod_simdarray<8, 2>(Bench &bench);

template void BM_dualprod_manual_tuple<6, 2>(Bench &bench);
template void BM_dualprod_manual_tuple<7, 2>(Bench &bench);
template void BM_dualprod_manual_tuple<8, 2>(Bench &bench);

template void BM_dualprod_simdarray_tuple<6, 2>(Bench &bench);
template void BM_dualprod_simdarray_tuple<7, 2>(Bench &bench);
template void BM_dualprod_simdarray_tuple<8, 2>(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualdivsum(Bench &bench, std::ptrdiff_t len) {
  std::mt19937_64 rng0;
  using D =
    std::conditional_t<(N > 0), Dual<Dual<double, M>, N>, Dual<double, M>>;
  math::Vector<std::array<D, 4>> x{math::length(len)};
  for (std::ptrdiff_t i = 0; i < len; ++i) {
    x[i] = {URand<D>{}(rng0), URand<D>{}(rng0), URand<D>{}(rng0),
            URand<D>{}(rng0)};
  }
  bench.run("BM_dualdivsum<" + std::to_string(M) + "," + std::to_string(N) +
              ">_len=" + std::to_string(len),
            [&] {
              D s{};
              for (auto a : x) s += (a[0] + a[1]) / (a[2] + a[3]);
              doNotOptimizeAway(s);
            });
}

template void BM_dualdivsum<7, 0>(Bench &bench, std::ptrdiff_t len);
template void BM_dualdivsum<8, 0>(Bench &bench, std::ptrdiff_t len);
template void BM_dualdivsum<7, 2>(Bench &bench, std::ptrdiff_t len);
template void BM_dualdivsum<8, 2>(Bench &bench, std::ptrdiff_t len);
template void BM_dualdivsum<7, 4>(Bench &bench, std::ptrdiff_t len);
template void BM_dualdivsum<8, 4>(Bench &bench, std::ptrdiff_t len);

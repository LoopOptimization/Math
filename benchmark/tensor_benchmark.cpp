
#include "LoopMacros.hxx"
import Array;
import BaseUtils;
import ManagedArray;
import Nanobench;
import RandDual;
import SIMD;
import StaticArray;
import std;
import UniformScaling;

using math::Dual, math::SquareMatrix, math::SquareDims, math::I, math::URand;

[[gnu::noinline]] static void A12pI120(auto &B, const auto &A) {
  B << 12.0 * A + 120.0 * I;
}
[[gnu::noinline]] static void BApI60(auto &C, const auto &A, const auto &B) {
  C << B * (A + 60.0 * I);
}

void BM_dual8x2dApI(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 8>, 2>;
  SquareMatrix<D> A{SquareDims{math::row(size)}};
  SquareMatrix<D> B{SquareDims{math::row(size)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  bench.run([&] {
    A12pI120(B, A);
    doNotOptimizeAway(B);
  });
}

void BM_dual8x2BmApI(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 8>, 2>;
  SquareMatrix<D> A{SquareDims{math::row(size)}};
  SquareMatrix<D> B{SquareDims{math::row(size)}};
  SquareMatrix<D> C{SquareDims{math::row(size)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  for (auto &&b : B) b = URand<D>{}(rng0);
  bench.run([&] {
    BApI60(C, A, B);
    doNotOptimizeAway(C);
  });
}

template <std::size_t M, std::size_t N>
void BtimesAplusdI(
  math::MutSquarePtrMatrix<std::array<std::array<double, M>, N>> C,
  math::MutSquarePtrMatrix<std::array<std::array<double, M>, N>> A,
  math::MutSquarePtrMatrix<std::array<std::array<double, M>, N>> B,
  double doffset) {
  using T = std::array<std::array<double, M>, N>;
  utils::invariant(C.numRow() == A.numRow());
  utils::invariant(C.numRow() == B.numRow());
  std::ptrdiff_t D = std::ptrdiff_t(C.numRow());
  POLYMATHNOVECTORIZE
  for (std::ptrdiff_t r = 0; r < D; ++r) {
    POLYMATHNOVECTORIZE
    for (std::ptrdiff_t c = 0; c < D; ++c) {
      T x{};
      POLYMATHNOVECTORIZE
      for (std::ptrdiff_t k = 0; k < D; ++k) {
        // x += B[r, k] * A[k, c];
        T &Brk = B[r, k];
        T &Akc = A[k, c];
        // x[0] += Brk[0] * Akc[0];
        x[0][0] += Brk[0][0] * (Akc[0][0] + doffset * (r == c));
        POLYMATHVECTORIZE
        for (std::size_t i = 1; i < M; ++i)
          x[0][i] += Brk[0][0] * Akc[0][i] + Brk[0][i] * Akc[0][0];
        POLYMATHNOVECTORIZE
        for (std::size_t o = 1; o < N; ++o) {
          // x[o] += Brk[0]*Akc[o];
          x[o][0] += Brk[0][0] * Akc[o][0];
          POLYMATHVECTORIZE
          for (std::size_t i = 1; i < M; ++i)
            x[o][i] += Brk[0][0] * Akc[o][i] + Brk[0][i] * Akc[o][0];
          // x[o] += Brk[o]*Akc[0];
          x[o][0] += Brk[o][0] * Akc[0][0];
          POLYMATHVECTORIZE
          for (std::size_t i = 1; i < M; ++i)
            x[o][i] += Brk[o][0] * Akc[0][i] + Brk[o][i] * Akc[0][0];
        }
      }
      C[r, c] = x;
    }
  }
}

void BM_dual8x2BmApI_manual(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = std::array<std::array<double, 9>, 3>;
  std::ptrdiff_t dim = size;
  SquareMatrix<D> A{SquareDims{math::row(dim)}};
  SquareMatrix<D> B{SquareDims{math::row(dim)}};
  SquareMatrix<D> C{SquareDims{math::row(dim)}};
  for (std::ptrdiff_t i = 0, L = dim * dim; i < L; ++i) {
    for (std::ptrdiff_t j = 0; j < 3; ++j) {
      for (std::ptrdiff_t k = 0; k < 9; ++k) {
        A[i][j][k] = URand<double>{}(rng0);
        B[i][j][k] = URand<double>{}(rng0);
      }
    }
  }
  bench.run([&] { BtimesAplusdI(C, A, B, 60.0); });
}

void BM_dual7x2dApI(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 7>, 2>;
  static_assert(utils::Compressible<Dual<double, 7>>);
  static_assert(utils::Compressible<D>);
  static_assert(sizeof(utils::compressed_t<D>) == (24 * sizeof(double)));
  static_assert(sizeof(D) == (24 * sizeof(double)));
  std::ptrdiff_t dim = size;
  SquareMatrix<D> A{SquareDims{math::row(dim)}};
  SquareMatrix<D> B{SquareDims{math::row(dim)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  bench.run([&] {
    A12pI120(B, A);
    doNotOptimizeAway(B);
  });
}

void BM_dual7x2BmApI(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 7>, 2>;
  std::ptrdiff_t dim = size;
  SquareMatrix<D> A{SquareDims{math::row(dim)}};
  SquareMatrix<D> B{SquareDims{math::row(dim)}};
  SquareMatrix<D> C{SquareDims{math::row(dim)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  for (auto &&b : B) b = URand<D>{}(rng0);
  bench.run([&] {
    BApI60(C, A, B);
    doNotOptimizeAway(C);
  });
}

void BM_dual7x2BmApI_manual(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = std::array<std::array<double, 8>, 3>;
  std::ptrdiff_t dim = size;
  SquareMatrix<D> A{SquareDims{math::row(dim)}};
  SquareMatrix<D> B{SquareDims{math::row(dim)}};
  SquareMatrix<D> C{SquareDims{math::row(dim)}};
  for (std::ptrdiff_t i = 0, L = dim * dim; i < L; ++i) {
    for (std::ptrdiff_t j = 0; j < 3; ++j) {
      for (std::ptrdiff_t k = 0; k < 8; ++k) {
        A[i][j][k] = URand<double>{}(rng0);
        B[i][j][k] = URand<double>{}(rng0);
      }
    }
  }
  bench.run([&] { BtimesAplusdI(C, A, B, 60.0); });
}

void BM_dual6x2dApI(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 6>, 2>;
  static_assert(utils::Compressible<Dual<double, 6>>);
  static_assert(utils::Compressible<D>);
  static_assert(sizeof(utils::compressed_t<D>) == (21 * sizeof(double)));
  static_assert(simd::SIMDSupported<double> ==
                (sizeof(math::SVector<double, 7>) == (8 * sizeof(double))));
  static_assert(simd::SIMDSupported<double> ==
                (sizeof(Dual<double, 6>) == (8 * sizeof(double))));
  static_assert(simd::SIMDSupported<double> ==
                (sizeof(Dual<double, 6>) == (8 * sizeof(double))));
  static_assert(simd::SIMDSupported<double> ==
                (sizeof(D) == (24 * sizeof(double))));
  // static_assert(sizeof(D) == sizeof(Dual<Dual<double, 8>, 2>));
  std::ptrdiff_t dim = size;
  SquareMatrix<D> A{SquareDims{math::row(dim)}};
  SquareMatrix<D> B{SquareDims{math::row(dim)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  bench.run([&] {
    A12pI120(B, A);
    doNotOptimizeAway(B);
  });
}

void BM_dual6x2BmApI(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  using D = Dual<Dual<double, 6>, 2>;
  std::ptrdiff_t dim = size;
  SquareMatrix<D> A{SquareDims{math::row(dim)}};
  SquareMatrix<D> B{SquareDims{math::row(dim)}};
  SquareMatrix<D> C{SquareDims{math::row(dim)}};
  for (auto &&a : A) a = URand<D>{}(rng0);
  for (auto &&b : B) b = URand<D>{}(rng0);
  bench.run([&] {
    BApI60(C, A, B);
    doNotOptimizeAway(C);
  });
}

void BM_dual6x2BmApI_manual(Bench& bench, std::ptrdiff_t size) {
  std::mt19937_64 rng0;
  constexpr std::size_t Dcount = 6;
  constexpr std::size_t N = Dcount + 1;
  using D = std::array<std::array<double, N>, 3>;
  std::ptrdiff_t dim = size;
  SquareMatrix<D> A{SquareDims{math::row(dim)}};
  SquareMatrix<D> B{SquareDims{math::row(dim)}};
  SquareMatrix<D> C{SquareDims{math::row(dim)}};
  for (std::ptrdiff_t i = 0, L = dim * dim; i < L; ++i) {
    for (std::ptrdiff_t j = 0; j < 3; ++j) {
      for (std::size_t k = 0; k < N; ++k) {
        A[i][j][k] = URand<double>{}(rng0);
        B[i][j][k] = URand<double>{}(rng0);
      }
    }
  }
  bench.run([&] { BtimesAplusdI(C, A, B, 60.0); });
}

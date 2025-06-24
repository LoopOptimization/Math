#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#include "Macros.hxx"
#ifndef USE_MODULE
#include "Alloc/Arena.cxx"
#include "Containers/Tuple.cxx"
#include "Math/Array.cxx"
#include "Math/ArrayConcepts.cxx"
#include "Math/AxisTypes.cxx"
#include "Math/Comparisons.cxx"
#include "Math/Constructors.cxx"
#include "Math/EmptyArrays.cxx"
#include "Math/GenericConstructors.cxx"
#include "Math/GreatestCommonDivisor.cxx"
#include "Math/Indexing.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/MatrixDimensions.cxx"
#include "Math/Ranges.cxx"
#include "Math/Rational.cxx"
#include "Math/VectorGreatestCommonDivisor.cxx"
#include "SIMD/Intrin.cxx"
#include "SIMD/Masks.cxx"
#include "SIMD/UnrollIndex.cxx"
#include "SIMD/Vec.cxx"
#include "Utilities/Invariant.cxx"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#else
export module NormalForm;
import Arena;
import Array;
import ArrayConcepts;
import ArrayConstructors;
import AxisTypes;
import Comparisons;
import EmptyMatrix;
import GenericArrayConstructors;
import Invariant;
import ManagedArray;
import Range;
import Rational;
import SIMD;
import std;
import Tuple;
import VGCD;
#endif

namespace math::NormalForm {
using alloc::Arena, containers::Tuple, containers::tie;

using namespace math;
namespace detail {

TRIVIAL constexpr auto gcdxScale(std::int64_t a, std::int64_t b)
  -> std::array<std::int64_t, 4> {
  if (constexpr_abs(a) == 1) return {a, 0, a, b};
  auto [g, p, q, adg, bdg] = dgcdx(a, b);
  return {p, q, adg, bdg};
}
// zero out below diagonal
TRIVIAL constexpr void zeroSupDiagonal(MutPtrMatrix<std::int64_t> A,
                                       MutSquarePtrMatrix<std::int64_t> K,
                                       std::ptrdiff_t i, Row<> M, Col<> N) {
  std::ptrdiff_t minMN = std::min(std::ptrdiff_t(M), std::ptrdiff_t(N));
  for (std::ptrdiff_t j = i + 1; j < M; ++j) {
    std::int64_t Aii = A[i, i];
    if (std::int64_t Aji = A[j, i]) {
      const auto [p, q, Aiir, Aijr] = gcdxScale(Aii, Aji);

      {
        MutPtrVector<std::int64_t> Ai{A[i, _(0, minMN)]}, Aj{A[j, _(0, minMN)]},
          Ki{K[i, _(0, minMN)]}, Kj{K[j, _(0, minMN)]};
        tie(Ai, Aj, Ki, Kj) << Tuple(p * Ai + q * Aj, Aiir * Aj - Aijr * Ai,
                                     p * Ki + q * Kj, Aiir * Kj - Aijr * Ki);
      }
      if (std::ptrdiff_t(M) > std::ptrdiff_t(N)) {
        MutPtrVector<std::int64_t> Ki{K[i, _(N, M)]}, Kj{K[j, _(N, M)]};
        tie(Ki, Kj) << Tuple(p * Ki + q * Kj, Aiir * Kj - Aijr * Ki);
      } else if (std::ptrdiff_t(N) > std::ptrdiff_t(M)) {
        MutPtrVector<std::int64_t> Ai{A[i, _(M, N)]}, Aj{A[j, _(M, N)]};
        tie(Ai, Aj) << Tuple(p * Ai + q * Aj, Aiir * Aj - Aijr * Ai);
      }
      // for (std::ptrdiff_t k = 0; k < minMN; ++k) {
      //   std::int64_t Aki = A[i, k];
      //   std::int64_t Akj = A[j, k];
      //   std::int64_t Kki = K[i, k];
      //   std::int64_t Kkj = K[j, k];
      //   // when k == i, then
      //   // p * Aii + q * Akj == r, so we set A(i,i) = r
      //   A[i, k] = p * Aki + q * Akj;
      //   // Aii/r * Akj - Aij/r * Aki = 0
      //   A[j, k] = Aiir * Akj - Aijr * Aki;
      //   // Mirror for K
      //   K[i, k] = p * Kki + q * Kkj;
      //   K[j, k] = Aiir * Kkj - Aijr * Kki;
      // }
      // for (auto k = std::ptrdiff_t(N); k < M; ++k) {
      //   std::int64_t Kki = K[i, k];
      //   std::int64_t Kkj = K[j, k];
      //   K[i, k] = p * Kki + q * Kkj;
      //   K[j, k] = Aiir * Kkj - Aijr * Kki;
      // }
      // for (auto k = std::ptrdiff_t(M); k < N; ++k) {
      //   std::int64_t Aki = A[i, k];
      //   std::int64_t Akj = A[j, k];
      //   A[i, k] = p * Aki + q * Akj;
      //   A[j, k] = Aiir * Akj - Aijr * Aki;
      // }
    }
  }
}
// This method is only called by orthogonalize, hence we can assume
// (Akk == 1) || (Akk == -1)
TRIVIAL constexpr void zeroSubDiagonal(MutPtrMatrix<std::int64_t> A,
                                       MutSquarePtrMatrix<std::int64_t> K,
                                       std::ptrdiff_t k, Row<> M, Col<> N) {
  std::int64_t Akk = A[k, k];
  if (Akk == -1) {
    for (std::ptrdiff_t m = 0; m < N; ++m) A[k, m] *= -1;
    for (std::ptrdiff_t m = 0; m < M; ++m) K[k, m] *= -1;
  } else {
    invariant(Akk == 1);
  }
  std::ptrdiff_t Mi = std::ptrdiff_t(M), Ni = std::ptrdiff_t(N),
                 minMN = std::min(Mi, Ni);
  for (std::ptrdiff_t z = 0; z < k; ++z) {
    // eliminate `A(k,z)`
    std::int64_t Akz = A[z, k];
    if (!Akz) continue;
    // A(k, k) == 1, so A(k,z) -= Akz * 1;
    // A(z,_) -= Akz * A(k,_);
    // K(z,_) -= Akz * K(k,_);
    Tuple(A[z, _(0, minMN)], K[z, _(0, minMN)]) -=
      Tuple(Akz * A[k, _(0, minMN)], Akz * K[k, _(0, minMN)]);
    if (Mi > Ni) K[z, _(Ni, Mi)] -= Akz * K[k, _(Ni, Mi)];
    else if (Ni > Mi) A[z, _(Mi, Ni)] -= Akz * A[k, _(Mi, Ni)];
    // for (std::ptrdiff_t i = 0; i < minMN; ++i) {
    //   A[z, i] -= Akz * A[k, i];
    //   K[z, i] -= Akz * K[k, i];
    // }
    // for (auto i = std::ptrdiff_t(N); i < M; ++i) K[z, i] -= Akz * K[k, i];
    // for (auto i = std::ptrdiff_t(M); i < N; ++i) A[z, i] -= Akz * A[k, i];
  }
}

TRIVIAL constexpr auto
pivotRowsPair(std::array<MutPtrMatrix<std::int64_t>, 2> AK, Col<> i, Row<> M,
              Row<> piv) -> bool {
  Row j = piv;
  while (AK[0][piv, i] == 0)
    if (++piv == M) return true;
  if (j != piv) {
    math::swap(AK[0], j, piv);
    math::swap(AK[1], j, piv);
  }
  return false;
}
TRIVIAL constexpr auto
pivotColsPair(std::array<MutPtrMatrix<std::int64_t>, 2> AK, Row<> i, Col<> N,
              Col<> piv) -> bool {
  Col j = piv;
  while (AK[0][i, piv] == 0)
    if (++piv == N) return true;
  if (j != piv) {
    math::swap(AK[0], j, piv);
    math::swap(AK[1], j, piv);
  }
  return false;
}
} // namespace detail
} // namespace math::NormalForm

#ifdef USE_MODULE
export namespace math::NormalForm {
#else
namespace math::NormalForm {
#endif
TRIVIAL constexpr auto pivotRows(MutPtrMatrix<std::int64_t> A,
                                 MutSquarePtrMatrix<std::int64_t> K,
                                 std::ptrdiff_t i, Row<> M) -> bool {
  MutPtrMatrix<std::int64_t> B = K;
  return detail::pivotRowsPair({A, B}, col(i), M, row(i));
}
TRIVIAL constexpr auto pivotRows(MutPtrMatrix<std::int64_t> A, Col<> i, Row<> M,
                                 Row<> piv) -> bool {
  Row j = piv;
  while (A[piv, i] == 0)
    if (++piv == std::ptrdiff_t(M)) return true;
  if (j != piv) swap(A, j, piv);
  return false;
}
TRIVIAL constexpr auto pivotRows(MutPtrMatrix<std::int64_t> A, std::ptrdiff_t i,
                                 Row<> N) -> bool {
  return pivotRows(A, col(i), N, row(i));
}
/// numNonZeroRows(PtrMatrix<std::int64_t> A) -> Row
/// Assumes some number of the trailing rows have been
/// zeroed out.  Returns the number of rows that are remaining.
TRIVIAL constexpr auto numNonZeroRows(PtrMatrix<std::int64_t> A) -> Row<> {
  Row newM = A.numRow();
  while (newM && allZero(A[std::ptrdiff_t(newM) - 1, _])) --newM;
  return newM;
}
} // namespace math::NormalForm
namespace math::NormalForm::detail {

TRIVIAL constexpr void dropCol(MutPtrMatrix<std::int64_t> A, std::ptrdiff_t i,
                               Row<> M, Col<> N) {
  // if any rows are left, we shift them up to replace it
  if (N <= i) return;
  A[_(0, M), _(i, N)] << A[_(0, M), _(i, N) + 1];
  // for (std::ptrdiff_t m = 0; m < M; ++m) A[m, _(i, N)] << A[m, _(i, N) + 1];
  // for (std::ptrdiff_t n = i; n < N; ++n) A[m, n] = A[m, n + 1];
}

TRIVIAL constexpr void zeroSupDiagonal(MutPtrMatrix<std::int64_t> A, Col<> c,
                                       Row<> r) {
  auto [M, N] = shape(A);
  for (std::ptrdiff_t j = std::ptrdiff_t(r) + 1; j < M; ++j) {
    std::int64_t Aii = A[r, c];
    if (std::int64_t Aij = A[j, c]) {
      const auto [p, q, Aiir, Aijr] = gcdxScale(Aii, Aij);
      MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
      tie(Ar, Aj) << Tuple(p * Ar + q * Aj, Aiir * Aj - Aijr * Ar);
    }
  }
}
TRIVIAL constexpr void
zeroSupDiagonal(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c,
                Row<> r) {
  auto [A, B] = AB;
  auto [M, N] = shape(A);
  invariant(M, std::ptrdiff_t(B.numRow()));
  for (std::ptrdiff_t j = std::ptrdiff_t(r) + 1; j < M; ++j) {
    std::int64_t Aii = A[r, c], Aij = A[j, c];
    if (!Aij) continue;
    const auto [p, q, Aiir, Aijr] = gcdxScale(Aii, Aij);
    MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    tie(Ar, Aj) << Tuple(p * Ar + q * Aj, Aiir * Aj - Aijr * Ar);
    MutPtrVector<std::int64_t> Br = B[r, _], Bj = B[j, _];
    tie(Br, Bj) << Tuple(p * Br + q * Bj, Aiir * Bj - Aijr * Br);
  }
}
TRIVIAL constexpr void reduceSubDiagonal(MutPtrMatrix<std::int64_t> A, Col<> c,
                                         Row<> r) {
  std::int64_t Akk = A[r, c];
  if (Akk < 0) {
    Akk = -Akk;
    A[r, _] *= -1;
  }
  for (std::ptrdiff_t z = 0; z < r; ++z) {
    // try to eliminate `A(k,z)`
    // if Akk == 1, then this zeros out Akz
    if (std::int64_t Azr = A[z, c]) {
      // we want positive but smaller subdiagonals
      // e.g., `Akz = 5, Akk = 2`, then in the loop below when `i=k`, we
      // set A(k,z) = A(k,z) - (A(k,z)/Akk) * Akk
      //        =   5 - 2*2 = 1
      // or if `Akz = -5, Akk = 2`, then in the loop below we get
      // A(k,z) = A(k,z) - ((A(k,z)/Akk) - ((A(k,z) % Akk) != 0) * Akk
      //        =  -5 - (-2 - 1)*2 = = 6 - 5 = 1
      // if `Akk = 1`, then
      // A(k,z) = A(k,z) - (A(k,z)/Akk) * Akk
      //        = A(k,z) - A(k,z) = 0
      // or if `Akz = -7, Akk = 39`, then in the loop below we get
      // A(k,z) = A(k,z) - ((A(k,z)/Akk) - ((A(k,z) % Akk) != 0) * Akk
      //        =  -7 - ((-7/39) - 1)*39 = = 6 - 5 = 1
      std::int64_t oAzr = Azr;
      Azr /= Akk;
      if (oAzr < 0) Azr -= (oAzr != (Azr * Akk));
      A[z, _] -= Azr * A[r, _];
    }
  }
}
TRIVIAL constexpr void reduceSubDiagonalStack(MutPtrMatrix<std::int64_t> A,
                                              MutPtrMatrix<std::int64_t> B,
                                              std::ptrdiff_t c,
                                              std::ptrdiff_t r) {
  std::int64_t Akk = A[r, c];
  if (Akk < 0) {
    Akk = -Akk;
    A[r, _] *= -1;
  }
  for (std::ptrdiff_t z = 0; z < r; ++z) {
    if (std::int64_t Akz = A[z, c]) {
      std::int64_t oAkz = Akz;
      Akz /= Akk;
      if (oAkz < 0) Akz -= (oAkz != (Akz * Akk));
      A[z, _] -= Akz * A[r, _];
    }
  }
  for (std::ptrdiff_t z = 0; z < B.numRow(); ++z) {
    if (std::int64_t Bzr = B[z, c]) {
      std::int64_t oBzr = Bzr;
      Bzr /= Akk;
      if (oBzr < 0) Bzr -= (oBzr != (Bzr * Akk));
      B[z, _] -= Bzr * A[r, _];
    }
  }
}
TRIVIAL constexpr void
reduceSubDiagonal(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c,
                  Row<> r) {
  auto [A, B] = AB;
  std::int64_t Akk = A[r, c];
  if (Akk < 0) {
    Akk = -Akk;
    A[r, _] *= -1;
    B[r, _] *= -1;
  }
  for (std::ptrdiff_t z = 0; z < r; ++z) {
    // try to eliminate `A(k,z)`
    std::int64_t Akz = A[z, c];
    if (!Akz) continue;
    // if Akk == 1, then this zeros out Akz
    if (Akk != 1) {
      // we want positive but smaller subdiagonals
      // e.g., `Akz = 5, Akk = 2`, then in the loop below when `i=k`,
      // we set A(k,z) = A(k,z) - (A(k,z)/Akk) * Akk
      //        =   5 - 2*2 = 1
      // or if `Akz = -5, Akk = 2`, then in the loop below we get
      // A(k,z) = A(k,z) - ((A(k,z)/Akk) - ((A(k,z) % Akk) != 0) * Akk
      //        =  -5 - (-2 - 1)*2 = = 6 - 5 = 1
      // if `Akk = 1`, then
      // A(k,z) = A(k,z) - (A(k,z)/Akk) * Akk
      //        = A(k,z) - A(k,z) = 0
      // or if `Akz = -7, Akk = 39`, then in the loop below we get
      // A(k,z) = A(k,z) - ((A(k,z)/Akk) - ((A(k,z) % Akk) != 0) * Akk
      //        =  -7 - ((-7/39) - 1)*39 = = 6 - 5 = 1
      std::int64_t oAkz = Akz;
      Akz /= Akk;
      if (oAkz < 0) Akz -= (oAkz != (Akz * Akk));
    }
    A[z, _] -= Akz * A[r, _];
    B[z, _] -= Akz * B[r, _];
  }
}

TRIVIAL constexpr void reduceColumn(MutPtrMatrix<std::int64_t> A, Col<> c,
                                    Row<> r) {
  zeroSupDiagonal(A, c, r);
  reduceSubDiagonal(A, c, r);
}

// NormalForm version assumes zero rows are sorted to end due to pivoting
TRIVIAL constexpr void removeZeroRows(MutDensePtrMatrix<std::int64_t> &A) {
  A.truncate(NormalForm::numNonZeroRows(A));
}

TRIVIAL constexpr void
reduceColumn(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c, Row<> r) {
  zeroSupDiagonal(AB, c, r);
  reduceSubDiagonal(AB, c, r);
}
/// multiplies `A` and `B` by matrix `X`, where `X` reduces `A` to a normal
/// form.
TRIVIAL constexpr void
zeroWithRowOperation(MutPtrMatrix<std::int64_t> A, Row<> i, Row<> j, Col<> k,
                     Range<std::ptrdiff_t, std::ptrdiff_t> skip) {
  if (std::int64_t Aik = A[i, k]) {
    std::int64_t Ajk = A[j, k];
    std::int64_t g = gcd(Aik, Ajk);
    Aik /= g;
    Ajk /= g;
    g = 0;
    for (std::ptrdiff_t l = 0; l < skip.b_; ++l) {
      std::int64_t Ail = Ajk * A[i, l] - Aik * A[j, l];
      A[i, l] = Ail;
      g = gcd(Ail, g);
    }
    for (std::ptrdiff_t l = skip.e_; l < A.numCol(); ++l) {
      std::int64_t Ail = Ajk * A[i, l] - Aik * A[j, l];
      A[i, l] = Ail;
      g = gcd(Ail, g);
    }
    if (g > 1) {
      for (std::ptrdiff_t l = 0; l < skip.b_; ++l)
        if (std::int64_t Ail = A[i, l]) A[i, l] = Ail / g;
      for (std::ptrdiff_t l = skip.e_; l < A.numCol(); ++l)
        if (std::int64_t Ail = A[i, l]) A[i, l] = Ail / g;
    }
  }
}

// use row `r` to zero the remaining rows of column `c`
TRIVIAL constexpr void
zeroColumnPair(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c, Row<> r) {
  auto [A, B] = AB;
  const Row M = A.numRow();
  invariant(M, B.numRow());
  for (std::ptrdiff_t j = 0; j < r; ++j) {
    std::int64_t Arc = A[r, c], Ajc = A[j, c];
    if (!Ajc) continue;
    std::int64_t g = gcd(Arc, Ajc), x = Arc / g, y = Ajc / g;
    // auto [x, y] = divgcd(Arc, Ajc);
    // MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    // Aj << x * Aj - y * Ar;
    // MutPtrVector<std::int64_t> Br = B[r, _], Bj = B[j, _];
    // Bj << x * Bj - y * Br;
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutPtrVector<std::int64_t> Ar = AB[i][r, _], Aj = AB[i][j, _];
      Aj << x * Aj - y * Ar;
    }
  }
  // greater rows in previous columns have been zeroed out
  // therefore it is safe to use them for row operations with this row
  for (auto j = std::ptrdiff_t(r); ++j < M;) {
    std::int64_t Arc = A[r, c], Ajc = A[j, c];
    if (!Ajc) continue;
    const auto [p, q, Arcr, Ajcr] = detail::gcdxScale(Arc, Ajc);
    // MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    // tie(Ar, Aj) << Tuple(q * Aj + p * Ar, Arcr * Aj - Ajcr * Ar);
    // MutPtrVector<std::int64_t> Br = B[r, _], Bj = B[j, _];
    // tie(Br, Bj) << Tuple(q * Bj + p * Br, Arcr * Bj - Ajcr * Br);
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutPtrVector<std::int64_t> Ar = AB[i][r, _], Aj = AB[i][j, _];
      tie(Ar, Aj) << Tuple(q * Aj + p * Ar, Arcr * Aj - Ajcr * Ar);
    }
  }
}
// use col `c` to zero the remaining cols of row `r`
TRIVIAL constexpr void
zeroColumnPair(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Row<> r, Col<> c) {
  auto [A, B] = AB;
  const Col N = A.numCol();
  invariant(N, B.numCol());
  for (std::ptrdiff_t j = 0; j < c; ++j) {
    std::int64_t Arc = A[r, c], Arj = A[r, j];
    if (!Arj) continue;
    // std::int64_t g = gcd(Arc, Arj), x = Arc / g, y = Arj / g;
    auto [x, y] = divgcd(Arc, Arj);
    // MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    // Aj << x * Aj - y * Ar;
    // MutPtrVector<std::int64_t> Br = B[r, _], Bj = B[j, _];
    // Bj << x * Bj - y * Br;
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutArray<std::int64_t, StridedRange<>> Ac = AB[i][_, c], Aj = AB[i][_, j];
      Aj << x * Aj - y * Ac;
    }
  }
  // greater rows in previous columns have been zeroed out
  // therefore it is safe to use them for row operations with this row
  for (auto j = std::ptrdiff_t(c); ++j < N;) {
    std::int64_t Arc = A[r, c], Arj = A[r, j];
    if (!Arj) continue;
    const auto [p, q, Arcr, Ajcr] = detail::gcdxScale(Arc, Arj);
    // MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    // tie(Ar, Aj) << Tuple(q * Aj + p * Ar, Arcr * Aj - Ajcr * Ar);
    // MutPtrVector<std::int64_t> Br = B[r, _], Bj = B[j, _];
    // tie(Br, Bj) << Tuple(q * Bj + p * Br, Arcr * Bj - Ajcr * Br);
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutArray<std::int64_t, StridedRange<>> Ac = AB[i][_, c], Aj = AB[i][_, j];
      tie(Ac, Aj) << Tuple(q * Aj + p * Ac, Arcr * Aj - Ajcr * Ac);
    }
  }
}
// use row `r` to zero the remaining rows of column `c`
TRIVIAL constexpr void zeroColumn(MutPtrMatrix<std::int64_t> A, Col<> c,
                                  Row<> r) {
  const Row M = A.numRow();
  for (std::ptrdiff_t j = 0; j < r; ++j) {
    std::int64_t Arc = A[r, c], Ajc = A[j, c];
    invariant(Arc != std::numeric_limits<std::int64_t>::min());
    invariant(Ajc != std::numeric_limits<std::int64_t>::min());
    if (!Ajc) continue;
    std::int64_t g = gcd(Arc, Ajc);
    A[j, _] << (Arc / g) * A[j, _] - (Ajc / g) * A[r, _];
  }
  // greater rows in previous columns have been zeroed out
  // therefore it is safe to use them for row operations with this row
  for (std::ptrdiff_t j = std::ptrdiff_t(r) + 1; j < M; ++j) {
    std::int64_t Arc = A[r, c], Ajc = A[j, c];
    if (!Ajc) continue;
    const auto [p, q, Arcr, Ajcr] = detail::gcdxScale(Arc, Ajc);
    MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    tie(Ar, Aj) << Tuple(q * Aj + p * Ar, Arcr * Aj - Ajcr * Ar);
  }
}

TRIVIAL constexpr auto pivotRowsBareiss(MutPtrMatrix<std::int64_t> A,
                                        std::ptrdiff_t i, Row<> M, Row<> piv)
  -> std::optional<std::ptrdiff_t> {
  Row j = piv;
  while (A[piv, i] == 0)
    if (++piv == M) return {};
  if (j != piv) swap(A, j, piv);
  return std::ptrdiff_t(piv);
}

TRIVIAL constexpr auto orthogonalizeBang(MutDensePtrMatrix<std::int64_t> &A)
  -> containers::Pair<SquareMatrix<std::int64_t>, Vector<unsigned>> {
  // we try to orthogonalize with respect to as many rows of `A` as we can
  // prioritizing earlier rows.
  auto [M, N] = shape(A);
  SquareMatrix<std::int64_t> K{
    identity(math::DefaultAlloc<std::int64_t>{}, std::ptrdiff_t(M))};
  Vector<unsigned> included;
  included.reserve(std::min(std::ptrdiff_t(M), std::ptrdiff_t(N)));
  for (std::ptrdiff_t i = 0, j = 0;
       i < std::min(std::ptrdiff_t(M), std::ptrdiff_t(N)); ++j) {
    // zero ith row
    if (NormalForm::pivotRows(A, K, i, row(M))) {
      // cannot pivot, this is a linear combination of previous
      // therefore, we drop the row
      dropCol(A, i, row(M), col(--N));
    } else {
      zeroSupDiagonal(A, K, i, row(M), col(N));
      std::int64_t Aii = A[i, i];
      if (math::constexpr_abs(Aii) != 1) {
        // including this row renders the matrix not unimodular!
        // therefore, we drop the row.
        dropCol(A, i, row(M), col(--N));
      } else {
        // we zero the sub diagonal
        zeroSubDiagonal(A, K, i++, row(M), col(N));
        included.push_back(j);
      }
    }
  }
  return {std::move(K), std::move(included)};
}
} // namespace math::NormalForm::detail

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif
namespace NormalForm {
// update a reduced matrix for a new row
// doesn't reduce last row (assumes you're solving for it)
TRIVIAL constexpr auto updateForNewRow(MutPtrMatrix<std::int64_t> A)
  -> std::ptrdiff_t {
  // use existing rows to reduce
  std::ptrdiff_t M = std::ptrdiff_t(A.numRow()), N = std::ptrdiff_t(A.numCol()),
                 MM = M - 1, NN = N - 1, n = 0, i,
                 j = std::numeric_limits<std::ptrdiff_t>::max();
  for (std::ptrdiff_t m = 0; m < MM; ++m) {
#ifndef NDEBUG
    if (!allZero(A[m, _(0, n)])) __builtin_trap();
#endif
    while (A[m, n] == 0) {
      if ((j > NN) && (A[MM, n] != 0)) {
        i = m;
        j = n;
      }
      invariant((++n) < NN);
    }
    if (std::int64_t Aln = A[MM, n]) {
      // use this to reduce last row
      auto [x, y] = divgcd(Aln, A[m, n]);
      A[MM, _] << A[MM, _] * y - A[m, _] * x;
      invariant(A[MM, n] == 0);
    }
    ++n;
  }
  // we've reduced the new row, now to use it...
  // swap A[i,_(j,end)] with A[MM,_(j,end)]
  if (j <= NN) { // we could do with a lot less copying...
    for (std::ptrdiff_t l = i; l < MM; ++l)
      for (std::ptrdiff_t k = j; k < N; ++k) std::swap(A[l, k], A[MM, k]);
  } else {
    // maybe there is a non-zero value
    j = n;
    for (; j < NN; ++j)
      if (A[MM, j] != 0) break;
    if (j == NN) return MM;
    i = MM;
  }
  // zero out A(_(0,i),j) using A(i,j)
  for (std::ptrdiff_t k = 0; k < i; ++k) {
    if (std::int64_t Akj = A[k, j]) {
      auto [x, y] = divgcd(Akj, A[i, j]);
      A[k, _] << A[k, _] * y - A[i, _] * x;
    }
  }
  return M;
}
TRIVIAL constexpr void
simplifySystemsImpl(std::array<MutPtrMatrix<std::int64_t>, 2> AB) {
  auto [M, N] = shape(AB[0]);
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c)
    if (!detail::pivotRowsPair(AB, col(c), row(M), row(r)))
      detail::reduceColumn(AB, col(c), row(r++));
}

// treats A as stacked on top of B
TRIVIAL constexpr void reduceColumnStack(MutPtrMatrix<std::int64_t> A,
                                         MutPtrMatrix<std::int64_t> B,
                                         std::ptrdiff_t c, std::ptrdiff_t r) {
  detail::zeroSupDiagonal(B, col(c), row(r));
  detail::reduceSubDiagonalStack(B, A, c, r);
}
// pass by value, returns number of rows to truncate
TRIVIAL constexpr auto simplifySystemImpl(MutPtrMatrix<std::int64_t> A,
                                          std::ptrdiff_t colInit = 0) -> Row<> {
  auto [M, N] = shape(A);
  for (std::ptrdiff_t r = 0, c = colInit; c < N && r < M; ++c)
    if (!pivotRows(A, col(c), row(M), row(r)))
      detail::reduceColumn(A, col(c), row(r++));
  return numNonZeroRows(A);
}

TRIVIAL constexpr void simplifySystem(EmptyMatrix<std::int64_t>,
                                      std::ptrdiff_t = 0) {}
TRIVIAL constexpr void simplifySystem(MutPtrMatrix<std::int64_t> &E,
                                      std::ptrdiff_t colInit = 0) {
  E.truncate(simplifySystemImpl(E, colInit));
}
// TODO: `const IntMatrix &` can be copied to `MutPtrMatrix<std::int64_t>`
// this happens via `const IntMatrix &` -> `const MutPtrMatrix<std::int64_t> &`
// -> `MutPtrMatrix<std::int64_t>`. Perhaps we should define `MutPtrMatrix(const
// MutPtrMatrix &) = delete;`?
//
// NOLINTNEXTLINE(performance-unnecessary-value-param)
TRIVIAL constexpr auto rank(Arena<> alloc, PtrMatrix<std::int64_t> A)
  -> std::ptrdiff_t {
  MutDensePtrMatrix<std::int64_t> B = matrix<std::int64_t>(&alloc, shape(A));
  return std::ptrdiff_t(simplifySystemImpl(B << A, 0));
}
template <MatrixDimension S0, MatrixDimension S1>
TRIVIAL constexpr void simplifySystem(MutArray<std::int64_t, S0> &A,
                                      MutArray<std::int64_t, S1> &B) {
  simplifySystemsImpl({A, B});
  if (Row newM = numNonZeroRows(A); newM < A.numRow()) {
    A.truncate(newM);
    B.truncate(newM);
  }
}
/// A is the matrix we factorize, `U` is initially uninitialized, but is
/// destination of the unimodular matrix.
TRIVIAL constexpr void hermite(MutPtrMatrix<std::int64_t> A,
                               MutSquarePtrMatrix<std::int64_t> U) {
  invariant(A.numRow() == U.numRow());
  U << 0;
  U.diag() << 1;
  simplifySystemsImpl({A, U});
}

#ifndef POLYMATHNOEXPLICITSIMDARRAY
/// use A[j,k] to zero A[i,k]
TRIVIAL constexpr auto zeroWithRowOp(MutPtrMatrix<std::int64_t> A, Row<> i,
                                     Row<> j, Col<> k, std::int64_t f)
  -> std::int64_t {
  std::int64_t Aik = A[i, k];
  if (!Aik) return f;
  std::int64_t Ajk = A[j, k];
  invariant(Ajk != 0);
  std::int64_t g = gcd(Aik, Ajk);
  if (g != 1) {
    Aik /= g;
    Ajk /= g;
  }
  std::int64_t ret = f * Ajk;
  constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
  simd::Vec<W, std::int64_t> vAjk = simd::vbroadcast<W, std::int64_t>(Ajk),
                             vAik = simd::vbroadcast<W, std::int64_t>(Aik),
                             vg = {ret};
  static constexpr simd::Vec<W, std::int64_t> one =
    simd::Vec<W, std::int64_t>{} + 1;
  PtrMatrix<std::int64_t> B = A; // const ref
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol()), l = 0;
  if (ret != 1) {
    for (;;) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      simd::Vec<W, std::int64_t> Ail =
        vAjk * B[i, u].vec_ - vAik * B[j, u].vec_;
      A[i, u] = Ail;
      vg = gcd<W>(Ail, vg);
      l += W;
      // if none are `> 1`, we break and take the route that skips gcd
      if (!bool(simd::cmp::gt<W, std::int64_t>(vg, one))) break;
    }
  }
  if (l < L) {
    // the gcd is 1
    for (;; l += W) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      A[i, u] = vAjk * B[i, u].vec_ - vAik * B[j, u].vec_;
    }
  } else if (simd::cmp::gt<W, std::int64_t>(vg, one)) {
    // gcd isn't one, so we can scale
    g = gcdreduce<W>(vg);
    if (g > 1) {
      for (std::ptrdiff_t ll = 0; ll < L; ++ll)
        if (std::int64_t Ail = A[i, ll]) A[i, ll] = Ail / g;
      std::int64_t r = ret / g;
      invariant(r * g, ret);
      ret = r;
    }
  }
  return ret;
}
TRIVIAL constexpr void zeroWithRowOp(MutPtrMatrix<std::int64_t> A, Row<> i,
                                     Row<> j, Col<> k) {
  std::int64_t Aik = A[i, k];
  if (!Aik) return;
  std::int64_t Ajk = A[j, k];
  invariant(Ajk != 0);
  std::int64_t g = gcd(Aik, Ajk);
  if (g != 1) {
    Aik /= g;
    Ajk /= g;
  }
  constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
  static constexpr simd::Vec<W, std::int64_t> one =
    simd::Vec<W, std::int64_t>{} + 1;
  simd::Vec<W, std::int64_t> vAjk = simd::vbroadcast<W, std::int64_t>(Ajk),
                             vAik = simd::vbroadcast<W, std::int64_t>(Aik),
                             vg = {Ajk};
  PtrMatrix<std::int64_t> B = A; // const ref
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol()), l = 0;
  if (Ajk != 1) {
    for (;;) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      simd::Vec<W, std::int64_t> Ail =
        vAjk * B[i, u].vec_ - vAik * B[j, u].vec_;
      A[i, u] = Ail;
      vg = gcd<W>(Ail, vg);
      l += W;
      // if none are `> 1`, we break and take the route that skips gcd
      if (!bool(simd::cmp::gt<W, std::int64_t>(vg, one))) break;
    }
  }
  if (l < L) {
    // requires we did not execute above branch, the gcd is 1
    for (;; l += W) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      A[i, u] = vAjk * B[i, u].vec_ - vAik * B[j, u].vec_;
    }
  } else if (simd::cmp::gt<W, std::int64_t>(vg, one)) {
    // gcd isn't one, so we can scale
    if (g = gcdreduce<W>(vg); g > 1) {
      for (std::ptrdiff_t ll = 0; ll < L; ++ll)
        if (std::int64_t Ail = A[i, ll]) A[i, ll] = Ail / g;
    }
  }
}
#else
/// use A[j,k] to zero A[i,k]
TRIVIAL constexpr auto zeroWithRowOp(MutPtrMatrix<std::int64_t> A, Row<> i,
                                     Row<> j, Col<> k, std::int64_t f)
  -> std::int64_t {
  std::int64_t Aik = A[i, k];
  if (!Aik) return f;
  std::int64_t Ajk = A[j, k];
  invariant(Ajk != 0);
  std::int64_t g = gcd(Aik, Ajk);
  Aik /= g;
  Ajk /= g;
  std::int64_t ret = f * Ajk, vg = ret;
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol()), l = 0;
  for (; (vg != 1) && (l < L); ++l) {
    std::int64_t Ail = A[i, l] = Ajk * A[i, l] - Aik * A[j, l];
    vg = gcd(Ail, vg);
  }
  if (l < L) {
    for (; l < L; ++l) A[i, l] = Ajk * A[i, l] - Aik * A[j, l];
  } else if (vg != 1) {
    for (std::ptrdiff_t ll = 0; ll < L; ++ll)
      if (std::int64_t Ail = A[i, ll]) A[i, ll] = Ail / vg;
    std::int64_t r = ret / vg;
    invariant(r * vg, ret);
    ret = r;
  }
  return ret;
}
TRIVIAL constexpr void zeroWithRowOp(MutPtrMatrix<std::int64_t> A, Row<> i,
                                     Row<> j, Col<> k) {
  std::int64_t Aik = A[i, k];
  if (!Aik) return;
  std::int64_t Ajk = A[j, k];
  invariant(Ajk != 0);
  std::int64_t g = gcd(Aik, Ajk), vg = 0;
  Aik /= g;
  Ajk /= g;
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol()), l = 0;
  for (; (vg != 1) && (l < L); ++l) {
    std::int64_t Ail = A[i, l] = Ajk * A[i, l] - Aik * A[j, l];
    vg = gcd(Ail, vg);
  }
  if (l < L) {
    for (; l < L; ++l) A[i, l] = Ajk * A[i, l] - Aik * A[j, l];
  } else if (vg != 1) {
    for (std::ptrdiff_t ll = 0; ll < L; ++ll)
      if (std::int64_t Ail = A[i, ll]) A[i, ll] = Ail / vg;
  }
  return;
}
#endif
TRIVIAL constexpr void bareiss(MutPtrMatrix<std::int64_t> A,
                               MutPtrVector<std::ptrdiff_t> pivots) {
  const auto [M, N] = shape(A);
  invariant(pivots.size(), std::min(M, N));
  std::int64_t prev = 1, piv_ind = 0;
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c) {
    if (std::optional<std::ptrdiff_t> piv =
          detail::pivotRowsBareiss(A, c, row(M), row(r))) {
      pivots[piv_ind++] = *piv;
      auto j{_(c + 1, N)};
      std::int64_t Arc = A[r, c];
      for (std::ptrdiff_t k = r + 1; k < M; ++k) {
        A[k, j] << (Arc * A[k, j] - A[k, c] * A[r, j]) / prev;
        A[k, r] = 0;
      }
      prev = A[r++, c];
    }
  }
}

TRIVIAL [[nodiscard]] constexpr auto bareiss(IntMatrix<> &A)
  -> Vector<std::ptrdiff_t> {
  Vector<std::ptrdiff_t> pivots(length(A.minRowCol()));
  bareiss(A, pivots);
  return pivots;
}

TRIVIAL constexpr void solveColumn(MutPtrMatrix<std::int64_t> A,
                                   MutPtrMatrix<std::int64_t> B,
                                   std::ptrdiff_t r, std::ptrdiff_t c) {
  const auto [M, N] = shape(A);
  utils::assume(B.numRow() == M);
  utils::invariant(r < M);
  utils::invariant(c < N);
  utils::invariant(c < B.numCol());
  utils::invariant(!detail::pivotRowsPair({A, B}, col(c), row(M), row(r)));
  detail::zeroColumnPair({A, B}, col(c), row(r));
}
/// void solveSystem(IntMatrix &A, IntMatrix &B)
/// Say we wanted to solve \f$\textbf{AX} = \textbf{B}\f$.
/// `solveSystem` left-multiplies both sides by
/// a matrix \f$\textbf{W}\f$ that diagonalizes \f$\textbf{A}\f$.
/// Once \f$\textbf{A}\f$ has been diagonalized, the solution is trivial.
/// Both inputs are overwritten with the product of the left multiplications.
TRIVIAL constexpr void solveSystem(MutPtrMatrix<std::int64_t> A,
                                   MutPtrMatrix<std::int64_t> B) {
  const auto [M, N] = shape(A);
  utils::assume(B.numRow() == M);
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c)
    if (!detail::pivotRowsPair({A, B}, col(c), row(M), row(r)))
      detail::zeroColumnPair({A, B}, col(c), row(r++));
}
// like solveSystem, except it right-multiplies.
// That is, given `XA = B`, it right-multiplies both sides by
// a matrix to diagonalize `A`.
TRIVIAL constexpr void solveSystemRight(MutPtrMatrix<std::int64_t> A,
                                        MutPtrMatrix<std::int64_t> B) {
  const auto [M, N] = shape(A);
  utils::assume(B.numCol() == N);
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++r)
    if (!detail::pivotColsPair({A, B}, row(r), col(N), col(c)))
      detail::zeroColumnPair({A, B}, row(r), col(c++));
}
// diagonalizes A(0:K,0:K)
TRIVIAL constexpr void solveSystem(MutPtrMatrix<std::int64_t> A,
                                   std::ptrdiff_t K) {
  Row M = A.numRow();
  for (std::ptrdiff_t r = 0, c = 0; c < K && r < M; ++c)
    if (!pivotRows(A, col(c), M, row(r)))
      detail::zeroColumn(A, col(c), row(r++));
}
// diagonalizes A(0:K,1:K+1)
TRIVIAL constexpr void solveSystemSkip(MutPtrMatrix<std::int64_t> A) {
  const auto [M, N] = shape(A);
  for (std::ptrdiff_t r = 0, c = 1; c < N && r < M; ++c)
    if (!pivotRows(A, col(c), row(M), row(r)))
      detail::zeroColumn(A, col(c), row(r++));
}

// returns `true` if the solve failed, `false` otherwise
// diagonals contain denominators.
// Assumes the last column is the vector to solve for.
TRIVIAL constexpr void solveSystem(MutPtrMatrix<std::int64_t> A) {
  solveSystem(A, std::ptrdiff_t(A.numCol()) - 1);
}

// /// bareissSolveSystem(A, B)
// /// Like solveSystem, but uses Bareiss algorithm instead of Hermite Normal Form
// /// for the initial triangular form. Operations done to A are also applied to B.
// /// Avoids allocating temporary arrays like pivot vectors.
// constexpr void bareissSolveSystem(MutPtrMatrix<std::int64_t> A,
//                                   MutPtrMatrix<std::int64_t> B) {
//   const auto [M, N] = shape(A);
//   utils::assume(B.numRow() == M);

//   // Step 1: Apply Bareiss algorithm to A, mirroring operations on B
//   std::int64_t prev = 1;
//   for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c) {
//     // Find pivot row without allocating array
//     std::ptrdiff_t piv_row = r;
//     while (piv_row < M && A[piv_row, c] == 0) ++piv_row;
//     if (piv_row == M) continue;
//     // Swap rows if needed (apply to both A and B)
//     if (r != piv_row) {
//       swap(A, row(r), row(piv_row));
//       swap(B, row(r), row(piv_row));
//     }

//     // Bareiss elimination: eliminate below pivot
//     std::int64_t Arc = A[r, c];
//     auto j{_(c + 1, N)};
//     for (std::ptrdiff_t k = r + 1; k < M; ++k) {
//       if (std::int64_t Akc = A[k, c]) {
//         if (prev == 1) {
//           // Apply Bareiss to A
//           A[k, j] << Arc * A[k, j] - Akc * A[r, j];
//           B[k, _] << Arc * B[k, _] - Akc * B[r, _];
//         } else {
//           A[k, j] << (Arc * A[k, j] - Akc * A[r, j]) / prev;
//           B[k, _] << (Arc * B[k, _] - Akc * B[r, _]) / prev;
//         }
//         A[k, c] = 0;
//       } else if (Arc != prev) {
//         A[k, j] << (Arc * A[k, j]) / prev;
//         B[k, _] << (Arc * B[k, _]) / prev;
//       }
//     }
//     prev = Arc;
//     ++r;
//   }
//   // perform a triangular solve
//   for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c) {
//     std::int64_t Arc = A[r, c];
//     if (!Arc) continue;
//     // zero rows < r; d is "diagonal"'s column
//     for (std::ptrdiff_t z = 0, d = 0; z < r; ++z) {
//       std::int64_t Azc = A[z, c];
//       if (!Azc) continue;
//       const auto [p, q, Arcr, Ajcr] = detail::gcdxScale(Arc, Azc);
//       MutPtrVector<std::int64_t> Ar{A[r, _(c + 1, end)]},
//         Az{A[z, _(c + 1, end)]};
//       // MutPtrVector<std::int64_t> Ar{A[r, _(c + 1, end)]},
//       //   Az{A[z, _(c + 1, end)]};
//       Az << Arcr * Az - Ajcr * Ar;
//       MutPtrVector<std::int64_t> Br{B[r, _]}, Bz{B[z, _]};
//       Bz << Arcr * Bz - Ajcr * Br;
//       A[z, c] = 0;
//       // Arc = q * Azc + p * Arc;
//       if (Arcr == 1) continue;
//       // update "diagonal"
//       int64_t Azd = A[z, d];
//       while (!Azd) Azd = A[z, ++d];
//       A[z, d++] = Arcr * Azd;
//     }
//     ++r;
//   }

//   // // Step 2: Back-solve the upper triangular system
//   // // Now A is upper triangular, solve A * X = B by back substitution
//   // // X = A \ B, overwrite B with X
//   // // B[i,j] = A[i,_(i,end)] * X[_(i,end),j].t()
//   // // B[i,j] = A[i,i] * X[i,j] + A[i,_(i+1,end)] * X[_(i+1,end),j].t()
//   // // X[i,j] = (B[i,j] - A[i,_(i+1,end)] * X[_(i+1,end),j].t()) / A[i,i]
//   // // Problem:
//   // std::ptrdiff_t rank = std::min(M, N);
//   // for (std::ptrdiff_t i = rank - 1; i >= 0; --i) {
//   //   if (A[i, i] != 0) {
//   //     // Eliminate entries above the diagonal
//   //     for (std::ptrdiff_t j = 0; j < i; ++j) {
//   //       if (A[j, i] != 0) {
//   //         std::int64_t Aii = A[i, i], Aji = A[j, i];
//   //         // Use integer arithmetic: B[j, :] = Aii * B[j, :] - Aji * B[i, :]
//   //         for (std::ptrdiff_t k = 0; k < B.numCol(); ++k)
//   //           B[j, k] = Aii * B[j, k] - Aji * B[i, k];
//   //         for (std::ptrdiff_t k = 0; k < N; ++k)
//   //           A[j, k] = Aii * A[j, k] - Aji * A[i, k];
//   //       }
//   //     }
//   //   }
//   // }
// }

// /// bareissSolveSystemRight(A, B)
// /// Like solveSystemRight, but uses Bareiss algorithm for the initial triangular
// /// form. Given XA = B, right-multiplies both sides to diagonalize A.
// TRIVIAL constexpr void bareissSolveSystemRight(MutPtrMatrix<std::int64_t> A,
//                                                MutPtrMatrix<std::int64_t> B) {
//   const auto [M, N] = shape(A);
//   utils::assume(B.numCol() == N);

//   // Step 1: Apply Bareiss algorithm to A, mirroring operations on B
//   // For right multiplication, we work on columns to get lower triangular form
//   std::int64_t prev = 1;
//   for (std::ptrdiff_t c = 0, r = 0; r < M && c < N; ++r) {
//     // Find pivot column without allocating array
//     std::ptrdiff_t piv_col = c;
//     while (piv_col < N && A[r, piv_col] == 0) ++piv_col;

//     if (piv_col < N) {
//       // Swap columns if needed (apply to both A and B)
//       if (c != piv_col) {
//         swap(A, col(c), col(piv_col));
//         swap(B, col(c), col(piv_col));
//       }

//       // Bareiss elimination: eliminate to the right of pivot
//       std::int64_t Arc = A[r, c];
//       for (std::ptrdiff_t k = c + 1; k < N; ++k) {
//         if (std::int64_t Ark = A[r, k]) {
//           // Apply Bareiss to A
//           for (std::ptrdiff_t i = r + 1; i < M; ++i)
//             A[i, k] = (Arc * A[i, k] - Ark * A[i, c]) / prev;
//           A[r, k] = 0;

//           // Apply same operations to B
//           for (std::ptrdiff_t i = 0; i < B.numRow(); ++i)
//             B[i, k] = (Arc * B[i, k] - Ark * B[i, c]) / prev;
//         }
//       }
//       prev = Arc;
//       ++c;
//     }
//   }

//   // Step 2: Back-solve the triangular system
//   // Now A has pivots on the diagonal, solve X * A = B by forward substitution
//   std::ptrdiff_t rank = std::min(M, N);
//   for (std::ptrdiff_t i = 0; i < rank; ++i) {
//     if (A[i, i] != 0) {
//       // Eliminate entries to the left of diagonal
//       for (std::ptrdiff_t j = i + 1; j < rank; ++j) {
//         if (A[i, j] != 0) {
//           std::int64_t Aii = A[i, i], Aij = A[i, j];
//           // Use integer arithmetic: B[:, j] = Aii * B[:, j] - Aij * B[:, i]
//           for (std::ptrdiff_t k = 0; k < B.numRow(); ++k)
//             B[k, j] = Aii * B[k, j] - Aij * B[k, i];
//           for (std::ptrdiff_t k = 0; k < M; ++k)
//             A[k, j] = Aii * A[k, j] - Aij * A[k, i];
//         }
//       }
//     }
//   }
// }

/// inv(A) -> (D, B)
/// Given a matrix \f$\textbf{A}\f$, returns two matrices \f$\textbf{D}\f$ and
/// \f$\textbf{B}\f$ so that \f$\textbf{D}^{-1}\textbf{B} = \textbf{A}^{-1}\f$,
/// and \f$\textbf{D}\f$ is diagonal.
/// NOTE: This function assumes non-singular
/// Mutates `A`
// NOLINTNEXTLINE(performance-unnecessary-value-param)
TRIVIAL [[nodiscard]] constexpr auto inv(Arena<> *alloc,
                                         MutSquarePtrMatrix<std::int64_t> A)
  -> MutSquarePtrMatrix<std::int64_t> {
  MutSquarePtrMatrix<std::int64_t> B =
    identity<std::int64_t>(alloc, std::ptrdiff_t(A.numCol()));
  solveSystem(A, B);
  return B;
}
// scaledInv(A, B) -> s
// reads and writes A, writes B
// B = s * inv(A) // A is diagonalized in the process
TRIVIAL [[nodiscard]] constexpr auto
scaledInv(MutSquarePtrMatrix<std::int64_t> A,
          MutSquarePtrMatrix<std::int64_t> B) -> std::int64_t {
  B.zero();
  B.diag() << 1;
  solveSystem(A, B);
  auto [s, nonUnity] = lcmNonUnity(A.diag());
  if (nonUnity)
    for (std::ptrdiff_t i = 0; i < A.numRow(); ++i) B[i, _] *= s / A[i, i];
  return s;
}
/// inv(A) -> (B, s)
/// Given a matrix \f$\textbf{A}\f$, returns a matrix \f$\textbf{B}\f$ and a
/// scalar \f$s\f$ such that \f$\frac{1}{s}\textbf{B} = \textbf{A}^{-1}\f$.
/// NOTE: This function assumes non-singular
/// D0 * B^{-1} = Binv0
/// (s/s) * D0 * B^{-1} = Binv0
/// s * B^{-1} = (s/D0) * Binv0
/// mutates `A`
TRIVIAL [[nodiscard]] constexpr auto
scaledInv(Arena<> *alloc, MutSquarePtrMatrix<std::int64_t> A)
  -> containers::Pair<MutSquarePtrMatrix<std::int64_t>, std::int64_t> {
  MutSquarePtrMatrix<std::int64_t> B =
    square_matrix<std::int64_t>(alloc, std::ptrdiff_t(A.numCol()));
  return {B, scaledInv(A, B)};
}

TRIVIAL constexpr auto nullSpace(MutDensePtrMatrix<std::int64_t> B,
                                 MutDensePtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  B << 0;
  B.diag() << 1;
  solveSystem(A, B);
  Row R = numNonZeroRows(A); // rank
  // we discard first R columns
  return B[_(R, end), _];
}
// one row per null dim
TRIVIAL constexpr auto nullSpace(Arena<> *alloc,
                                 MutDensePtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  Row M = A.numRow();
  return nullSpace(matrix<std::int64_t>(alloc, M, math::col(std::ptrdiff_t(M))),
                   A);
}

// FIXME: why do we have two?
TRIVIAL constexpr auto orthogonalize(IntMatrix<> A)
  -> containers::Pair<SquareMatrix<std::int64_t>, Vector<unsigned>> {
  return detail::orthogonalizeBang(A);
}

} // namespace NormalForm

// FIXME: why do we have two?
TRIVIAL constexpr auto orthogonalize(Arena<> *alloc,
                                     MutDensePtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  if ((A.numCol() < 2) || (A.numRow() == 0)) return A;
  normalizeByGCD(A[0, _]);
  if (A.numRow() == 1) return A;
  auto s = alloc->scope();
  MutPtrVector<Rational> buff{
    vector<Rational>(alloc, std::ptrdiff_t(A.numCol()))};
  // Vector<Rational, 8> buff;
  // buff.resizeForOverwrite(std::ptrdiff_t(A.numCol()));
  std::ptrdiff_t offset = 0;
  for (std::ptrdiff_t i = 1; i < A.numRow(); ++i) {
    buff << A[i, _];
    for (std::ptrdiff_t j = 0; j < i - offset; ++j) {
      std::int64_t n = 0;
      std::int64_t d = 0;
      for (std::ptrdiff_t k = 0; k < A.numCol(); ++k) {
        n += A[i, k] * A[j, k];
        d += A[j, k] * A[j, k];
      }
      for (std::ptrdiff_t k = 0; k < A.numCol(); ++k)
        buff[k] -= Rational::createPositiveDenominator(A[j, k] * n, d);
    }
    bool all_zero = true;
    for (std::ptrdiff_t k = 0; k < A.numCol(); ++k) {
      if (buff[k].numerator_) {
        all_zero = false;
        break;
      }
    }
    if (all_zero) {
      ++offset;
      continue;
    }
    std::int64_t lm = 1;
    for (std::ptrdiff_t k = 0; k < A.numCol(); ++k)
      lm = lcm(lm, buff[k].denominator_);
    for (std::ptrdiff_t k = 0; k < A.numCol(); ++k)
      A[i - offset, k] = buff[k].numerator_ * (lm / buff[k].denominator_);
  }
  return A[_(end - offset), _];
}

TRIVIAL [[nodiscard]] constexpr auto
orthogonalNullSpace(Arena<> *alloc, MutDensePtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  return orthogonalize(alloc, NormalForm::nullSpace(alloc, A));
}

} // namespace math

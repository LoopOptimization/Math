module;
#include "Macros.hxx"
module NormalForm;
// Implementation unit for NormalForm module

using namespace math;
using containers::Pair, containers::Tuple;

namespace math::NormalForm {
auto gcdxScale(std::int64_t a, std::int64_t b) -> std::array<std::int64_t, 4> {
  if (constexpr_abs(a) == 1) return {a, 0, a, b};
  auto [g, p, q, adg, bdg] = dgcdx(a, b);
  return {p, q, adg, bdg};
}
} // namespace math::NormalForm

namespace {
// zero out below diagonal
void zeroSupDiagonal(MutPtrMatrix<std::int64_t> A,
                     MutSquarePtrMatrix<std::int64_t> K, std::ptrdiff_t i,
                     Row<> M, Col<> N) {
  std::ptrdiff_t minMN = std::min(std::ptrdiff_t(M), std::ptrdiff_t(N));
  for (std::ptrdiff_t j = i + 1; j < M; ++j) {
    std::int64_t Aii = A[i, i];
    if (std::int64_t Aji = A[j, i]) {
      const auto [p, q, Aiir, Aijr] = NormalForm::gcdxScale(Aii, Aji);

      {
        MutPtrVector<std::int64_t> Ai{A[i, _(0, minMN)]}, Aj{A[j, _(0, minMN)]},
          Ki{K[i, _(0, minMN)]}, Kj{K[j, _(0, minMN)]};
        Tuple(Ai, Aj, Ki, Kj) << Tuple(p * Ai + q * Aj, Aiir * Aj - Aijr * Ai,
                                       p * Ki + q * Kj, Aiir * Kj - Aijr * Ki);
      }
      if (std::ptrdiff_t(M) > std::ptrdiff_t(N)) {
        MutPtrVector<std::int64_t> Ki{K[i, _(N, M)]}, Kj{K[j, _(N, M)]};
        Pair(Ki, Kj) << Tuple(p * Ki + q * Kj, Aiir * Kj - Aijr * Ki);
      } else if (std::ptrdiff_t(N) > std::ptrdiff_t(M)) {
        MutPtrVector<std::int64_t> Ai{A[i, _(M, N)]}, Aj{A[j, _(M, N)]};
        Pair(Ai, Aj) << Tuple(p * Ai + q * Aj, Aiir * Aj - Aijr * Ai);
      }
    }
  }
}
// This method is only called by orthogonalize, hence we can assume
// (Akk == 1) || (Akk == -1)
void zeroSubDiagonal(MutPtrMatrix<std::int64_t> A,
                     MutSquarePtrMatrix<std::int64_t> K, std::ptrdiff_t k,
                     Row<> M, Col<> N) {
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
  }
}

auto pivotRowsPair(std::array<MutPtrMatrix<std::int64_t>, 2> AK, Col<> i,
                   Row<> M, Row<> piv) -> bool {
  Row j = piv;
  while (AK[0][piv, i] == 0)
    if (++piv == M) return true;
  if (j != piv) {
    math::swap(AK[0], j, piv);
    math::swap(AK[1], j, piv);
  }
  return false;
}
auto pivotColsPair(std::array<MutPtrMatrix<std::int64_t>, 2> AK, Row<> i,
                   Col<> N, Col<> piv) -> bool {
  Col j = piv;
  while (AK[0][i, piv] == 0)
    if (++piv == N) return true;
  if (j != piv) {
    math::swap(AK[0], j, piv);
    math::swap(AK[1], j, piv);
  }
  return false;
}
void dropCol(MutPtrMatrix<std::int64_t> A, std::ptrdiff_t i, Row<> M, Col<> N) {
  // if any rows are left, we shift them up to replace it
  if (N <= i) return;
  A[_(0, M), _(i, N)] << A[_(0, M), _(i, N) + 1];
}

void zeroSupDiagonal(MutPtrMatrix<std::int64_t> A, Col<> c, Row<> r) {
  auto [M, N] = shape(A);
  for (std::ptrdiff_t j = std::ptrdiff_t(r) + 1; j < M; ++j) {
    std::int64_t Aii = A[r, c];
    if (std::int64_t Aij = A[j, c]) {
      const auto [p, q, Aiir, Aijr] = NormalForm::gcdxScale(Aii, Aij);
      MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
      Pair(Ar, Aj) << Tuple(p * Ar + q * Aj, Aiir * Aj - Aijr * Ar);
    }
  }
}
void zeroSupDiagonal(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c,
                     Row<> r) {
  auto [A, B] = AB;
  auto [M, N] = shape(A);
  invariant(M, std::ptrdiff_t(B.numRow()));
  for (std::ptrdiff_t j = std::ptrdiff_t(r) + 1; j < M; ++j) {
    std::int64_t Aii = A[r, c], Aij = A[j, c];
    if (!Aij) continue;
    const auto [p, q, Aiir, Aijr] = NormalForm::gcdxScale(Aii, Aij);
    MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    Pair(Ar, Aj) << Tuple(p * Ar + q * Aj, Aiir * Aj - Aijr * Ar);
    MutPtrVector<std::int64_t> Br = B[r, _], Bj = B[j, _];
    Pair(Br, Bj) << Tuple(p * Br + q * Bj, Aiir * Bj - Aijr * Br);
  }
}
void reduceSubDiagonal(MutPtrMatrix<std::int64_t> A, Col<> c, Row<> r) {
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
void reduceSubDiagonalStack(MutPtrMatrix<std::int64_t> A,
                            MutPtrMatrix<std::int64_t> B, std::ptrdiff_t c,
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
void reduceSubDiagonal(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c,
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

void reduceColumn(MutPtrMatrix<std::int64_t> A, Col<> c, Row<> r) {
  zeroSupDiagonal(A, c, r);
  reduceSubDiagonal(A, c, r);
}

// NormalForm version assumes zero rows are sorted to end due to pivoting
void removeZeroRows(MutDensePtrMatrix<std::int64_t> &A) {
  A.truncate(NormalForm::numNonZeroRows(A));
}

void reduceColumn(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c,
                  Row<> r) {
  zeroSupDiagonal(AB, c, r);
  reduceSubDiagonal(AB, c, r);
}
/// multiplies `A` and `B` by matrix `X`, where `X` reduces `A` to a normal
/// form.
void zeroWithRowOperation(MutPtrMatrix<std::int64_t> A, Row<> i, Row<> j,
                          Col<> k, Range<std::ptrdiff_t, std::ptrdiff_t> skip) {
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
void zeroColumnPair(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Col<> c,
                    Row<> r) {
  auto [A, B] = AB;
  const Row M = A.numRow();
  invariant(M, B.numRow());
  std::int64_t Arc = A[r, c];
  for (std::ptrdiff_t j = 0; j < r; ++j) {
    std::int64_t Ajc = A[j, c];
    if (!Ajc) continue;
    std::int64_t g = gcd(Arc, Ajc), x = Arc / g, y = Ajc / g;
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutPtrVector<std::int64_t> Ar = AB[i][r, _], Aj = AB[i][j, _];
      Aj << x * Aj - y * Ar;
    }
  }
  // greater rows in previous columns have been zeroed out
  // therefore it is safe to use them for row operations with this row
  for (auto j = std::ptrdiff_t(r); ++j < M;) {
    std::int64_t Ajc = A[j, c];
    if (!Ajc) continue;
    const auto [p, q, Arcr, Ajcr] = NormalForm::gcdxScale(Arc, Ajc);
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutPtrVector<std::int64_t> Ar = AB[i][r, _], Aj = AB[i][j, _];
      Pair(Ar, Aj) << Tuple(q * Aj + p * Ar, Arcr * Aj - Ajcr * Ar);
    }
    Arc = A[r, c];
  }
}
// use col `c` to zero the remaining cols of row `r`
void zeroColumnPair(std::array<MutPtrMatrix<std::int64_t>, 2> AB, Row<> r,
                    Col<> c) {
  auto [A, B] = AB;
  const Col N = A.numCol();
  invariant(N, B.numCol());
  std::int64_t Arc = A[r, c];
  for (std::ptrdiff_t j = 0; j < c; ++j) {
    std::int64_t Arj = A[r, j];
    if (!Arj) continue;
    auto [x, y] = divgcd(Arc, Arj);
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutArray<std::int64_t, StridedRange<>> Ac = AB[i][_, c], Aj = AB[i][_, j];
      Aj << x * Aj - y * Ac;
    }
  }
  // greater rows in previous columns have been zeroed out
  // therefore it is safe to use them for row operations with this row
  for (auto j = std::ptrdiff_t(c); ++j < N;) {
    std::int64_t Arj = A[r, j];
    if (!Arj) continue;
    const auto [p, q, Arcr, Ajcr] = NormalForm::gcdxScale(Arc, Arj);
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
      MutArray<std::int64_t, StridedRange<>> Ac = AB[i][_, c], Aj = AB[i][_, j];
      Pair(Ac, Aj) << Tuple(q * Aj + p * Ac, Arcr * Aj - Ajcr * Ac);
    }
    Arc = A[r, c]; // reload in case it changed
  }
}
// use row `r` to zero the remaining rows of column `c`
void zeroColumn(MutPtrMatrix<std::int64_t> A, Col<> c, Row<> r) {
  const Row M = A.numRow();
  std::int64_t Arc = A[r, c];
  for (std::ptrdiff_t j = 0; j < r; ++j) {
    std::int64_t Ajc = A[j, c];
    invariant(Arc != std::numeric_limits<std::int64_t>::min());
    invariant(Ajc != std::numeric_limits<std::int64_t>::min());
    if (!Ajc) continue;
    std::int64_t g = gcd(Arc, Ajc);
    A[j, _] << (Arc / g) * A[j, _] - (Ajc / g) * A[r, _];
  }
  // greater rows in previous columns have been zeroed out
  // therefore it is safe to use them for row operations with this row
  for (std::ptrdiff_t j = std::ptrdiff_t(r); ++j < M;) {
    std::int64_t Ajc = A[j, c];
    if (!Ajc) continue;
    const auto [p, q, Arcr, Ajcr] = NormalForm::gcdxScale(Arc, Ajc);
    MutPtrVector<std::int64_t> Ar = A[r, _], Aj = A[j, _];
    Pair(Ar, Aj) << Pair(q * Aj + p * Ar, Arcr * Aj - Ajcr * Ar);
    Arc = A[r, c]; // reload in case it changed
  }
}

auto pivotRowsBareiss(MutPtrMatrix<std::int64_t> A, std::ptrdiff_t i, Row<> M,
                      Row<> piv) -> std::optional<std::ptrdiff_t> {
  Row j = piv;
  while (A[piv, i] == 0)
    if (++piv == M) return {};
  if (j != piv) swap(A, j, piv);
  return std::ptrdiff_t(piv);
}

auto orthogonalizeBang(MutDensePtrMatrix<std::int64_t> &A)
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
auto nullSpace(MutDensePtrMatrix<std::int64_t> B, MutPtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  B << 0;
  B.diag() << 1;
  math::NormalForm::solveSystem(A, B);
  Row R = math::NormalForm::numNonZeroRows(A); // rank
  // we discard first R columns
  return B[_(R, end), _];
}

} // namespace
namespace math {
auto gcd(PtrVector<std::int64_t> x) -> std::int64_t {
  const std::ptrdiff_t N = x.size();
  if (!N) return 0;
  std::int64_t g = constexpr_abs(x[0]);
  for (std::ptrdiff_t n = 1; (n < N) & (g != 1); ++n) g = gcd(g, x[n]);
  return g;
}
void normalizeByGCD(MutPtrVector<std::int64_t> x) {
  std::ptrdiff_t N = x.size();
  switch (N) {
  case 0: return;
  case 1: x[0] = 1; return;
  default:
    std::int64_t g = gcd(x[0], x[1]);
    for (std::ptrdiff_t n = 2; (n < N) & (g != 1); ++n) g = gcd(g, x[n]);
    if (g > 1) x /= g;
  }
}

namespace NormalForm {
using alloc::Arena;

auto pivotRows(MutPtrMatrix<std::int64_t> A, MutSquarePtrMatrix<std::int64_t> K,
               std::ptrdiff_t i, Row<> M) -> bool {
  MutPtrMatrix<std::int64_t> B = K;
  return ::pivotRowsPair({A, B}, col(i), M, row(i));
}
auto pivotRows(MutPtrMatrix<std::int64_t> A, Col<> i, Row<> M, Row<> piv)
  -> bool {
  Row j = piv;
  while (A[piv, i] == 0)
    if (++piv == std::ptrdiff_t(M)) return true;
  if (j != piv) swap(A, j, piv);
  return false;
}
auto pivotRows(MutPtrMatrix<std::int64_t> A, std::ptrdiff_t i, Row<> N)
  -> bool {
  return pivotRows(A, col(i), N, row(i));
}
/// numNonZeroRows(PtrMatrix<std::int64_t> A) -> Row
/// Assumes some number of the trailing rows have been
/// zeroed out.  Returns the number of rows that are remaining.
auto numNonZeroRows(PtrMatrix<std::int64_t> A) -> Row<> {
  Row newM = A.numRow();
  while (newM && allZero(A[std::ptrdiff_t(newM) - 1, _])) --newM;
  return newM;
}
auto updateForNewRow(MutPtrMatrix<std::int64_t> A) -> std::ptrdiff_t {
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
void simplifySystemsImpl(std::array<MutPtrMatrix<std::int64_t>, 2> AB) {
  auto [M, N] = shape(AB[0]);
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c)
    if (!::pivotRowsPair(AB, col(c), row(M), row(r)))
      ::reduceColumn(AB, col(c), row(r++));
}

// treats A as stacked on top of B
void reduceColumnStack(MutPtrMatrix<std::int64_t> A,
                       MutPtrMatrix<std::int64_t> B, std::ptrdiff_t c,
                       std::ptrdiff_t r) {
  ::zeroSupDiagonal(B, col(c), row(r));
  ::reduceSubDiagonalStack(B, A, c, r);
}
// pass by value, returns number of rows to truncate
auto simplifySystemImpl(MutPtrMatrix<std::int64_t> A, std::ptrdiff_t colInit)
  -> Row<> {
  auto [M, N] = shape(A);
  for (std::ptrdiff_t r = 0, c = colInit; c < N && r < M; ++c)
    if (!pivotRows(A, col(c), row(M), row(r)))
      ::reduceColumn(A, col(c), row(r++));
  return numNonZeroRows(A);
}

void simplifySystem(MutPtrMatrix<std::int64_t> &E, std::ptrdiff_t colInit) {
  E.truncate(simplifySystemImpl(E, colInit));
}
// TODO: `const IntMatrix &` can be copied to `MutPtrMatrix<std::int64_t>`
// this happens via `const IntMatrix &` -> `const MutPtrMatrix<std::int64_t> &`
// -> `MutPtrMatrix<std::int64_t>`. Perhaps we should define `MutPtrMatrix(const
// MutPtrMatrix &) = delete;`?
//
// NOLINTNEXTLINE(performance-unnecessary-value-param)
auto rank(Arena<> alloc, PtrMatrix<std::int64_t> A) -> std::ptrdiff_t {
  MutDensePtrMatrix<std::int64_t> B = matrix<std::int64_t>(&alloc, shape(A));
  return std::ptrdiff_t(simplifySystemImpl(B << A, 0));
}
/// A is the matrix we factorize, `U` is initially uninitialized, but is
/// destination of the unimodular matrix.
void hermite(MutPtrMatrix<std::int64_t> A, MutSquarePtrMatrix<std::int64_t> U) {
  invariant(A.numRow() == U.numRow());
  U << 0;
  U.diag() << 1;
  simplifySystemsImpl({A, U});
}

// SIMD optimized functions
/// use A[j,k] to zero A[i,k]
auto zeroWithRowOp(MutPtrMatrix<std::int64_t> A, Row<> i, Row<> j, Col<> k,
                   std::int64_t f) -> std::int64_t {
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
  static constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
  simd::Vec<W, std::int64_t> vAjk = simd::vbroadcast<W, std::int64_t>(Ajk),
                             vAik = simd::vbroadcast<W, std::int64_t>(Aik),
                             vg = {ret};
  static simd::Vec<W, std::int64_t> one = simd::Vec<W, std::int64_t>{} + 1;
  PtrMatrix<std::int64_t> B = A; // const ref
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol()), l = 0;
  if (ret != 1) {
    for (;;) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      simd::Vec<W, std::int64_t> Ail =
        (vAjk * B[i, u].vec_) - (vAik * B[j, u].vec_);
      A[i, u] = Ail;
      vg = gcd<W>(Ail, vg);
      l += W;
      // if none are `> 1`, we break and take the route that skips gcd
      if (!simd::cmp::gt<W, std::int64_t>(vg, one).any()) break;
    }
  }
  if (l < L) {
    // the gcd is 1
    for (;; l += W) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      A[i, u] = (vAjk * B[i, u].vec_) - (vAik * B[j, u].vec_);
    }
  } else if (simd::cmp::gt<W, std::int64_t>(vg, one).any()) {
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
void zeroWithRowOp(MutPtrMatrix<std::int64_t> A, Row<> i, Row<> j, Col<> k) {
  std::int64_t Aik = A[i, k];
  if (!Aik) return;
  std::int64_t Ajk = A[j, k];
  invariant(Ajk != 0);
  std::int64_t g = gcd(Aik, Ajk);
  if (g != 1) {
    Aik /= g;
    Ajk /= g;
  }
  static constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
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
        (vAjk * B[i, u].vec_) - (vAik * B[j, u].vec_);
      A[i, u] = Ail;
      vg = gcd<W>(Ail, vg);
      l += W;
      // if none are `> 1`, we break and take the route that skips gcd
      if (!simd::cmp::gt<W, std::int64_t>(vg, one).any()) break;
    }
  }
  if (l < L) {
    // requires we did not execute above branch, the gcd is 1
    for (;; l += W) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      A[i, u] = (vAjk * B[i, u].vec_) - (vAik * B[j, u].vec_);
    }
  } else if (simd::cmp::gt<W, std::int64_t>(vg, one).any()) {
    // gcd isn't one, so we can scale
    if (g = gcdreduce<W>(vg); g > 1) {
      for (std::ptrdiff_t ll = 0; ll < L; ++ll)
        if (std::int64_t Ail = A[i, ll]) A[i, ll] = Ail / g;
    }
  }
}

// Batched versions for improved memory bandwidth

// Collect up to B eligible rows starting from start_row
// Returns {count, zero_idx}
// count = number of rows collected (0 to B)
// zero_idx = index of row 0 in batch, or -1 if not present
// start_row is updated to the next row to process
template <std::ptrdiff_t B>
auto collectBatch(PtrMatrix<std::int64_t> A, std::ptrdiff_t num_rows,
                  std::ptrdiff_t pivot_row, Col<> col_k,
                  std::ptrdiff_t &start_row,
                  std::array<std::ptrdiff_t, B> &rows,
                  std::array<std::int64_t, B> &coeffs)
  -> std::pair<std::ptrdiff_t, std::ptrdiff_t> {
  std::ptrdiff_t count = 0, zero_idx = -1;
  std::ptrdiff_t i = start_row;
  for (; i < num_rows && count < B; ++i) {
    if (i == pivot_row) continue;
    std::int64_t Aik = A[i, col_k];
    if (Aik == 0) continue;
    if (i == 0) zero_idx = count;
    rows[count] = i;
    coeffs[count] = Aik;
    ++count;
  }
  start_row = i;
  return {count, zero_idx};
}

// Batched elimination without factor tracking
template <std::ptrdiff_t B>
void zeroWithRowOp(MutPtrMatrix<std::int64_t> A,
                   std::array<std::ptrdiff_t, B> rows,
                   std::array<std::int64_t, B> coeffs, Row<> j, Col<> k) {
  static constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
  static_assert(B == W, "Batch size must equal SIMD width");
  static constexpr simd::Vec<W, std::int64_t> one =
    simd::Vec<W, std::int64_t>{} + 1;

  std::int64_t Ajk = A[j, k];
  invariant(Ajk != 0);

  // Pre-compute GCD-reduced coefficients using SIMD
  // Load coeffs using unaligned SIMD load (B == W by design)
  simd::Vec<W, std::int64_t> v_aik_coeff =
    simd::load(coeffs.data(), simd::mask::None<W>{});
  simd::Vec<W, std::int64_t> v_ajk_broadcast =
    simd::vbroadcast<W, std::int64_t>(Ajk);

  // Parallel GCD computation for all B pairs
  simd::Vec<W, std::int64_t> v_gcd = gcd<W>(v_aik_coeff, v_ajk_broadcast);

  // Parallel division with mask for g != 1
  auto needs_reduce = simd::cmp::ne<W, std::int64_t>(v_gcd, one);
  simd::Vec<W, std::int64_t> v_aik_reduced =
    needs_reduce.select(v_aik_coeff / v_gcd, v_aik_coeff);
  simd::Vec<W, std::int64_t> v_ajk_reduced =
    needs_reduce.select(v_ajk_broadcast / v_gcd, v_ajk_broadcast);

  // Broadcast coefficients to SIMD vectors for row operations
  std::array<simd::Vec<W, std::int64_t>, B> vAjk, vAik, vg;
#pragma clang loop unroll(full)
  for (std::ptrdiff_t b = 0; b < B; ++b) {
    vAjk[b] = simd::vbroadcast<W, std::int64_t>(v_ajk_reduced[b]);
    vAik[b] = simd::vbroadcast<W, std::int64_t>(v_aik_reduced[b]);
    vg[b] = simd::Vec<W, std::int64_t>{v_ajk_reduced[b]};
  }

  PtrMatrix<std::int64_t> C = A; // const ref
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol());
  std::ptrdiff_t global_l = 0;

  // SIMD mask for GCD tracking: active if Ajk_reduced != 1
  auto active_mask = simd::cmp::ne<W, std::int64_t>(v_ajk_reduced, one);

  // Main loop: process columns with GCD tracking
  if (active_mask.any()) {
    for (; global_l < L; global_l += W) {
      auto u{simd::index::unrollmask<1, W>(L, global_l)};
      if (!u) break;

      // Load pivot row once
      simd::Vec<W, std::int64_t> Ajl = C[j, u].vec_;

// Apply to all B target rows
#pragma clang loop unroll(full)
      for (std::ptrdiff_t b = 0; b < B; ++b) {
        std::ptrdiff_t row_i = rows[b];
        simd::Vec<W, std::int64_t> Ail =
          (vAjk[b] * C[row_i, u].vec_) - (vAik[b] * Ajl);
        A[row_i, u] = Ail;

        // GCD tracking for active lanes (check via intmask bit)
        if (active_mask.intmask() & (1ULL << b)) {
          vg[b] = gcd<W>(Ail, vg[b]);
          // Check if this lane's GCD is still > 1
          if (!simd::cmp::gt<W, std::int64_t>(vg[b], one).any()) {
            // Clear this lane from active_mask
            active_mask =
              decltype(active_mask){active_mask.intmask() & ~(1ULL << b)};
          }
        }
      }
      if (!active_mask.any()) {
        global_l += W; // Move past the column we just processed
        break;
      }
    }
  }

  // Continue without GCD tracking for remaining columns
  for (; global_l < L; global_l += W) {
    auto u{simd::index::unrollmask<1, W>(L, global_l)};
    if (!u) break;

    simd::Vec<W, std::int64_t> Ajl = C[j, u].vec_;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t b = 0; b < B; ++b) {
      std::ptrdiff_t row_i = rows[b];
      A[row_i, u] = (vAjk[b] * C[row_i, u].vec_) - (vAik[b] * Ajl);
    }
  }

  // Post-loop: apply GCD reduction to each row
  // Check active lanes via intmask
  std::uint64_t remaining_mask = active_mask.intmask();
#pragma clang loop unroll(full)
  for (std::ptrdiff_t b = 0; b < B; ++b) {
    if ((remaining_mask & (1ULL << b)) &&
        simd::cmp::gt<W, std::int64_t>(vg[b], one).any()) {
      std::int64_t g = gcdreduce<W>(vg[b]);
      if (g > 1) {
        std::ptrdiff_t row_i = rows[b];
        for (std::ptrdiff_t ll = 0; ll < L; ++ll)
          if (std::int64_t Ail = A[row_i, ll]) A[row_i, ll] = Ail / g;
      }
    }
  }
}

// Batched elimination with factor tracking
// zero_idx: index in batch where row 0 is (-1 if not present)
template <std::ptrdiff_t B>
auto zeroWithRowOp(MutPtrMatrix<std::int64_t> A,
                   std::array<std::ptrdiff_t, B> rows,
                   std::array<std::int64_t, B> coeffs, Row<> j, Col<> k,
                   std::int64_t f, std::ptrdiff_t zero_idx) -> std::int64_t {
  static constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
  static_assert(B == W, "Batch size must equal SIMD width");
  static constexpr simd::Vec<W, std::int64_t> one =
    simd::Vec<W, std::int64_t>{} + 1;

  std::int64_t Ajk = A[j, k];
  invariant(Ajk != 0);

  // Pre-compute GCD-reduced coefficients using SIMD
  // Load coeffs using unaligned SIMD load (B == W by design)
  simd::Vec<W, std::int64_t> v_aik_coeff =
    simd::load(coeffs.data(), simd::mask::None<W>{});
  simd::Vec<W, std::int64_t> v_ajk_broadcast =
    simd::vbroadcast<W, std::int64_t>(Ajk);

  // Parallel GCD computation for all B pairs
  simd::Vec<W, std::int64_t> v_gcd = gcd<W>(v_aik_coeff, v_ajk_broadcast);

  // Parallel division with mask for g != 1
  auto needs_reduce = simd::cmp::ne<W, std::int64_t>(v_gcd, one);
  simd::Vec<W, std::int64_t> v_aik_reduced =
    needs_reduce.select(v_aik_coeff / v_gcd, v_aik_coeff);
  simd::Vec<W, std::int64_t> v_ajk_reduced =
    needs_reduce.select(v_ajk_broadcast / v_gcd, v_ajk_broadcast);

  // Compute ret = f * Ajk_reduced[zero_idx] if zero_idx is valid
  std::int64_t ret = f;
  if (zero_idx >= 0) ret = f * v_ajk_reduced[zero_idx];

  // Broadcast coefficients to SIMD vectors for row operations
  std::array<simd::Vec<W, std::int64_t>, B> vAjk, vAik, vg;
#pragma clang loop unroll(full)
  for (std::ptrdiff_t b = 0; b < B; ++b) {
    vAjk[b] = simd::vbroadcast<W, std::int64_t>(v_ajk_reduced[b]);
    vAik[b] = simd::vbroadcast<W, std::int64_t>(v_aik_reduced[b]);
    vg[b] = (b == zero_idx) ? simd::Vec<W, std::int64_t>{ret}
                            : simd::Vec<W, std::int64_t>{v_ajk_reduced[b]};
  }

  PtrMatrix<std::int64_t> C = A; // const ref
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol());
  std::ptrdiff_t global_l = 0;

  // SIMD mask for GCD tracking
  // Build check vector: use ret for zero_idx, Ajk_reduced otherwise
  simd::Vec<W, std::int64_t> check_vec = v_ajk_reduced;
  if (zero_idx >= 0) check_vec[zero_idx] = ret;
  auto active_mask = simd::cmp::ne<W, std::int64_t>(check_vec, one);

  // Main loop: process columns with GCD tracking
  if (active_mask.any()) {
    for (; global_l < L; global_l += W) {
      auto u{simd::index::unrollmask<1, W>(L, global_l)};
      if (!u) break;

      // Load pivot row once
      simd::Vec<W, std::int64_t> Ajl = C[j, u].vec_;

      // Apply to all B target rows
#pragma clang loop unroll(full)
      for (std::ptrdiff_t b = 0; b < B; ++b) {
        std::ptrdiff_t row_i = rows[b];
        simd::Vec<W, std::int64_t> Ail =
          (vAjk[b] * C[row_i, u].vec_) - (vAik[b] * Ajl);
        A[row_i, u] = Ail;

        // GCD tracking for active lanes (check via intmask bit)
        if (active_mask.intmask() & (1ULL << b)) {
          vg[b] = gcd<W>(Ail, vg[b]);
          // Check if this lane's GCD is still > 1
          if (!simd::cmp::gt<W, std::int64_t>(vg[b], one).any()) {
            // Clear this lane from active_mask
            active_mask =
              decltype(active_mask){active_mask.intmask() & ~(1ULL << b)};
          }
        }
      }
      if (!active_mask.any()) {
        global_l += W; // Move past the column we just processed
        break;
      }
    }
  }

  // Continue without GCD tracking for remaining columns
  for (; global_l < L; global_l += W) {
    auto u{simd::index::unrollmask<1, W>(L, global_l)};
    if (!u) break;

    simd::Vec<W, std::int64_t> Ajl = C[j, u].vec_;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t b = 0; b < B; ++b) {
      std::ptrdiff_t row_i = rows[b];
      A[row_i, u] = (vAjk[b] * C[row_i, u].vec_) - (vAik[b] * Ajl);
    }
  }

  // Post-loop: apply GCD reduction to each row
  // Check active lanes via intmask
  std::uint64_t remaining_mask = active_mask.intmask();
  for (std::ptrdiff_t b = 0; b < B; ++b) {
    if ((remaining_mask & (1ULL << b)) &&
        simd::cmp::gt<W, std::int64_t>(vg[b], one).any()) {
      std::int64_t g = gcdreduce<W>(vg[b]);
      if (g > 1) {
        std::ptrdiff_t row_i = rows[b];
        for (std::ptrdiff_t ll = 0; ll < L; ++ll)
          if (std::int64_t Ail = A[row_i, ll]) A[row_i, ll] = Ail / g;
        if (b == zero_idx) {
          std::int64_t r = ret / g;
          invariant(r * g, ret);
          ret = r;
        }
      }
    }
  }
  return ret;
}

// Explicit template instantiations
template auto
collectBatch<DEFAULT_BATCH>(PtrMatrix<std::int64_t> A, std::ptrdiff_t num_rows,
                            std::ptrdiff_t pivot_row, Col<> col_k,
                            std::ptrdiff_t &start_row,
                            std::array<std::ptrdiff_t, DEFAULT_BATCH> &rows,
                            std::array<std::int64_t, DEFAULT_BATCH> &coeffs)
  -> std::pair<std::ptrdiff_t, std::ptrdiff_t>;

template void zeroWithRowOp<DEFAULT_BATCH>(
  MutPtrMatrix<std::int64_t> A, std::array<std::ptrdiff_t, DEFAULT_BATCH> rows,
  std::array<std::int64_t, DEFAULT_BATCH> coeffs, Row<> j, Col<> k);

template auto zeroWithRowOp<DEFAULT_BATCH>(
  MutPtrMatrix<std::int64_t> A, std::array<std::ptrdiff_t, DEFAULT_BATCH> rows,
  std::array<std::int64_t, DEFAULT_BATCH> coeffs, Row<> j, Col<> k,
  std::int64_t f, std::ptrdiff_t zero_idx) -> std::int64_t;

void zeroColumn(MutPtrMatrix<std::int64_t> A, Row<> pivot_row,
                Col<> pivot_col) {
  std::ptrdiff_t M = std::ptrdiff_t(A.numRow());
  std::ptrdiff_t pivot = std::ptrdiff_t(pivot_row);

  constexpr std::ptrdiff_t B = DEFAULT_BATCH;
  std::array<std::ptrdiff_t, B> rows;
  std::array<std::int64_t, B> coeffs;
  std::ptrdiff_t i = 0;

  while (i < M) {
    auto [count, zero_idx] =
      collectBatch<B>(A, M, pivot, pivot_col, i, rows, coeffs);
    if (count == 0) break;
    if (count == B) {
      zeroWithRowOp<B>(A, rows, coeffs, pivot_row, pivot_col);
    } else {
      for (std::ptrdiff_t j = 0; j < count; ++j)
        zeroWithRowOp(A, row(rows[j]), pivot_row, pivot_col);
      break;
    }
  }
}

auto zeroColumn(MutPtrMatrix<std::int64_t> A, Row<> pivot_row, Col<> pivot_col,
                std::int64_t f) -> std::int64_t {
  std::ptrdiff_t M = std::ptrdiff_t(A.numRow());
  std::ptrdiff_t pivot = std::ptrdiff_t(pivot_row);

  constexpr std::ptrdiff_t B = DEFAULT_BATCH;
  std::array<std::ptrdiff_t, B> rows;
  std::array<std::int64_t, B> coeffs;
  std::ptrdiff_t i = 0;

  while (i < M) {
    auto [count, zero_idx] =
      collectBatch<B>(A, M, pivot, pivot_col, i, rows, coeffs);
    if (count == 0) break;
    if (count == B) {
      if (zero_idx >= 0)
        f =
          zeroWithRowOp<B>(A, rows, coeffs, pivot_row, pivot_col, f, zero_idx);
      else zeroWithRowOp<B>(A, rows, coeffs, pivot_row, pivot_col);
    } else {
      for (std::ptrdiff_t j = 0; j < count; ++j)
        if (rows[j] == 0) f = zeroWithRowOp(A, row(0), pivot_row, pivot_col, f);
        else zeroWithRowOp(A, row(rows[j]), pivot_row, pivot_col);
      break;
    }
  }
  return f;
}

void bareiss(MutPtrMatrix<std::int64_t> A,
             MutPtrVector<std::ptrdiff_t> pivots) {
  const auto [M, N] = shape(A);
  invariant(pivots.size(), std::min(M, N));
  std::int64_t prev = 1, piv_ind = 0;
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c) {
    if (std::optional<std::ptrdiff_t> piv =
          ::pivotRowsBareiss(A, c, row(M), row(r))) {
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

[[nodiscard]] auto bareiss(IntMatrix<> &A) -> Vector<std::ptrdiff_t> {
  Vector<std::ptrdiff_t> pivots(length(A.minRowCol()));
  bareiss(A, pivots);
  return pivots;
}

void solveColumn(MutPtrMatrix<std::int64_t> A, MutPtrMatrix<std::int64_t> B,
                 std::ptrdiff_t r, std::ptrdiff_t c) {
  const auto [M, N] = shape(A);
  utils::assume(B.numRow() == M);
  utils::invariant(r < M);
  utils::invariant(c < N);
  utils::invariant(c < B.numCol());
  utils::invariant(!::pivotRowsPair({A, B}, col(c), row(M), row(r)));
  ::zeroColumnPair({A, B}, col(c), row(r));
}
/// void solveSystem(IntMatrix &A, IntMatrix &B)
/// Say we wanted to solve \f$\textbf{AX} = \textbf{B}\f$.
/// `solveSystem` left-multiplies both sides by
/// a matrix \f$\textbf{W}\f$ that diagonalizes \f$\textbf{A}\f$.
/// Once \f$\textbf{A}\f$ has been diagonalized, the solution is trivial.
/// Both inputs are overwritten with the product of the left multiplications.
void solveSystem(MutPtrMatrix<std::int64_t> A, MutPtrMatrix<std::int64_t> B) {
  const auto [M, N] = shape(A);
  utils::assume(B.numRow() == M);
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++c)
    if (!::pivotRowsPair({A, B}, col(c), row(M), row(r)))
      ::zeroColumnPair({A, B}, col(c), row(r++));
}
// like solveSystem, except it right-multiplies.
// That is, given `XA = B`, it right-multiplies both sides by
// a matrix to diagonalize `A`.
void solveSystemRight(MutPtrMatrix<std::int64_t> A,
                      MutPtrMatrix<std::int64_t> B) {
  const auto [M, N] = shape(A);
  utils::assume(B.numCol() == N);
  for (std::ptrdiff_t r = 0, c = 0; c < N && r < M; ++r)
    if (!::pivotColsPair({A, B}, row(r), col(N), col(c)))
      ::zeroColumnPair({A, B}, row(r), col(c++));
}
// diagonalizes A(0:K,0:K)
void solveSystem(MutPtrMatrix<std::int64_t> A, std::ptrdiff_t K) {
  Row M = A.numRow();
  for (std::ptrdiff_t r = 0, c = 0; c < K && r < M; ++c)
    if (!pivotRows(A, col(c), M, row(r))) ::zeroColumn(A, col(c), row(r++));
}
// diagonalizes A(0:K,1:K+1)
void solveSystemSkip(MutPtrMatrix<std::int64_t> A) {
  const auto [M, N] = shape(A);
  for (std::ptrdiff_t r = 0, c = 1; c < N && r < M; ++c)
    if (!pivotRows(A, col(c), row(M), row(r)))
      ::zeroColumn(A, col(c), row(r++));
}

// returns `true` if the solve failed, `false` otherwise
// diagonals contain denominators.
// Assumes the last column is the vector to solve for.
void solveSystem(MutPtrMatrix<std::int64_t> A) {
  solveSystem(A, std::ptrdiff_t(A.numCol()) - 1);
}

/// inv(A) -> (D, B)
/// Given a matrix \f$\textbf{A}\f$, returns two matrices \f$\textbf{D}\f$ and
/// \f$\textbf{B}\f$ so that \f$\textbf{D}^{-1}\textbf{B} = \textbf{A}^{-1}\f$,
/// and \f$\textbf{D}\f$ is diagonal.
/// NOTE: This function assumes non-singular
/// Mutates `A`
// NOLINTNEXTLINE(performance-unnecessary-value-param)
[[nodiscard]] auto inv(Arena<> *alloc, MutSquarePtrMatrix<std::int64_t> A)
  -> MutSquarePtrMatrix<std::int64_t> {
  MutSquarePtrMatrix<std::int64_t> B =
    identity<std::int64_t>(alloc, std::ptrdiff_t(A.numCol()));
  solveSystem(A, B);
  return B;
}
// scaledInv(A, B) -> s
// reads and writes A, writes B
// B = s * inv(A) // A is diagonalized in the process
[[nodiscard]] auto scaledInv(MutSquarePtrMatrix<std::int64_t> A,
                             MutSquarePtrMatrix<std::int64_t> B)
  -> std::int64_t {
  B.zero();
  B.diag() << 1;
  solveSystem(A, B);
  auto [s, nonUnity] = lcmNonUnity(A.diag());
  if (nonUnity) B *= s / A.diag();
  // for (std::ptrdiff_t i = 0; i < A.numRow(); ++i) B[i, _] *= s / A[i, i];
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
[[nodiscard]] auto scaledInv(Arena<> *alloc, MutSquarePtrMatrix<std::int64_t> A)
  -> containers::Pair<MutSquarePtrMatrix<std::int64_t>, std::int64_t> {
  MutSquarePtrMatrix<std::int64_t> B =
    square_matrix<std::int64_t>(alloc, std::ptrdiff_t(A.numCol()));
  return {B, scaledInv(A, B)};
}

auto nullSpace(MutDensePtrMatrix<std::int64_t> B,
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
auto nullSpace(Arena<> *alloc, MutDensePtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  Row M = A.numRow();
  return nullSpace(matrix<std::int64_t>(alloc, M, math::col(std::ptrdiff_t(M))),
                   A);
}
auto nullSpace(Arena<> *alloc, MutPtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  Row M = A.numRow();
  return ::nullSpace(
    matrix<std::int64_t>(alloc, M, math::col(std::ptrdiff_t(M))), A);
}

// FIXME: why do we have two?
auto orthogonalize(IntMatrix<> A)
  -> containers::Pair<SquareMatrix<std::int64_t>, Vector<unsigned>> {
  return ::orthogonalizeBang(A);
}

} // namespace NormalForm

// FIXME: why do we have two?
auto orthogonalize(Arena<> *alloc, MutDensePtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  std::ptrdiff_t num_rows = std::ptrdiff_t(A.numRow()),
                 num_cols = std::ptrdiff_t(A.numCol());
  if ((num_cols < 2) || (num_rows == 0)) return A;
  normalizeByGCD(A[0, _]);
  if (num_rows == 1) return A;
  auto s = alloc->scope();
  MutPtrVector<Rational> buff{vector<Rational>(alloc, num_cols)};
  std::ptrdiff_t offset = 0;
  while (allZero(A[0, _])) A[0, _] << A[--num_rows, _];
  for (std::ptrdiff_t i = 1; i < num_rows; ++i) {
    buff << A[i, _];
    for (std::ptrdiff_t j = 0; j < i - offset; ++j) {
      std::int64_t n = 0;
      std::int64_t d = 0;
      for (std::ptrdiff_t k = 0; k < num_cols; ++k) {
        n += A[i, k] * A[j, k];
        d += A[j, k] * A[j, k];
      }
      for (std::ptrdiff_t k = 0; k < num_cols; ++k)
        buff[k] -= Rational::createPositiveDenominator(A[j, k] * n, d);
    }
    bool all_zero = true;
    for (std::ptrdiff_t k = 0; k < num_cols; ++k) {
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
    for (std::ptrdiff_t k = 0; k < num_cols; ++k)
      lm = lcm(lm, buff[k].denominator_);
    for (std::ptrdiff_t k = 0; k < num_cols; ++k)
      A[i - offset, k] = buff[k].numerator_ * (lm / buff[k].denominator_);
  }
  return A[_(num_rows - offset), _];
}

[[nodiscard]] auto orthogonalNullSpace(Arena<> *alloc,
                                       MutDensePtrMatrix<std::int64_t> A)
  -> MutDensePtrMatrix<std::int64_t> {
  return orthogonalize(alloc, NormalForm::nullSpace(alloc, A));
}

} // namespace math

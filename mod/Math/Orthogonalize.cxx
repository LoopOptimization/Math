module;

#include <cstddef>
#include <cstdint>

export module Orthogonalize;

import Array;
import ManagedArray;
import NormalForm;
import Pair;
import Rational;
import VGCD;

namespace math {
constexpr auto orthogonalizeBang(MutDensePtrMatrix<int64_t> &A)
  -> containers::Pair<SquareMatrix<int64_t>, Vector<unsigned>> {
  // we try to orthogonalize with respect to as many rows of `A` as we can
  // prioritizing earlier rows.
  auto [M, N] = shape(A);
  SquareMatrix<int64_t> K{identity(math::DefaultAlloc<int64_t>{}, unsigned(M))};
  Vector<unsigned> included;
  included.reserve(std::min(ptrdiff_t(M), ptrdiff_t(N)));
  for (ptrdiff_t i = 0, j = 0; i < std::min(ptrdiff_t(M), ptrdiff_t(N)); ++j) {
    // zero ith row
    if (pivotRows(A, K, i, row(M))) {
      // cannot pivot, this is a linear combination of previous
      // therefore, we drop the row
      dropCol(A, i, row(M), col(--N));
    } else {
      zeroSupDiagonal(A, K, i, row(M), col(N));
      int64_t Aii = A[i, i];
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
} // namespace math
// end unexported
// begin exported
export namespace math {
constexpr auto orthogonalize(IntMatrix<> A)
  -> containers::Pair<SquareMatrix<int64_t>, Vector<unsigned>> {
  return orthogonalizeBang(A);
}

[[nodiscard]] constexpr auto
orthogonalize(DenseMatrix<int64_t> A) -> DenseMatrix<int64_t> {
  if ((A.numCol() < 2) || (A.numRow() == 0)) return A;
  normalizeByGCD(A[0, _]);
  if (A.numRow() == 1) return A;
  Vector<Rational, 8> buff;
  buff.resizeForOverwrite(ptrdiff_t(A.numCol()));
  for (ptrdiff_t i = 1; i < A.numRow(); ++i) {
    for (ptrdiff_t j = 0; j < A.numCol(); ++j) buff[j] = A[i, j];
    for (ptrdiff_t j = 0; j < i; ++j) {
      int64_t n = 0;
      int64_t d = 0;
      for (ptrdiff_t k = 0; k < A.numCol(); ++k) {
        n += A[i, k] * A[j, k];
        d += A[j, k] * A[j, k];
      }
      for (ptrdiff_t k = 0; k < A.numCol(); ++k)
        buff[k] -= Rational::createPositiveDenominator(A[j, k] * n, d);
    }
    int64_t lm = 1;
    for (ptrdiff_t k = 0; k < A.numCol(); ++k)
      lm = lcm(lm, buff[k].denominator);
    for (ptrdiff_t k = 0; k < A.numCol(); ++k)
      A[i, k] = buff[k].numerator * (lm / buff[k].denominator);
  }
  return A;
}

[[nodiscard]] constexpr auto
orthogonalNullSpace(DenseMatrix<int64_t> A) -> DenseMatrix<int64_t> {
  return orthogonalize(NormalForm::nullSpace(std::move(A)));
}
} // namespace math

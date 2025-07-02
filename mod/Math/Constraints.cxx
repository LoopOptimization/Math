module;
#include "Macros.hxx"
module Constraints;
import ArrayConcepts;
import AxisTypes;
import BitSet;
import Comparisons;
import CorePrint;
import EmptyMatrix;
import GCD;
import GenericArrayConstructors;
import ManagedArray;
import NormalForm;
import std;

namespace math {

auto substituteEqualityImpl(MutDensePtrMatrix<std::int64_t> E,
                            const std::ptrdiff_t i) -> Row<> {
  const auto [numConstraints, numVar] = shape(E);
  std::ptrdiff_t minNonZero = numVar + 1;
  std::ptrdiff_t rowMinNonZero = numConstraints;
  for (std::ptrdiff_t j = 0; j < numConstraints; ++j)
    if (E[j, i]) {
      std::ptrdiff_t nonZero = 0;
      for (std::ptrdiff_t v = 0; v < numVar; ++v) nonZero += (E[j, v] != 0);
      if (nonZero < minNonZero) {
        minNonZero = nonZero;
        rowMinNonZero = j;
      }
    }
  if (rowMinNonZero == numConstraints) return row(rowMinNonZero);
  auto Es = E[rowMinNonZero, _];
  std::int64_t Eis = Es[i];
  // we now substitute the equality expression with the minimum number
  // of terms.
  if (constexpr_abs(Eis) == 1) {
    for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
      if (j == rowMinNonZero) continue;
      if (std::int64_t Eij = E[j, i]) E[j, _] << Eis * E[j, _] - Eij * Es;
    }
  } else {
    for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
      if (j == rowMinNonZero) continue;
      if (std::int64_t Eij = E[j, i]) {
        std::int64_t g = gcd(Eij, Eis);
        E[j, _] << (Eis / g) * E[j, _] - (Eij / g) * Es;
      }
    }
  }
  return row(rowMinNonZero);
}

auto substituteEquality(DenseMatrix<std::int64_t> &E, const std::ptrdiff_t i)
  -> bool {
  Row minNonZero = substituteEqualityImpl(E, i);
  if (minNonZero == E.numRow()) return true;
  eraseConstraint(E, minNonZero);
  return false;
}

auto substituteEqualityPairImpl(
  std::array<MutDensePtrMatrix<std::int64_t>, 2> AE, std::ptrdiff_t i)
  -> Row<> {
  auto [A, E] = AE;
  const auto [numConstraints, numVar] = shape(E);
  std::ptrdiff_t minNonZero = numVar + 1;
  std::ptrdiff_t rowMinNonZero = numConstraints;
  for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
    if (E[j, i]) {
      std::ptrdiff_t nonZero = 0;
      for (std::ptrdiff_t v = 0; v < numVar; ++v) nonZero += (E[j, v] != 0);
      if (nonZero < minNonZero) {
        minNonZero = nonZero;
        rowMinNonZero = j;
      }
    }
  }
  if (rowMinNonZero == numConstraints) return row(rowMinNonZero);
  auto Es = E[rowMinNonZero, _];
  std::int64_t Eis = Es[i], s = 2 * (Eis > 0) - 1;
  // we now substitute the equality expression with the minimum number
  // of terms.
  if (constexpr_abs(Eis) == 1) {
    for (std::ptrdiff_t j = 0; j < A.numRow(); ++j)
      if (std::int64_t Aij = A[j, i])
        A[j, _] << (s * Eis) * A[j, _] - (s * Aij) * Es;
    for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
      if (j == rowMinNonZero) continue;
      if (std::int64_t Eij = E[j, i]) E[j, _] << Eis * E[j, _] - Eij * Es;
    }
  } else {
    for (std::ptrdiff_t j = 0; j < A.numRow(); ++j) {
      if (std::int64_t Aij = A[j, i]) {
        std::int64_t g = gcd(Aij, Eis);
        invariant(g > 0);
        // `A` contains inequalities; flipping signs is illegal
        A[j, _] << ((s * Eis) / g) * A[j, _] - ((s * Aij) / g) * Es;
      }
    }
    for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
      if (j == rowMinNonZero) continue;
      if (std::int64_t Eij = E[j, i]) {
        std::int64_t g = gcd(Eij, Eis);
        E[j, _] << (Eis / g) * E[j, _] - (Eij / g) * Es;
      }
    }
  }
  return row(rowMinNonZero);
}

auto substituteEquality(MutDensePtrMatrix<std::int64_t> &A,
                        MutDensePtrMatrix<std::int64_t> &E,
                        const std::ptrdiff_t i) -> bool {

  Row minNonZero = substituteEqualityPairImpl({A, E}, i);
  if (minNonZero == E.numRow()) return true;
  eraseConstraint(E, minNonZero);
  return false;
}

void fourierMotzkinCore(DenseMatrix<std::int64_t> &A, std::ptrdiff_t v,
                        std::array<std::ptrdiff_t, 2> negPos) {
  auto [numNeg, numPos] = negPos;
  // we need one extra, as on the last overwrite, we still need to
  // read from two constraints we're deleting; we can't write into
  // both of them. Thus, we use a little extra memory here,
  // and then truncate.
  const Row numRowsOld = A.numRow(),
            numRowsNew = row(std::ptrdiff_t(numRowsOld) - numNeg - numPos +
                             numNeg * numPos + 1);
  A.resize(numRowsNew);
  // plan is to replace
  for (std::ptrdiff_t i = 0, numRows = std::ptrdiff_t(numRowsOld),
                      posCount = numPos;
       posCount; ++i) {
    std::int64_t Aiv = A[i, v];
    if (Aiv <= 0) continue;
    --posCount;
    for (std::ptrdiff_t negCount = numNeg, j = 0; negCount; ++j) {
      std::int64_t Ajv = A[j, v];
      if (Ajv >= 0) continue;
      // for the last `negCount`, we overwrite `A(i, k)`
      // last posCount does not get overwritten
      --negCount;
      std::ptrdiff_t c =
        posCount ? (negCount ? numRows++ : std::ptrdiff_t(i)) : j;
      std::int64_t Ai = Aiv, Aj = Ajv;
      std::int64_t g = gcd(Aiv, Ajv);
      if (g != 1) {
        Ai /= g;
        Aj /= g;
      }
      bool allZero = true;
      for (std::ptrdiff_t k = 0; k < A.numCol(); ++k) {
        std::int64_t Ack = Ai * A[j, k] - Aj * A[i, k];
        A[c, k] = Ack;
        allZero &= (Ack == 0);
      }
      if (allZero) {
        eraseConstraint(A, row(c));
        if (posCount)
          if (negCount) --numRows;
          else --i;
        else --j;
      }
    }
    if (posCount == 0) // last posCount not overwritten, so we erase
      eraseConstraint(A, row(i));
  }
}

void fourierMotzkin(DenseMatrix<std::int64_t> &A, std::ptrdiff_t v) {
  invariant(v < A.numCol());
  const auto [numNeg, numPos] = countNonZeroSign(A, v);
  const Row numRowsOld = A.numRow();
  if ((numNeg == 0) | (numPos == 0)) {
    if ((numNeg == 0) & (numPos == 0)) return;
    for (Row i = numRowsOld; i != 0;)
      if (A[--i, v]) eraseConstraint(A, i);
    return;
  }
  fourierMotzkinCore(A, v, {numNeg, numPos});
}

auto removeRedundantRows(MutDensePtrMatrix<std::int64_t> A,
                         MutDensePtrMatrix<std::int64_t> B)
  -> std::array<Row<>, 2> {
  auto [M, N] = shape(B);
  for (std::ptrdiff_t r = 0, c = 0; c++ < N && r < M;)
    if (!NormalForm::pivotRows(B, col(c == N ? 0 : c), row(M), row(r)))
      NormalForm::reduceColumnStack(A, B, c == N ? 0 : c, r++);
  // scan duplicate rows in `A`
  for (Row r = A.numRow(); r != 0;)
    if (!uniqueConstraint(A, --r)) eraseConstraint(A, r);
  return {NormalForm::numNonZeroRows(A), NormalForm::numNonZeroRows(B)};
}

void printConstraints(DensePtrMatrix<std::int64_t> A, bool inequality) {
  Row numConstraints = A.numRow();
  for (std::ptrdiff_t c = 0; c < numConstraints; ++c) {
    printConstraint(A[c, _], 1, inequality);
    utils::print('\n');
  }
}

void eraseConstraintImpl(MutDensePtrMatrix<std::int64_t> A, std::ptrdiff_t _i,
                         std::ptrdiff_t _j) {
  invariant(_i != _j);
  std::ptrdiff_t i = std::min(_i, _j), j = std::max(_i, _j);
  const auto [M, N] = shape(A);
  std::ptrdiff_t lastRow = M - 1;
  std::ptrdiff_t penuRow = lastRow - 1;
  if (j == penuRow) {
    // then we only need to copy one column (i to lastCol)
    eraseConstraintImpl(A, row(i));
  } else if ((i != penuRow) && (i != lastRow)) {
    // if i == penuCol, then j == lastCol
    // and we thus don't need to copy
    for (std::ptrdiff_t n = 0; n < N; ++n) {
      A[i, n] = A[penuRow, n];
      A[j, n] = A[lastRow, n];
    }
  }
}

void slackEqualityConstraints(MutPtrMatrix<std::int64_t> C,
                              PtrMatrix<std::int64_t> A,
                              PtrMatrix<std::int64_t> B) {
  const Col numVar = A.numCol();
  invariant(numVar, B.numCol());
  const Row numSlack = A.numRow(), numStrict = B.numRow();
  invariant(C.numRow(), numSlack + numStrict);
  std::ptrdiff_t slackAndVar =
    std::ptrdiff_t(numSlack) + std::ptrdiff_t(numVar);
  invariant(std::ptrdiff_t(C.numCol()), slackAndVar);
  // [I A]
  for (std::ptrdiff_t s = 0; s < numSlack; ++s) {
    C[s, _(begin, numSlack)] << 0;
    C[s, s] = 1;
    C[s, _(numSlack, slackAndVar)] << A[s, _(begin, numVar)];
  }
  // [0 B]
  for (std::ptrdiff_t s = 0; s < numStrict; ++s) {
    C[s + std::ptrdiff_t(numSlack), _(begin, numSlack)] << 0;
    C[s + std::ptrdiff_t(numSlack), _(numSlack, slackAndVar)]
      << B[s, _(begin, numVar)];
  }
}

void slackEqualityConstraints(MutPtrMatrix<std::int64_t> C,
                              PtrMatrix<std::int64_t> A) {
  const Col numVar = A.numCol();
  const Row numSlack = A.numRow();
  std::ptrdiff_t slackAndVar =
    std::ptrdiff_t(numSlack) + std::ptrdiff_t(numVar);
  invariant(std::ptrdiff_t(C.numCol()), slackAndVar);
  // [I A]
  for (std::ptrdiff_t s = 0; s < numSlack; ++s) {
    C[s, _(begin, numSlack)] << 0;
    C[s, s] = 1;
    C[s, _(numSlack, slackAndVar)] << A[s, _(begin, numVar)];
  }
}

void removeZeroRows(MutDensePtrMatrix<std::int64_t> &A) {
  for (Row i = A.numRow(); i;)
    if (allZero(A[--i, _])) eraseConstraint(A, i);
}

template <bool NonNegative>
auto fourierMotzkinCore(MutDensePtrMatrix<std::int64_t> B,
                        DensePtrMatrix<std::int64_t> A, std::ptrdiff_t v,
                        std::array<containers::BitSet64, 3> znp) -> Row<> {
  const auto &[zero, neg, pos] = znp;
  // we have the additional v >= 0
  if constexpr (NonNegative)
    invariant(B.numRow() == std::ptrdiff_t(A.numRow()) - pos.size() +
                              std::ptrdiff_t(neg.size()) * pos.size());
  else
    invariant(B.numRow() == std::ptrdiff_t(A.numRow()) - pos.size() -
                              neg.size() +
                              std::ptrdiff_t(neg.size()) * pos.size());
  invariant(++auto{B.numCol()}, A.numCol());
  std::ptrdiff_t r = 0;
  // x - v >= 0 -> x >= v
  // x + v >= 0 -> v >= -x
  for (auto i : neg) {
    // we  have implicit v >= 0, matching x >= v
    if constexpr (NonNegative) {
      B[r, _(0, v)] << A[i, _(0, v)];
      B[r, _(v, end)] << A[i, _(v + 1, end)];
      r += anyNEZero(B[r, _(0, end)]);
    }
    std::int64_t Aiv = A[i, v];
    invariant(Aiv < 0);
    for (auto j : pos) {
      std::int64_t Ajv = A[j, v];
      invariant(Ajv > 0);
      auto [ai, aj] = divgcd(Aiv, Ajv);
      B[r, _(0, v)] << aj * A[i, _(0, v)] - ai * A[j, _(0, v)];
      B[r, _(v, end)] << aj * A[i, _(v + 1, end)] - ai * A[j, _(v + 1, end)];
      r += anyNEZero(B[r, _(0, end)]);
    }
  }
  for (auto i : zero) {
    B[r, _(0, v)] << A[i, _(0, v)];
    B[r, _(v, end)] << A[i, _(v + 1, end)];
    r += anyNEZero(B[r, _(0, end)]);
  }
  return row(r);
}

// template <bool nonnegative>
// auto fouriermotzkin(alloc<std::int64_t> auto &alloc,
//                     denseptrmatrix<std::int64_t> a, std::ptrdiff_t v)
//   -> mutdenseptrmatrix<std::int64_t> {

//   auto znp = indszeronegpos(a[_, v]);
//   auto &[zero, neg, pos] = znp;
//   std::ptrdiff_t r = std::ptrdiff_t(a.numrow()) - pos.size() +
//                      std::ptrdiff_t(neg.size()) * pos.size();
//   if constexpr (!nonnegative) r -= neg.size();
//   auto b = matrix(alloc, row(r), --auto{a.numcol()});
//   b.truncate(fouriermotzkincore<nonnegative>(b, a, v, znp));
//   return b;
// }

// Explicit instantiation definitions
template auto fourierMotzkinCore<true>(MutDensePtrMatrix<std::int64_t> B,
                                       DensePtrMatrix<std::int64_t> A,
                                       std::ptrdiff_t v,
                                       std::array<containers::BitSet64, 3> znp)
  -> Row<>;

template auto fourierMotzkinCore<false>(MutDensePtrMatrix<std::int64_t> B,
                                        DensePtrMatrix<std::int64_t> A,
                                        std::ptrdiff_t v,
                                        std::array<containers::BitSet64, 3> znp)
  -> Row<>;

// template auto fourierMotzkin<true>(Alloc<std::int64_t> auto &alloc,
//                                    DensePtrMatrix<std::int64_t> A,
//                                    std::ptrdiff_t v)
//   -> MutDensePtrMatrix<std::int64_t>;

// template auto fourierMotzkin<false>(Alloc<std::int64_t> auto &alloc,
//                                     DensePtrMatrix<std::int64_t> A,
//                                     std::ptrdiff_t v)
//   -> MutDensePtrMatrix<std::int64_t>;

} // namespace math

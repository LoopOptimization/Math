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

auto substituteEqualityImpl(MutDensePtrMatrix<std::int64_t> E, std::ptrdiff_t i)
  -> Row<> {
  const auto [numConstraints, numVar] = shape(E);
  std::ptrdiff_t min_non_zero = numVar + 1;
  std::ptrdiff_t row_min_non_zero = numConstraints;
  for (std::ptrdiff_t j = 0; j < numConstraints; ++j)
    if (E[j, i]) {
      std::ptrdiff_t non_zero = 0;
      for (std::ptrdiff_t v = 0; v < numVar; ++v) non_zero += (E[j, v] != 0);
      if (non_zero < min_non_zero) {
        min_non_zero = non_zero;
        row_min_non_zero = j;
      }
    }
  if (row_min_non_zero == numConstraints) return row(row_min_non_zero);
  auto Es = E[row_min_non_zero, _];
  std::int64_t Eis = Es[i];
  // we now substitute the equality expression with the minimum number
  // of terms.
  if (constexpr_abs(Eis) == 1) {
    for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
      if (j == row_min_non_zero) continue;
      if (std::int64_t Eij = E[j, i]) E[j, _] << Eis * E[j, _] - Eij * Es;
    }
  } else {
    for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
      if (j == row_min_non_zero) continue;
      if (std::int64_t Eij = E[j, i]) {
        std::int64_t g = gcd(Eij, Eis);
        E[j, _] << (Eis / g) * E[j, _] - (Eij / g) * Es;
      }
    }
  }
  return row(row_min_non_zero);
}

auto substituteEquality(DenseMatrix<std::int64_t> &E, std::ptrdiff_t i)
  -> bool {
  Row min_non_zero = substituteEqualityImpl(E, i);
  if (min_non_zero == E.numRow()) return true;
  eraseConstraint(E, min_non_zero);
  return false;
}

auto substituteEqualityPairImpl(
  std::array<MutDensePtrMatrix<std::int64_t>, 2> AE, std::ptrdiff_t i)
  -> Row<> {
  auto [A, E] = AE;
  const auto [numConstraints, numVar] = shape(E);
  std::ptrdiff_t min_non_zero = numVar + 1;
  std::ptrdiff_t row_min_non_zero = numConstraints;
  for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
    if (E[j, i]) {
      std::ptrdiff_t non_zero = 0;
      for (std::ptrdiff_t v = 0; v < numVar; ++v) non_zero += (E[j, v] != 0);
      if (non_zero < min_non_zero) {
        min_non_zero = non_zero;
        row_min_non_zero = j;
      }
    }
  }
  if (row_min_non_zero == numConstraints) return row(row_min_non_zero);
  auto Es = E[row_min_non_zero, _];
  std::int64_t Eis = Es[i], s = (2 * (Eis > 0)) - 1;
  // we now substitute the equality expression with the minimum number
  // of terms.
  if (constexpr_abs(Eis) == 1) {
    for (std::ptrdiff_t j = 0; j < A.numRow(); ++j)
      if (std::int64_t Aij = A[j, i])
        A[j, _] << (s * Eis) * A[j, _] - (s * Aij) * Es;
    for (std::ptrdiff_t j = 0; j < numConstraints; ++j) {
      if (j == row_min_non_zero) continue;
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
      if (j == row_min_non_zero) continue;
      if (std::int64_t Eij = E[j, i]) {
        std::int64_t g = gcd(Eij, Eis);
        E[j, _] << (Eis / g) * E[j, _] - (Eij / g) * Es;
      }
    }
  }
  return row(row_min_non_zero);
}

auto substituteEquality(MutDensePtrMatrix<std::int64_t> &A,
                        MutDensePtrMatrix<std::int64_t> &E,
                        const std::ptrdiff_t i) -> bool {

  Row min_non_zero = substituteEqualityPairImpl({A, E}, i);
  if (min_non_zero == E.numRow()) return true;
  eraseConstraint(E, min_non_zero);
  return false;
}

auto fourierMotzkinCore(MutDensePtrMatrix<std::int64_t> A, std::ptrdiff_t v,
                        std::array<std::ptrdiff_t, 2> negPos,
                        std::ptrdiff_t num_rows) -> Row<> {
  auto [numNeg, numPos] = negPos;
  // we need one extra, as on the last overwrite, we still need to
  // read from two constraints we're deleting; we can't write into
  // both of them. Thus, we use a little extra memory here,
  // and then truncate.
  // plan is to replace
  for (std::ptrdiff_t i = 0, pos_count = numPos; pos_count; ++i) {
    std::int64_t Aiv = A[i, v];
    if (Aiv <= 0) continue;
    --pos_count;
    for (std::ptrdiff_t neg_count = numNeg, j = 0; neg_count; ++j) {
      std::int64_t Ajv = A[j, v];
      if (Ajv >= 0) continue;
      // for the last `negCount`, we overwrite `A(i, k)`
      // last posCount does not get overwritten
      --neg_count;
      std::ptrdiff_t c = pos_count ? (neg_count ? num_rows++ : i) : j;
      std::int64_t Ai = Aiv, Aj = Ajv, g = gcd(Aiv, Ajv);
      if (g != 1) Ai /= g, Aj /= g;
      bool all_zero = true;
      for (std::ptrdiff_t k = 0; k < A.numCol(); ++k) {
        std::int64_t Ack = (Ai * A[j, k]) - (Aj * A[i, k]);
        A[c, k] = Ack;
        all_zero &= (Ack == 0);
      }
      if (all_zero) {
        eraseConstraint(A, row(c));
        if (pos_count)
          if (neg_count) --num_rows;
          else --i;
        else --j;
      }
    }
    if (pos_count == 0) // last posCount not overwritten, so we erase
      eraseConstraint(A, row(i));
  }
  return A.numRow();
}
void fourierMotzkinCore(DenseMatrix<std::int64_t> &A, std::ptrdiff_t v,
                        std::array<std::ptrdiff_t, 2> negPos) {
  auto [numNeg, numPos] = negPos;
  // we need one extra, as on the last overwrite, we still need to
  // read from two constraints we're deleting; we can't write into
  // both of them. Thus, we use a little extra memory here,
  // and then truncate.
  Row num_rows_old = A.numRow(),
      num_rows_new = row(std::ptrdiff_t(num_rows_old) - numNeg - numPos +
                         (numNeg * numPos) + 1);
  A.resize(num_rows_new);
  A.resize(fourierMotzkinCore(A, v, negPos, std::ptrdiff_t(num_rows_old)));
}

void fourierMotzkin(DenseMatrix<std::int64_t> &A, std::ptrdiff_t v) {
  invariant(v < A.numCol());
  const auto [numNeg, numPos] = countNonZeroSign(A, v);
  const Row num_rows_old = A.numRow();
  if ((numNeg == 0) | (numPos == 0)) {
    if ((numNeg == 0) & (numPos == 0)) return;
    for (Row i = num_rows_old; i != 0;)
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
  Row num_constraints = A.numRow();
  for (std::ptrdiff_t c = 0; c < num_constraints; ++c) {
    printConstraint(A[c, _], 1, inequality);
    utils::print('\n');
  }
}

void eraseConstraintImpl(MutDensePtrMatrix<std::int64_t> A, std::ptrdiff_t _i,
                         std::ptrdiff_t _j) {
  invariant(_i != _j);
  std::ptrdiff_t i = std::min(_i, _j), j = std::max(_i, _j);
  const auto [M, N] = shape(A);
  std::ptrdiff_t last_row = M - 1, penu_row = last_row - 1;
  if (j == penu_row) {
    // then we only need to copy one column (i to lastCol)
    eraseConstraintImpl(A, row(i));
  } else if ((i != penu_row) && (i != last_row)) {
    // if i == penuCol, then j == lastCol
    // and we thus don't need to copy
    for (std::ptrdiff_t n = 0; n < N; ++n) {
      A[i, n] = A[penu_row, n];
      A[j, n] = A[last_row, n];
    }
  }
}

void slackEqualityConstraints(MutPtrMatrix<std::int64_t> C,
                              PtrMatrix<std::int64_t> A,
                              PtrMatrix<std::int64_t> B) {
  const Col num_var = A.numCol();
  invariant(num_var, B.numCol());
  const Row num_slack = A.numRow(), num_strict = B.numRow();
  invariant(C.numRow(), num_slack + num_strict);
  std::ptrdiff_t slack_and_var =
    std::ptrdiff_t(num_slack) + std::ptrdiff_t(num_var);
  invariant(std::ptrdiff_t(C.numCol()), slack_and_var);
  // [I A]
  for (std::ptrdiff_t s = 0; s < num_slack; ++s) {
    C[s, _(begin, num_slack)] << 0;
    C[s, s] = 1;
    C[s, _(num_slack, slack_and_var)] << A[s, _(begin, num_var)];
  }
  // [0 B]
  for (std::ptrdiff_t s = 0; s < num_strict; ++s) {
    C[s + std::ptrdiff_t(num_slack), _(begin, num_slack)] << 0;
    C[s + std::ptrdiff_t(num_slack), _(num_slack, slack_and_var)]
      << B[s, _(begin, num_var)];
  }
}

void slackEqualityConstraints(MutPtrMatrix<std::int64_t> C,
                              PtrMatrix<std::int64_t> A) {
  const Col num_var = A.numCol();
  const Row num_slack = A.numRow();
  std::ptrdiff_t slack_and_var =
    std::ptrdiff_t(num_slack) + std::ptrdiff_t(num_var);
  invariant(std::ptrdiff_t(C.numCol()), slack_and_var);
  // [I A]
  for (std::ptrdiff_t s = 0; s < num_slack; ++s) {
    C[s, _(begin, num_slack)] << 0;
    C[s, s] = 1;
    C[s, _(num_slack, slack_and_var)] << A[s, _(begin, num_var)];
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
                              (std::ptrdiff_t(neg.size()) * pos.size()));
  else
    invariant(B.numRow() == std::ptrdiff_t(A.numRow()) - pos.size() -
                              neg.size() +
                              (std::ptrdiff_t(neg.size()) * pos.size()));
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

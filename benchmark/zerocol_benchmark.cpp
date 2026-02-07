import Nanobench;
import Arena;
import ManagedArray;
import MatrixBinaryIO;
import NormalForm;
import GCD;
import SIMD;
import std;

namespace {

using namespace math;

// Find a (pivot_row, pivot_col) where the column has many nonzero entries
auto findPivot(PtrMatrix<std::int64_t> A) -> std::array<std::ptrdiff_t, 2> {
  auto [M, N] = shape(A);
  std::ptrdiff_t best_row = 0, best_col = 0, best_count = 0;
  for (std::ptrdiff_t c = 0; c < N; ++c) {
    std::ptrdiff_t count = 0;
    std::ptrdiff_t pivot_row = -1;
    for (std::ptrdiff_t r = 0; r < M; ++r) {
      if (A[r, c] != 0) {
        ++count;
        if (pivot_row < 0) pivot_row = r;
      }
    }
    if (count > best_count && pivot_row >= 0) {
      best_count = count;
      best_row = pivot_row;
      best_col = c;
    }
  }
  return {best_row, best_col};
}

// --- Double zeroColumn implementation ---

static constexpr std::ptrdiff_t W = simd::Width<double>;
// double and int64_t have the same width
static_assert(W == simd::Width<std::int64_t>);
static constexpr std::ptrdiff_t B = W;

template <std::ptrdiff_t BB>
auto collectBatchDouble(PtrMatrix<double> A, std::ptrdiff_t num_rows,
                        std::ptrdiff_t pivot_row, Col<> col_k,
                        std::ptrdiff_t &start_row,
                        std::array<std::ptrdiff_t, BB> &rows,
                        std::array<double, BB> &coeffs)
  -> std::array<std::ptrdiff_t, 2> {
  std::ptrdiff_t count = 0, zero_idx = -1;
  std::ptrdiff_t i = start_row;
  for (; i < num_rows && count < BB; ++i) {
    if (i == pivot_row) continue;
    double Aik = A[i, col_k];
    if (Aik == 0.0) continue;
    if (i == 0) zero_idx = count;
    rows[count] = i;
    coeffs[count] = Aik;
    ++count;
  }
  start_row = i;
  return {count, zero_idx};
}

// Scalar zeroWithRowOp for double
void zeroWithRowOpDouble(MutPtrMatrix<double> A, Row<> i, Row<> j, Col<> k) {
  double Aik = A[i, k];
  if (Aik == 0.0) return;
  double Ajk = A[j, k];

  // Convert to int64 for GCD
  auto iAik = static_cast<std::int64_t>(Aik);
  auto iAjk = static_cast<std::int64_t>(Ajk);
  std::int64_t g = gcd(iAik, iAjk);
  if (g != 1) {
    Aik /= static_cast<double>(g);
    Ajk /= static_cast<double>(g);
  }

  simd::Vec<W, double> vAjk = simd::vbroadcast<W, double>(Ajk),
                       vAik = simd::vbroadcast<W, double>(Aik);

  // GCD tracking in int64
  static constexpr simd::Vec<W, std::int64_t> ione =
    simd::Vec<W, std::int64_t>{} + 1;
  simd::Vec<W, std::int64_t> vg{static_cast<std::int64_t>(Ajk)};

  PtrMatrix<double> C = A; // const ref
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol()), l = 0;
  if (Ajk != 1.0) {
    for (;;) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      // fnma: c - a*b = Ajk*C[i,u] - Aik*C[j,u]
      simd::Vec<W, double> Ail =
        simd::fnma(vAik, C[j, u].vec_, vAjk * C[i, u].vec_);
      A[i, u] = Ail;
      // Convert to int64 for GCD tracking
      simd::Vec<W, std::int64_t> iAil =
        __builtin_convertvector(Ail, simd::Vec<W, std::int64_t>);
      vg = gcd<W>(iAil, vg);
      l += W;
      if (!simd::cmp::gt<W, std::int64_t>(vg, ione).any()) break;
    }
  }
  if (l < L) {
    for (;; l += W) {
      auto u{simd::index::unrollmask<1, W>(L, l)};
      if (!u) break;
      A[i, u] = simd::fnma(vAik, C[j, u].vec_, vAjk * C[i, u].vec_);
    }
  } else if (simd::cmp::gt<W, std::int64_t>(vg, ione).any()) {
    if (g = gcdreduce<W>(vg); g > 1) {
      double dg = static_cast<double>(g);
      for (std::ptrdiff_t ll = 0; ll < L; ++ll)
        if (A[i, ll] != 0.0) A[i, ll] /= dg;
    }
  }
}

// Batched zeroWithRowOp for double
template <std::ptrdiff_t BB>
void zeroWithRowOpDouble(MutPtrMatrix<double> A,
                         std::array<std::ptrdiff_t, BB> rows,
                         std::array<double, BB> coeffs, Row<> j, Col<> k) {
  static_assert(BB == W, "Batch size must equal SIMD width");
  static constexpr simd::Vec<W, std::int64_t> ione =
    simd::Vec<W, std::int64_t>{} + 1;

  double Ajk = A[j, k];

  // Convert coefficients to int64 for GCD pre-computation
  simd::Vec<W, std::int64_t> v_icoeffs;
  simd::Vec<W, double> v_dcoeffs;
  std::memcpy(&v_dcoeffs, coeffs.data(), sizeof(coeffs));
  v_icoeffs = __builtin_convertvector(v_dcoeffs, simd::Vec<W, std::int64_t>);

  auto iAjk = static_cast<std::int64_t>(Ajk);
  simd::Vec<W, std::int64_t> v_iajk = simd::vbroadcast<W, std::int64_t>(iAjk);

  // Parallel GCD
  simd::Vec<W, std::int64_t> v_gcd = gcd<W>(v_icoeffs, v_iajk);
  auto needs_reduce = simd::cmp::ne<W, std::int64_t>(v_gcd, ione);
  simd::Vec<W, std::int64_t> v_iaik_reduced =
    needs_reduce.select(v_icoeffs / v_gcd, v_icoeffs);
  simd::Vec<W, std::int64_t> v_iajk_reduced =
    needs_reduce.select(v_iajk / v_gcd, v_iajk);

  // Convert back to double for row ops
  simd::Vec<W, double> v_daik_reduced =
    __builtin_convertvector(v_iaik_reduced, simd::Vec<W, double>);
  simd::Vec<W, double> v_dajk_reduced =
    __builtin_convertvector(v_iajk_reduced, simd::Vec<W, double>);

  // Broadcast per-row coefficients
  std::array<simd::Vec<W, double>, BB> vAjk_d, vAik_d;
  std::array<simd::Vec<W, std::int64_t>, BB> vg;
#pragma clang loop unroll(full)
  for (std::ptrdiff_t b = 0; b < BB; ++b) {
    vAjk_d[b] = simd::vbroadcast<W, double>(v_dajk_reduced[b]);
    vAik_d[b] = simd::vbroadcast<W, double>(v_daik_reduced[b]);
    vg[b] = simd::Vec<W, std::int64_t>{v_iajk_reduced[b]};
  }

  PtrMatrix<double> C = A; // const ref
  std::ptrdiff_t L = std::ptrdiff_t(A.numCol());
  std::ptrdiff_t global_l = 0;

  // Active mask for GCD tracking
  auto active_mask = simd::cmp::ne<W, std::int64_t>(v_iajk_reduced, ione);

  // Main loop with GCD tracking
  if (active_mask.any()) {
    for (; global_l < L; global_l += W) {
      auto u{simd::index::unrollmask<1, W>(L, global_l)};
      if (!u) break;

      simd::Vec<W, double> Ajl = C[j, u].vec_;

#pragma clang loop unroll(full)
      for (std::ptrdiff_t b = 0; b < BB; ++b) {
        std::ptrdiff_t row_i = rows[b];
        // fnma: Ajk*C[row_i,u] - Aik*Ajl
        simd::Vec<W, double> Ail =
          simd::fnma(vAik_d[b], Ajl, vAjk_d[b] * C[row_i, u].vec_);
        A[row_i, u] = Ail;

        if (active_mask.intmask() & (1ULL << b)) {
          simd::Vec<W, std::int64_t> iAil =
            __builtin_convertvector(Ail, simd::Vec<W, std::int64_t>);
          vg[b] = gcd<W>(iAil, vg[b]);
          if (!simd::cmp::gt<W, std::int64_t>(vg[b], ione).any())
            active_mask =
              decltype(active_mask){active_mask.intmask() & ~(1ULL << b)};
        }
      }
      if (!active_mask.any()) {
        global_l += W;
        break;
      }
    }
  }

  // Continue without GCD tracking
  for (; global_l < L; global_l += W) {
    auto u{simd::index::unrollmask<1, W>(L, global_l)};
    if (!u) break;

    simd::Vec<W, double> Ajl = C[j, u].vec_;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t b = 0; b < BB; ++b) {
      std::ptrdiff_t row_i = rows[b];
      A[row_i, u] = simd::fnma(vAik_d[b], Ajl, vAjk_d[b] * C[row_i, u].vec_);
    }
  }

  // Post-loop GCD reduction
  std::uint64_t remaining_mask = active_mask.intmask();
#pragma clang loop unroll(full)
  for (std::ptrdiff_t b = 0; b < BB; ++b) {
    if ((remaining_mask & (1ULL << b)) &&
        simd::cmp::gt<W, std::int64_t>(vg[b], ione).any()) {
      std::int64_t g = gcdreduce<W>(vg[b]);
      if (g > 1) {
        double dg = static_cast<double>(g);
        std::ptrdiff_t row_i = rows[b];
        for (std::ptrdiff_t ll = 0; ll < L; ++ll)
          if (A[row_i, ll] != 0.0) A[row_i, ll] /= dg;
      }
    }
  }
}

// Driver: zeroColumn for double matrices
void zeroColumnDouble(MutPtrMatrix<double> A, Row<> pivot_row,
                      Col<> pivot_col) {
  std::ptrdiff_t M = std::ptrdiff_t(A.numRow());
  std::ptrdiff_t pivot = std::ptrdiff_t(pivot_row);

  std::array<std::ptrdiff_t, B> batch_rows;
  std::array<double, B> batch_coeffs;
  std::ptrdiff_t i = 0;

  while (i < M) {
    auto [count, zero_idx] = collectBatchDouble<B>(A, M, pivot, pivot_col, i,
                                                   batch_rows, batch_coeffs);
    if (count == 0) break;
    if (count == B) {
      zeroWithRowOpDouble<B>(A, batch_rows, batch_coeffs, pivot_row, pivot_col);
    } else {
      for (std::ptrdiff_t j = 0; j < count; ++j)
        zeroWithRowOpDouble(A, row(batch_rows[j]), pivot_row, pivot_col);
      break;
    }
  }
}

} // namespace

void BM_ZeroCol_Int64_0(Bench &bench) {
  auto tableau =
    utils::readMatrixBinary(MATH_DATA_DIR "/simplex_tableau_0.binmat");
  auto [pivot_row, pivot_col] = findPivot(tableau);

  math::DenseMatrix<std::int64_t> work(shape(tableau));
  bench.run("zeroColumn int64 (tableau_0)", [&] {
    work << tableau;
    math::NormalForm::zeroColumn(work, math::row(pivot_row),
                                 math::col(pivot_col));
    doNotOptimizeAway(work);
  });
}

void BM_ZeroCol_Double_0(Bench &bench) {
  auto tableau =
    utils::readMatrixBinary(MATH_DATA_DIR "/simplex_tableau_0.binmat");
  std::ptrdiff_t M = std::ptrdiff_t(tableau.numRow()),
                 N = std::ptrdiff_t(tableau.numCol());
  auto [pivot_row, pivot_col] = findPivot(tableau);

  // Convert to double
  math::DenseMatrix<double> backup(math::DenseDims(math::row(M), math::col(N)));
  for (std::ptrdiff_t r = 0; r < M; ++r)
    for (std::ptrdiff_t c = 0; c < N; ++c)
      backup[r, c] = static_cast<double>(tableau[r, c]);

  math::DenseMatrix<double> work(math::DenseDims(math::row(M), math::col(N)));
  bench.run("zeroColumn double (tableau_0)", [&] {
    work << backup;
    zeroColumnDouble(work, math::row(pivot_row), math::col(pivot_col));
    doNotOptimizeAway(work);
  });
}

void BM_ZeroCol_Int64_1(Bench &bench) {
  auto tableau =
    utils::readMatrixBinary(MATH_DATA_DIR "/simplex_tableau_1.binmat");
  auto [pivot_row, pivot_col] = findPivot(tableau);

  math::DenseMatrix<std::int64_t> work(shape(tableau));
  bench.run("zeroColumn int64 (tableau_1)", [&] {
    work << tableau;
    math::NormalForm::zeroColumn(work, math::row(pivot_row),
                                 math::col(pivot_col));
    doNotOptimizeAway(work);
  });
}

void BM_ZeroCol_Double_1(Bench &bench) {
  auto tableau =
    utils::readMatrixBinary(MATH_DATA_DIR "/simplex_tableau_1.binmat");
  std::ptrdiff_t M = std::ptrdiff_t(tableau.numRow()),
                 N = std::ptrdiff_t(tableau.numCol());
  auto [pivot_row, pivot_col] = findPivot(tableau);

  // Convert to double
  math::DenseMatrix<double> backup(math::DenseDims(math::row(M), math::col(N)));
  for (std::ptrdiff_t r = 0; r < M; ++r)
    for (std::ptrdiff_t c = 0; c < N; ++c)
      backup[r, c] = static_cast<double>(tableau[r, c]);

  math::DenseMatrix<double> work(math::DenseDims(math::row(M), math::col(N)));
  bench.run("zeroColumn double (tableau_1)", [&] {
    work << backup;
    zeroColumnDouble(work, math::row(pivot_row), math::col(pivot_col));
    doNotOptimizeAway(work);
  });
}

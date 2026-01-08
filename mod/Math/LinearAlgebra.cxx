module LinearAlgebra;

namespace math::LU {
auto ldivrat(SquarePtrMatrix<Rational> F, PtrVector<unsigned> ipiv,
             MutPtrMatrix<Rational> rhs) -> bool {
  auto [M, N] = shape(rhs);
  invariant(std::ptrdiff_t(F.numRow()), std::ptrdiff_t(M));
  // permute rhs
  for (std::ptrdiff_t i = 0; i < M; ++i)
    if (std::ptrdiff_t ip = ipiv[i]; i != ip)
      for (std::ptrdiff_t j = 0; j < M; ++j)
        std::ranges::swap(rhs[ip, j], rhs[i, j]);

  // LU x = rhs
  // L y = rhs // L is UnitLowerTriangular
  for (std::ptrdiff_t m = 0; m < M; ++m) {
    for (std::ptrdiff_t n = 0; n < N; ++n) {
      Rational Ymn = rhs[m, n];
      for (std::ptrdiff_t k = 0; k < m; ++k)
        if (Ymn.fnmadd(F[m, k], rhs[k, n])) return true;
      rhs[m, n] = Ymn;
    }
  }
  // U x = y
  for (std::ptrdiff_t m = M; m--;) {
    for (std::ptrdiff_t n = 0; n < N; ++n) {
      Rational Ymn = rhs[m, n];
      for (std::ptrdiff_t k = m + 1; k < M; ++k)
        if (Ymn.fnmadd(F[m, k], rhs[k, n])) return true;
      if (auto div = Ymn.safeDiv(F[m, m])) rhs[m, n] = *div;
      else return true;
    }
  }
  return false;
}

auto rdivrat(SquarePtrMatrix<Rational> F, PtrVector<unsigned> ipiv,
             MutPtrMatrix<Rational> rhs) -> bool {
  auto [M, N] = shape(rhs);
  invariant(std::ptrdiff_t(F.numCol()), std::ptrdiff_t(N));
  // PA = LU
  // x LU = rhs
  // y U = rhs
  for (std::ptrdiff_t m = 0; m < M; ++m) {
    for (std::ptrdiff_t n = 0; n < N; ++n) {
      Rational Ymn = rhs[m, n];
      for (std::ptrdiff_t k = 0; k < n; ++k)
        if (Ymn.fnmadd(rhs[m, k], F[k, n])) return true;
      if (auto div = Ymn.safeDiv(F[n, n])) rhs[m, n] = *div;
      else return true;
    }
  }
  // x L = y
  for (std::ptrdiff_t m = 0; m < M; ++m) {
    for (std::ptrdiff_t n = N; n--;) {
      Rational Xmn = rhs[m, n];
      for (std::ptrdiff_t k = n + 1; k < N; ++k)
        if (Xmn.fnmadd(rhs[m, k], F[k, n])) return true;
      rhs[m, n] = Xmn;
    }
  }
  // permute rhs
  for (std::ptrdiff_t j = N; j--;)
    if (std::ptrdiff_t jp = ipiv[j]; j != jp)
      for (std::ptrdiff_t i = 0; i < M; ++i) std::swap(rhs[i, jp], rhs[i, j]);

  return false;
}
} // namespace math::LU

using namespace math;
auto factImplRat(MutSquarePtrMatrix<Rational> A, MutPtrVector<unsigned> ipiv)
  -> bool {
  Row M = A.numRow();
  invariant(ipiv.size(), std::ptrdiff_t(M));
  for (std::ptrdiff_t k = 0;; ++k) {
    std::ptrdiff_t kp = k;
    for (;; ++kp) {
      if (kp == M) return true;
      if (A[kp, k] == 0) continue;
      ipiv[k] = kp;
      break;
    }
    if (kp != k)
      for (std::ptrdiff_t j = 0; j < M; ++j) std::swap(A[kp, j], A[k, j]);
    if (k + 1 == M) break;
    Rational invAkk = A[k, k].inv();
    for (std::ptrdiff_t i = k + 1; i < M; ++i)
      if (std::optional<Rational> Aik = A[i, k].safeMul(invAkk)) A[i, k] = *Aik;
      else return true;
    for (std::ptrdiff_t i = k + 1; i < M; ++i) {
      for (std::ptrdiff_t j = k + 1; j < M; ++j) {
        if (std::optional<Rational> kAij = A[i, k].safeMul(A[k, j])) {
          if (std::optional<Rational> Aij = A[i, j].safeSub(*kAij)) {
            A[i, j] = *Aij;
            continue;
          }
        }
        return true;
      }
    }
  }
  return false;
}

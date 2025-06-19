#ifdef USE_MODULE
module;
#else
#pragma once
#endif
#ifndef USE_MODULE
#include "Containers/Tuple.cxx"
#include "Math/Array.cxx"
#include "Math/ArrayConcepts.cxx"
#include "Math/AxisTypes.cxx"
#include "Math/GenericConstructors.cxx"
#include "Math/Indexing.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/Rational.cxx"
#include "Utilities/ArrayPrint.cxx"
#include "Utilities/Parameters.cxx"
#include "Utilities/TypeCompression.cxx"
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#else
export module LinearAlgebra;

import Array;
import ArrayConcepts;
import ArrayConstructors;
import AxisTypes;
import GenericArrayConstructors;
import ManagedArray;
import Param;
import Rational;
import std;
import Tuple;
import TypeCompression;
#endif

template <typename T>
concept TrivialVec = utils::TriviallyCopyable<T> && math::AbstractVector<T>;
// template <typename T>
// concept TrivialMat = utils::TriviallyCopyable<T> && AbstractMatrix<T>;
//
#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif
namespace LU {
[[nodiscard]] constexpr auto ldivrat(SquarePtrMatrix<Rational> F,
                                     PtrVector<unsigned> ipiv,
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
template <class S>
constexpr void ldiv(SquarePtrMatrix<S> F, PtrVector<unsigned> ipiv,
                    MutPtrMatrix<S> R) {
  auto [M, N] = shape(R);
  invariant(std::ptrdiff_t(F.numRow()), M);
  invariant(M > 0);
  // permute rhs
  for (std::ptrdiff_t i = 0; i < M; ++i)
    if (std::ptrdiff_t ip = ipiv[i]; i != ip)
      for (std::ptrdiff_t j = 0; j < M; ++j)
        std::ranges::swap(R[ip, j], R[i, j]);

  // LU x = rhs
  // L y = rhs // L is UnitLowerTriangular
  for (std::ptrdiff_t m = 1; m < M; ++m)
    R[m, _] -= F[m, _(0, m)] * R[_(0, m), _];
  // U x = y
  R[last, _] /= F[last, last];
  for (std::ptrdiff_t m = M - 1; m--;)
    R[m, _] << (R[m, _] - F[m, _(m + 1, end)] * R[_(m + 1, end), _]) / F[m, m];
}

[[nodiscard]] constexpr auto rdivrat(SquarePtrMatrix<Rational> F,
                                     PtrVector<unsigned> ipiv,
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
template <class S>
constexpr void rdiv(SquarePtrMatrix<S> F, PtrVector<unsigned> ipiv,
                    MutPtrMatrix<S> rhs) {
  auto [M, N] = shape(rhs);
  invariant(std::ptrdiff_t(F.numCol()), N);
  invariant(N > 0);
  // PA = LU
  // x LU = rhs
  // y U = rhs
  rhs[_, 0] /= F[0, 0];
  for (std::ptrdiff_t n = 1; n < N; ++n)
    rhs[_, n] << (rhs[_, n] - rhs[_, _(0, n)] * F[_(0, n), n]) / F[n, n];
  // x L = y
  for (std::ptrdiff_t n = N - 1; n--;)
    rhs[_, n] -= rhs[_, _(n + 1, end)] * F[_(n + 1, end), n];
  // permute rhs
  for (std::ptrdiff_t j = N; j--;)
    if (std::ptrdiff_t jp = ipiv[j]; j != jp)
      for (std::ptrdiff_t i = 0; i < M; ++i)
        std::ranges::swap(rhs[i, jp], rhs[i, j]);
}

template <class T, std::ptrdiff_t L> class Fact {
  SquareMatrix<T, L> F;
  Vector<unsigned> ipiv;

public:
  constexpr void ldiv(MutPtrMatrix<T> rhs) const { LU::ldiv(F, ipiv, rhs); }
  constexpr void rdiv(MutPtrMatrix<T> rhs) const { LU::rdiv(F, ipiv, rhs); }
  constexpr auto ldivrat(MutPtrMatrix<T> rhs) const -> bool {
    return LU::ldivrat(F, ipiv, rhs);
  }
  constexpr auto rdivrat(MutPtrMatrix<T> rhs) const -> bool {
    return LU::rdivrat(F, ipiv, rhs);
  }
  constexpr Fact(SquareMatrix<T, L> f, Vector<unsigned> ip)
    : F(std::move(f)), ipiv(std::move(ip)) {
    invariant(std::ptrdiff_t(F.numRow()), ipiv.size());
  }

  [[nodiscard]] constexpr auto inv() const
    -> std::optional<SquareMatrix<Rational, L>> {
    SquareMatrix<Rational, L> A{
      SquareMatrix<Rational, L>::identity(std::ptrdiff_t(F.numCol()))};
    if (!ldivrat(A)) return A;
    return {};
  }
  [[nodiscard]] constexpr auto det() const -> std::optional<Rational> {
    Rational d = F(0, 0);
    for (std::ptrdiff_t i = 1; i < F.numCol(); ++i)
      if (auto di = d.safeMul(F(i, i))) d = *di;
      else return {};
    return d;
  }
  [[nodiscard]] constexpr auto perm() const -> Vector<unsigned> {
    Col M = F.numCol();
    Vector<unsigned> perm{M};
    for (std::ptrdiff_t m = 0; m < M; ++m) perm[m] = m;
    for (std::ptrdiff_t m = 0; m < M; ++m) std::swap(perm[m], perm[ipiv[m]]);
    return perm;
  }
  void print() const {
    utils::print("LU fact:\n");
    F.print();
    utils::print("\nperm = \n");
    utils::printVector(ipiv.begin(), ipiv.end());
    utils::print('\n');
  }
};
template <std::ptrdiff_t L>
[[nodiscard]] constexpr auto fact(const SquareMatrix<std::int64_t, L> &B)
  -> std::optional<Fact<Rational, L>> {
  Row M = B.numRow();
  SquareMatrix<Rational, L> A{B};
  // auto ipiv = Vector<unsigned>{.s = unsigned(M)};
  auto ipiv{vector(math::DefaultAlloc<unsigned>{}, std::ptrdiff_t(M))};
  // Vector<unsigned> ipiv{.s = unsigned(M)};
  invariant(ipiv.size(), std::ptrdiff_t(M));
  for (std::ptrdiff_t k = 0;; ++k) {
    std::ptrdiff_t kp = k;
    for (;; ++kp) {
      if (kp == M) return {};
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
      else return {};
    for (std::ptrdiff_t i = k + 1; i < M; ++i) {
      for (std::ptrdiff_t j = k + 1; j < M; ++j) {
        if (std::optional<Rational> kAij = A[i, k].safeMul(A[k, j])) {
          if (std::optional<Rational> Aij = A[i, j].safeSub(*kAij)) {
            A[i, j] = *Aij;
            continue;
          }
        }
        return {};
      }
    }
  }
  return Fact<Rational, L>{std::move(A), std::move(ipiv)};
}
template <typename S> constexpr auto factImpl(MutSquarePtrMatrix<S> A) {
  using V = decltype(extractvalue(S{}));
  Row M = A.numRow();
  auto ipiv{vector(math::DefaultAlloc<unsigned>{}, std::ptrdiff_t(M))};
  invariant(std::ptrdiff_t(ipiv.size()), std::ptrdiff_t(M));
  for (std::ptrdiff_t k = 0;; ++k) {
    containers::Pair<std::ptrdiff_t, V> mi{-1, {}};
    for (std::ptrdiff_t i = k; i < M; ++i)
      if (V v = std::abs(extractvalue(A[i, k])); v > mi._1) mi = {i, v};
    invariant(mi._0 >= 0); // TODO: return info?
    ipiv[k] = mi._0;
    if (mi._0 != k)
      for (std::ptrdiff_t j = 0; j < M; ++j)
        std::ranges::swap(A[mi._0, j], A[k, j]);
    if (k + 1 == M) break;
    S invAkk = 1.0 / A[k, k];
    for (std::ptrdiff_t i = k + 1; i < M; ++i)
      A[i, _(k + 1, end)] -= (A[i, k] *= invAkk) * A[k, _(k + 1, end)];
  }
  return ipiv;
}
template <class S, std::ptrdiff_t L>
[[nodiscard]] constexpr auto fact(SquareMatrix<S, L> A) -> Fact<S, L> {
  auto &&ipiv{factImpl(A)};
  return Fact<S, L>{std::move(A), std::move(ipiv)};
}
/// ldiv(A, B)
/// computes A \ B, modifying A and B
/// Note that this computes an LU factorization;
/// if you are performing more than one division,
/// it would be more efficient to precompute an
/// `auto F = LU::fact(A)`, and use this for multiple
/// `F.ldiv(B)` calls.
template <typename T>
constexpr void ldiv(MutSquarePtrMatrix<T> A, MutPtrMatrix<T> B) {
  auto ipiv{factImpl(A)};
  ldiv(A, ipiv, B);
}
/// rdiv(A, B)
/// Computes B / A, modifying A and B
/// Note that this computes an LU factorization;
/// if you are performing more than one division,
/// it would be more efficient to precompute an
/// `auto F = LU::fact(A)`, and use this for multiple
/// `F.rdiv(B)` calls.
template <typename T>
constexpr void rdiv(MutSquarePtrMatrix<T> A, MutPtrMatrix<T> B) {
  auto ipiv{factImpl(A)};
  rdiv(A, ipiv, B);
}

} // namespace LU

/// factorizes symmetric full-rank (but not necessarily positive-definite)
/// matrix A into LD^-1L', where L is lower-triangular with 1s on the
/// diagonal
/// Only uses the lower triangle of A, overwriting it.
/// `D` is stored into the diagonal of `A`.
namespace LDL {

/// NOT OWNING
/// TODO: make the API consistent between LU and LDL
template <typename T> class Fact {
  MutSquarePtrMatrix<T> F;

public:
  constexpr Fact(MutSquarePtrMatrix<T> A) : F{A} {};

  constexpr void ldiv(MutPtrMatrix<T> R) {
    std::ptrdiff_t M = std::ptrdiff_t(R.numRow());
    invariant(std::ptrdiff_t(F.numRow()), M);
    // LD^-1L' x = rhs
    // L y = rhs // L is UnitLowerTriangular
    for (std::ptrdiff_t m = 1; m < M; ++m)
      R[m, _] -= F[m, _(0, m)] * R[_(0, m), _];
    // D^-1 L' x = y
    // L' x = D y
    R[last, _] *= F[last, last];
    for (std::ptrdiff_t m = M - 1; m--;)
      R[m, _] << R[m, _] * F[m, m] - F[_(m + 1, M), m].t() * R[_(m + 1, M), _];
  }
  constexpr void ldiv(MutPtrVector<T> R) {
    std::ptrdiff_t M = R.size();
    invariant(std::ptrdiff_t(F.numRow()), M);
    // LD^-1L' x = rhs
    // L y = rhs // L is UnitLowerTriangular
    for (std::ptrdiff_t m = 1; m < M; ++m)
      R[m] -= R[_(0, m)] * F[m, _(0, m)].t();
    // D^-1 L' x = y
    // L' x = D y
    R[last] *= F[last, last];
    for (std::ptrdiff_t m = M - 1; m--;)
      R[m] = R[m] * F[m, m] - R[_(m + 1, M)] * F[_(m + 1, M), m];
  }
  constexpr void ldiv(MutPtrVector<T> dst, TrivialVec auto src) {
    std::ptrdiff_t M = dst.size();
    invariant(M, std::ptrdiff_t(src.size()));
    invariant(std::ptrdiff_t(F.numRow()), M);
    // LD^-1L' x = rhs
    // L y = rhs // L is UnitLowerTriangular
    dst[0] = src[0];
    for (std::ptrdiff_t m = 1; m < M; ++m)
      dst[m] = src[m] - dst[_(0, m)] * F[m, _(0, m)].t();
    // D^-1 L' x = y
    // L' x = D y
    dst[last] *= F[last, last];
    for (std::ptrdiff_t m = M - 1; m--;)
      dst[m] = dst[m] * F[m, m] - dst[_(m + 1, M)] * F[_(m + 1, M), m];
  }
};

template <bool ForcePD = false, typename T>
constexpr auto factorize(MutSquarePtrMatrix<T> A) -> Fact<T> {
  Row M = A.numRow();
  invariant(std::ptrdiff_t(M), std::ptrdiff_t(A.numCol()));
  for (std::ptrdiff_t k = 0;; ++k) {
    T Akk = A[k, k];
    if constexpr (ForcePD) Akk = std::max(Akk, T(0.000001));
    T invAkk = A[k, k] = 1.0 / Akk;
    if (k + 1 == M) break;
    A[_(k + 1, M), k] *= invAkk;
    for (std::ptrdiff_t i = k + 1; i < M; ++i)
      A[i, _(k + 1, i + 1)] -= (A[i, k] * Akk) * A[_(k + 1, i + 1), k].t();
  }
  return Fact{A};
}

template <bool ForcePD = false, typename T>
constexpr void ldiv(MutSquarePtrMatrix<T> A, MutPtrMatrix<T> B) {
  factorize<ForcePD>(A).ldiv(B);
}
template <bool ForcePD = false, typename T>
constexpr void ldiv(MutSquarePtrMatrix<T> A, MutPtrVector<T> B) {
  factorize<ForcePD>(A).ldiv(B);
}
template <bool ForcePD = false, typename T>
constexpr void ldiv(MutSquarePtrMatrix<T> A, MutPtrVector<T> B,
                    TrivialVec auto C) {
  factorize<ForcePD>(A).ldiv(B, C);
}

} // namespace LDL
} // namespace math

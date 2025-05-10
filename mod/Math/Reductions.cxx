#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#include "LoopMacros.hxx"

#ifndef USE_MODULE
#include "Math/ArrayConcepts.cxx"
#include <cstddef>
#else
export module Reductions;
import ArrayConcepts;
import std;
#endif

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif
constexpr auto abs2(auto x) { return x * x; }
template <AbstractTensor B> constexpr auto norm2(const B &A) {
  utils::eltype_t<B> s = 0;
  if constexpr (!LinearlyIndexable<B, utils::eltype_t<B>>) {
    for (std::ptrdiff_t i = 0; i < A.numRow(); ++i) {
      for (std::ptrdiff_t j = 0; j < A.numCol(); ++j) {
        POLYMATHFAST
        s += abs2(A[i, j]);
      }
    }
  } else
    for (std::ptrdiff_t j = 0, L = std::ptrdiff_t(A.size()); j < L; ++j) {
      POLYMATHFAST
      s += abs2(A[j]);
    }
  return s;
}

constexpr auto norm2(const auto &a) {
  decltype(a[0] * a[0] + a[1] * a[1]) s{};
  for (auto x : a) {
    POLYMATHFAST
    s += abs2(x, x);
  }
  return s;
}
// TODO: SIMD?
constexpr auto dot(const auto &a, const auto &b) {
  std::ptrdiff_t L = a.size();
  invariant(L, b.size());
  decltype(a[0] * b[0] + a[1] * b[1]) s{};
  for (std::ptrdiff_t i = 0; i < L; ++i) {
    POLYMATHFAST
    s += a[i] * b[i];
  }
  return s;
}
} // namespace math

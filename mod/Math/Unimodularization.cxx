#ifdef USE_MODULE
module;
#else
#pragma once
#endif
#ifndef USE_MODULE
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>

#include "Containers/Pair.cxx"
#include "Math/Array.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/MatrixDimensions.cxx"
#include "Math/NormalForm.cxx"
#else
export module Unimodularization;

import Array;
import ManagedArray;
import MatDim;
import NormalForm;
import Pair;
import std;
#endif

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif
// if `A` can be unimodularized, returns the inverse of the unimodularized `A`
[[nodiscard]] inline auto unimodularize(IntMatrix<> A)
  -> std::optional<SquareMatrix<std::int64_t>> {
  SquareMatrix<std::int64_t> U{SquareDims{A.numRow()}};
  NormalForm::hermite(A, U);
  for (std::ptrdiff_t m = 0; m < A.numCol(); ++m)
    if (A[m, m] != 1) return {};
  return U;
}
} // namespace math

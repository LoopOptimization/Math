#include <gtest/gtest.h>
#ifndef USE_MODULE
#include "Alloc/Arena.cxx"
#include "Macros.hxx"
#include "Math/Array.cxx"
#include "Math/AxisTypes.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/MatrixDimensions.cxx"
#include "Math/NormalForm.cxx"
#include <cstddef>
#include <cstdint>
#include <random>
#else

import Arena;
import Array;
import ManagedArray;
import MatDim;
import NormalForm;
import std;
#endif

using math::DenseMatrix, math::DenseDims, math::row, math::col;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(OrthogonalizeMatricesTest, BasicAssertions) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(-2, 2);

  constexpr std::ptrdiff_t M = 7;
  constexpr std::ptrdiff_t N = 7;
  std::int64_t mem[M * N * 2];
  alloc::OwningArena<> alloc;
  math::MutDensePtrMatrix<std::int64_t> A{mem, DenseDims<>{row(M), col(N)}};
  constexpr std::size_t iters = 1000;
  for (std::size_t i = 0; i < iters; ++i) {
    for (auto &&a : A) a = distrib(gen);
    std::ptrdiff_t r = math::NormalForm::rank(alloc, A);
    auto orth = math::orthogonalize(&alloc, A);
    ASSERT_EQ(orth.numCol(), A.numCol());
    ASSERT_EQ(orth.numRow(), r);
    math::MutDensePtrMatrix<std::int64_t> B{mem + (M * N),
                                            DenseDims<>{row(r), col(r)}};
    // note, A'A is not diagonal
    // but AA' is
    B << orth * orth.t();
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
    for (std::ptrdiff_t m = 0; m < r; ++m)
      for (std::ptrdiff_t n = 0; n < r; ++n)
        if (m != n) ASSERT_EQ((B[m, n]), 0);
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  }
}

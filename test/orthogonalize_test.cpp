#include <gtest/gtest.h>
#ifndef USE_MODULE
#include "Alloc/Arena.cxx"
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
import STL;
#endif

using math::DenseMatrix, math::DenseDims, math::row, math::col;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(OrthogonalizeMatricesTest, BasicAssertions) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(-3, 3);

  constexpr ptrdiff_t M = 7;
  constexpr ptrdiff_t N = 7;
  int64_t mem[M * N * 2];
  alloc::OwningArena<> alloc;
  math::MutDensePtrMatrix<int64_t> A{mem,
                                     DenseDims<>{math::row(M), math::col(N)}};
  math::MutDensePtrMatrix<int64_t> B{mem + (M * N),
                                     DenseDims<>{math::row(M), math::col(N)}};
  constexpr size_t iters = 1000;
  for (size_t i = 0; i < iters; ++i) {
    for (auto &&a : A) a = distrib(gen);
    // std::cout << "Random A =\n" << A << "\n";
    math::orthogonalize(&alloc, A);
    // std::cout << "Orthogonal A =\n" << A << "\n";
    // note, A'A is not diagonal
    // but AA' is
    B << A * A.t();
    // std::cout << "A'A =\n" << B << "\n";
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
    for (size_t m = 0; m < M; ++m)
      for (size_t n = 0; n < N; ++n)
        if (m != n) EXPECT_EQ((B[m, n]), 0);
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  }
}

#include <gtest/gtest.h>
import Arena;
import Array;
import ArrayConcepts;
import ArrayConstructors;
import ArrayParse;
import Comparisons;
import CorePrint;
import LinearAlgebra;
import ManagedArray;
import MatDim;
import NormalForm;
import std;
import UniformScaling;

using namespace math;
using utils::operator""_mat;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(OrthogonalizationTest, BasicAssertions) {
  SquareMatrix<std::int64_t> A(SquareDims<>{math::row(4)});
  utils::print("\n\n\n========\n========\n========\n\n");
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(-10, 10);
  std::ptrdiff_t orth_any_count = 0;
  std::ptrdiff_t orth_max_count = 0;
  std::ptrdiff_t orth_count = 0;
  std::ptrdiff_t lu_failed_count = 0;
  std::ptrdiff_t inv_failed_count = 0;
  std::ptrdiff_t num_iters = 1000;
  IntMatrix<> B(DenseDims<>{row(4), col(8)});
  SquareMatrix<std::int64_t> I4 = SquareMatrix<std::int64_t>::identity(4);
  for (std::ptrdiff_t i = 0; i < num_iters; ++i) {
    for (std::ptrdiff_t n = 0; n < 4; ++n)
      for (std::ptrdiff_t m = 0; m < 8; ++m) B[n, m] = distrib(gen);
    auto [K, included] = NormalForm::orthogonalize(B);
    orth_count += included.size();
    orth_any_count += (!included.empty());
    orth_max_count += (included.size() == 4);
    if (included.size() == 4) {
      for (std::ptrdiff_t n = 0; n < 4; ++n) {
        std::ptrdiff_t m = 0;
        for (auto mb : included) A[n, m++] = B[n, mb];
      }
      utils::print("K=\n");
      K.print();
      utils::print('\n');
      utils::print("A=\n");
      A.print();
      utils::print('\n');
      EXPECT_EQ(K * A, I4);
      EXPECT_EQ(K * A, I);
    } else if (auto optlu = LU::fact(K)) {
      if (auto opt_a2 = (*optlu).inv()) {
        auto &A2 = *opt_a2;
        for (std::ptrdiff_t n = 0; n < 4; ++n)
          for (std::ptrdiff_t j = 0; j < included.size(); ++j)
            EXPECT_EQ((A2[n, j]), (B[n, included[j]]));
      } else {
        ++inv_failed_count;
      }
    } else {
      ++lu_failed_count;
      utils::print("B = ");
      B.print();
      utils::print("\nK = ");
      K.print();
      utils::print('\n');
      continue;
    }
  }
  utils::print("Mean orthogonalized: ");
  utils::print(double(orth_count) / double(num_iters));
  utils::print("\nOrthogonalization succeeded on at least one: ");
  utils::print(orth_any_count);
  utils::print(" / ");
  utils::print(num_iters);
  utils::print("\nOrthogonalization succeeded on 4: ");
  utils::print(orth_max_count);
  utils::print(" / ");
  utils::print(num_iters);
  utils::print("\nLU fact failed count: ");
  utils::print(lu_failed_count);
  utils::print(" / ");
  utils::print(num_iters);
  utils::print("\nInv fact failed count: ");
  utils::print(inv_failed_count);
  utils::print(" / ");
  utils::print(num_iters);
  utils::print('\n');

  B[0, 0] = 1;
  B[1, 0] = 0;
  B[2, 0] = 1;
  B[3, 0] = 0;
  B[0, 1] = 0;
  B[1, 1] = 1;
  B[2, 1] = 0;
  B[3, 1] = 1;
  B[0, 2] = 1;
  B[1, 2] = 0;
  B[2, 2] = 0;
  B[3, 2] = 0;
  B[0, 3] = 0;
  B[1, 3] = 1;
  B[2, 3] = 0;
  B[3, 3] = 0;
  B[0, 4] = 0;
  B[1, 4] = 0;
  B[2, 4] = 1;
  B[3, 4] = 0;
  B[0, 5] = 0;
  B[1, 5] = 0;
  B[2, 5] = 0;
  B[3, 5] = 1;
  utils::print("B_orth_motivating_example = ");
  B.print();
  utils::print('\n');
  auto [K, included] = NormalForm::orthogonalize(B);
  EXPECT_EQ(included.size(), 4);
  for (std::ptrdiff_t i = 0; i < 4; ++i) EXPECT_EQ(included[i], i);
  for (std::ptrdiff_t n = 0; n < 4; ++n) {
    std::ptrdiff_t m = 0;
    for (auto mb : included) {
      A[n, m] = B[n, mb];
      ++m;
    }
  }
  IntMatrix<> KA{K * A};
  utils::print("A = ");
  A.print();
  utils::print("\nA * K = ");
  KA.print();
  utils::print('\n');
  EXPECT_TRUE(KA == I4);
}

namespace {
auto isHNF(PtrMatrix<std::int64_t> A) -> bool {
  const auto [M, N] = shape(A);
  // l is lead
  Col<> l = {};
  for (std::ptrdiff_t m = 0; m < M; ++m) {
    // all entries must be 0
    for (std::ptrdiff_t n = 0; n < l; ++n)
      if (A[m, n]) return false;
    // now search for next lead
    while ((l < N) && A[m, l] == 0) ++l;
    if (l == N) continue;
    std::int64_t Aml = A[m, l];
    if (Aml < 0) return false;
    for (std::ptrdiff_t r = 0; r < m; ++r) {
      std::int64_t Arl = A[r, l];
      if ((Arl >= Aml) || (Arl < 0)) return false;
    }
  }
  return true;
}
} // namespace

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(Hermite, BasicAssertions) {
  {
    IntMatrix<> A43{DenseDims<>{row(4), col(3)}};
    A43[0, 0] = 2;
    A43[1, 0] = 3;
    A43[2, 0] = 6;
    A43[3, 0] = 2;
    A43[0, 1] = 5;
    A43[1, 1] = 6;
    A43[2, 1] = 1;
    A43[3, 1] = 6;
    A43[0, 2] = 8;
    A43[1, 2] = 3;
    A43[2, 2] = 1;
    A43[3, 2] = 1;
    utils::print("A=\n");
    A43.print();
    utils::print('\n');
    IntMatrix<> H = A43;
    SquareMatrix<std::int64_t> U{SquareDims{H.numRow()}};
    NormalForm::hermite(H, U);
    utils::print("H=\n");
    H.print();
    utils::print("\nU=\n");
    U.print();
    utils::print('\n');

    EXPECT_TRUE(isHNF(H));
    EXPECT_TRUE(H == U * A43);

    for (std::ptrdiff_t i = 0; i < 3; ++i) A43[2, i] = A43[0, i] + A43[1, i];
    utils::print("\n\n\n=======\n\nA=\n");
    A43.print();
    utils::print('\n');
    H << A43;
    NormalForm::hermite(H, U);
    utils::print("H=\n");
    H.print();
    utils::print("\nU=\n");
    U.print();
    utils::print('\n');
    EXPECT_TRUE(isHNF(H));
    EXPECT_TRUE(H == U * A43);
  }
  {
    SquareMatrix<std::int64_t> A(SquareDims<>{math::row(4)});
    A[0, 0] = 3;
    A[1, 0] = -6;
    A[2, 0] = 7;
    A[3, 0] = 7;
    A[0, 1] = 7;
    A[1, 1] = -8;
    A[2, 1] = 10;
    A[3, 1] = 6;
    A[0, 2] = -5;
    A[1, 2] = 8;
    A[2, 2] = 7;
    A[3, 2] = 3;
    A[0, 3] = -5;
    A[1, 3] = -6;
    A[2, 3] = 8;
    A[3, 3] = -1;
    IntMatrix<> H = A;
    SquareMatrix<std::int64_t> U{SquareDims{H.numRow()}};
    NormalForm::hermite(H, U);
    utils::print("\n\n\n====\n\nH=\n");
    H.print();
    utils::print("\nU=\n");
    U.print();
    utils::print('\n');
    EXPECT_TRUE(isHNF(H));
    EXPECT_TRUE(H == U * A);
  }
  {
    IntMatrix<> A{"[1 -3 0 -2 0 0 -1 -1 0 0 -1 0 0 0 0 0 0 "
                  "0 0 0 0 0; 0 1 0 1 0 0 0 1 0 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 "
                  "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0; 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 "
                  "0 0 0 0; 0 -1 1 -1 1 0 0 -1 1 "
                  "0 0 0 0 0 0 0 0 0 0 0 0 0; 0 -1 1 0 0 1 "
                  "-1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
                  "0 0; 0 -1 1 -1 1 0 0 0 0 1 -1 0 0 0 0 0 "
                  "0 0 0 0 0 0; -1 0 0 0 0 0 0 0 "
                  "0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 -1 0 0 0 "
                  "0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 "
                  "0 0; 0 0 -1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 "
                  "0 0 0 0 0; 0 0 0 -1 0 0 0 0 0 "
                  "0 0 0 0 0 1 0 0 0 0 0 0 0; 0 0 0 0 -1 0 "
                  "0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 "
                  "0; 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 1 0 "
                  "0 0 0 0; 0 0 0 0 0 0 -1 0 0 0 "
                  "0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 "
                  "-1 0 0 0 0 0 0 0 0 0 0 1 0 0 "
                  "0; 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 "
                  "0 1 0 0; 0 0 0 0 0 0 0 0 0 -1 "
                  "0 0 0 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 "
                  "0 0 0 -1 0 0 0 0 0 0 0 0 0 0 "
                  "1]"_mat};
    IntMatrix<> H = A;
    SquareMatrix<std::int64_t> U{SquareDims{H.numRow()}};
    NormalForm::hermite(H, U);
    utils::print("\n\n\n====\n\nH=");
    H.print();
    utils::print("\nU=");
    U.print();
    utils::print('\n');
    EXPECT_TRUE(isHNF(H));
    EXPECT_TRUE(H == U * A);
  }
  {
    IntMatrix<> A = "[-3 -1 1; 0 0 -2]"_mat;
    IntMatrix<> H = A;
    SquareMatrix<std::int64_t> U{SquareDims{H.numRow()}};
    NormalForm::hermite(H, U);
    EXPECT_TRUE(isHNF(H));
    EXPECT_TRUE(U * A == H);
    utils::print("A = \n");
    A.print();
    utils::print("\nH =\n");
    H.print();
    utils::print("\nU =\n");
    U.print();
    utils::print('\n');
  }
  {
    IntMatrix<> A =
      "[3 3 -3 1 0 -1 -2 1 1 2 -1; 3 3 -3 1 1 -3 2 0 3 0 -3; 2 -3 -2 -1 1 -2 3 3 3 3 -3]"_mat;
    IntMatrix<> H = A;
    SquareMatrix<std::int64_t> U{SquareDims{H.numRow()}};
    NormalForm::hermite(H, U);
    EXPECT_TRUE(isHNF(H));
    EXPECT_TRUE(U * A == H);
    utils::print("A = \n");
    A.print();
    utils::print("\nH =\n");
    H.print();
    utils::print("\nU =\n");
    U.print();
    utils::print('\n');
  }
}

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(NullSpaceTests, BasicAssertions) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(-10, 100);

  std::ptrdiff_t num_iters = 1;
  constexpr std::ptrdiff_t max_col = 10;
  std::int64_t mem[8 * (2 * 8z + 3 * max_col)];
  MutDensePtrMatrix<std::int64_t> ns_init{mem, DenseDims<>{row(8), col(8)}};
  for (std::ptrdiff_t num_col = 2; num_col <= max_col; num_col += 2) {
    MutDensePtrMatrix<std::int64_t> B{mem + (8z * 8),
                                      DenseDims<>{row(8), col(num_col)}},
      Bc{mem + (8z * (8 + max_col)), DenseDims<>{row(8), col(num_col)}},
      Y{mem + (8z * (8 + 2 * max_col)), DenseDims<>{row(8), col(num_col)}};
    std::ptrdiff_t null_dim = 0;
    for (std::ptrdiff_t i = 0; i < num_iters; ++i) {
      for (auto &&b : B) {
        b = distrib(gen);
        b = b > 10 ? 0 : b;
      }
      MutDensePtrMatrix<std::int64_t> NS =
        NormalForm::nullSpace(ns_init, Bc << B);
      null_dim += std::ptrdiff_t(NS.numRow());
      MutDensePtrMatrix<std::int64_t> Z{Y[_(NS.numRow()), _]};
      Z << NS * B;
      if (!allZero(Z)) {
        utils::print("B = \n");
        B.print();
        utils::print("\nNS = \n");
        NS.print();
        utils::print("\nZ = \n");
        Z.print();
        utils::print('\n');
      }
      for (auto &z : Z) EXPECT_EQ(z, 0);
      MutDensePtrMatrix<std::int64_t> ZN{
        mem + (8z * (8 + 3 * max_col)),
        DenseDims<>{math::row(null_dim), math::col(null_dim)}};
      EXPECT_EQ(NormalForm::nullSpace(ZN, NS).numRow(), 0);
    }
    utils::print("Average tested null dim = ");
    utils::print(double(null_dim) / double(num_iters));
    utils::print('\n');
  }
}

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(SimplifySystemTests, BasicAssertions) {
  IntMatrix<> A = "[2 4 5 5 -5; -4 3 -4 -3 -1; 1 0 -2 1 -4; -4 -2 3 -2 -1]"_mat,
              B = "[-6 86 -27 46 0 -15; -90 -81 91 44 -2 78; 4 -54 -98 "
                  "80 -10 82; -98 -15 -28 98 82 87]"_mat,
              At = A.t(), Bt = B.t();
  NormalForm::solveSystemRight(At, Bt);
  NormalForm::solveSystem(A, B);
  IntMatrix<> sA = "[-3975 0 0 0 -11370; 0 -1325 0 0 -1305; "
                   "0 0 -265 0 -347; 0 0 0 265 -1124]"_mat;
  IntMatrix<> true_b =
    "[-154140 -128775 -205035 317580 83820 299760; -4910 -21400 -60890 "
    "44820 14480 43390; -1334 -6865 -7666 8098 -538 9191; -6548 -9165 "
    "-24307 26176 4014 23332]"_mat;

  EXPECT_EQ(sA, A);
  EXPECT_EQ(true_b, B);
  EXPECT_EQ(sA, At.t());
  EXPECT_EQ(true_b, Bt.t());

  IntMatrix<> C = "[1 1 0; 0 1 1; 1 2 1]"_mat;
  IntMatrix<> D = "[1 0 0; 0 1 0; 0 0 1]"_mat;
  NormalForm::simplifySystem(C, D);
  IntMatrix<> true_c = "[1 0 -1; 0 1 1]"_mat;
  IntMatrix<> true_d = "[1 -1 0; 0 1 0]"_mat;
  EXPECT_EQ(true_c, C);
  EXPECT_EQ(true_d, D);
}

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(BareissTests, BasicAssertions) {
  IntMatrix<> A =
    "[-4 3 -2 2 -5; -5 1 -1 2 -5; -1 0 5 -3 2; -4 5 -4 -2 -4]"_mat;
  auto piv = NormalForm::bareiss(A);
  IntMatrix<> B =
    "[-4 3 -2 2 -5; 0 11 -6 2 -5; 0 0 56 -37 32; 0 0 0 -278 136]"_mat;
  EXPECT_EQ(A, B);
  Vector<std::ptrdiff_t> true_piv{std::array{0, 1, 2, 3}};
  EXPECT_EQ(piv, true_piv);

  IntMatrix<> C = "[-2 -2 -1 -2 -1; 1 1 2 2 -2; -2 2 2 -1 "
                  "-1; 0 0 -2 1 -1; -1 -2 2 1 -1]"_mat;
  IntMatrix<> D = "[-2 -2 -1 -2 -1; 0 -8 -6 -2 0; 0 0 -12 -8 "
                  "20; 0 0 0 -28 52; 0 0 0 0 -142]"_mat;
  auto pivots = NormalForm::bareiss(C);
  EXPECT_EQ(C, D);
  auto true_pivots = Vector<std::ptrdiff_t, 16>{"[0 2 2 3 4]"_mat};
  EXPECT_EQ(pivots, true_pivots);
}

// // NOLINTNEXTLINE(modernize-use-trailing-return-type)
// TEST(BareissSolveSystemTests, BasicAssertions) {
//   // Test with the same matrices as SimplifySystemTests to ensure
//   compatibility IntMatrix<> A = "[2 4 5 5 -5; -4 3 -4 -3 -1; 1 0 -2 1 -4; -4
//   -2 3 -2 -1]"_mat,
//               B = "[-6 86 -27 46 0 -15; -90 -81 91 44 -2 78; 4 -54 -98 "
//                   "80 -10 82; -98 -15 -28 98 82 87]"_mat,
//               At = A.t(), Bt = B.t();

//   // Store original matrices for comparison
//   IntMatrix<> A_orig = A, B_orig = B, At_orig = At, Bt_orig = Bt;

//   // Test bareissSolveSystemRight
//   NormalForm::bareissSolveSystemRight(At, Bt);
//   // Test bareissSolveSystem
//   NormalForm::bareissSolveSystem(A, B);

//   // Expected results (same as original solveSystem tests)
//   IntMatrix<> sA = "[-3975 0 0 0 -11370; 0 -1325 0 0 -1305; "
//                    "0 0 -265 0 -347; 0 0 0 265 -1124]"_mat;
//   IntMatrix<> true_b =
//     "[-154140 -128775 -205035 317580 83820 299760; -4910 -21400 -60890 "
//     "44820 14480 43390; -1334 -6865 -7666 8098 -538 9191; -6548 -9165 "
//     "-24307 26176 4014 23332]"_mat;
//   EXPECT_EQ(sA, A);
//   EXPECT_EQ(true_b, B);
//   EXPECT_EQ(sA, At.t());
//   EXPECT_EQ(true_b, Bt.t());

//   // Test smaller system
//   IntMatrix<> C = "[1 1 0; 0 1 1; 1 2 1]"_mat;
//   IntMatrix<> D = "[1 0 0; 0 1 0; 0 0 1]"_mat;
//   NormalForm::bareissSolveSystem(C, D);
//   IntMatrix<> true_c = "[1 0 -1; 0 1 1]"_mat;
//   IntMatrix<> true_d = "[1 -1 0; 0 1 0]"_mat;
//   EXPECT_EQ(true_c, C);
//   EXPECT_EQ(true_d, D);

//   // Test that bareiss and HNF versions produce equivalent results on fresh
//   // matrices
//   IntMatrix<> A1 = A_orig, B1 = B_orig;
//   IntMatrix<> A2 = A_orig, B2 = B_orig;

//   NormalForm::solveSystem(A1, B1);        // Original HNF version
//   NormalForm::bareissSolveSystem(A2, B2); // New Bareiss version

//   EXPECT_EQ(A1, A2); // Both should produce same A result
//   EXPECT_EQ(B1, B2); // Both should produce same B result

//   // Test right versions as well
//   IntMatrix<> At1 = At_orig, Bt1 = Bt_orig;
//   IntMatrix<> At2 = At_orig, Bt2 = Bt_orig;

//   NormalForm::solveSystemRight(At1, Bt1);        // Original HNF version
//   NormalForm::bareissSolveSystemRight(At2, Bt2); // New Bareiss version

//   EXPECT_EQ(At1, At2); // Both should produce same A result
//   EXPECT_EQ(Bt1, Bt2); // Both should produce same B result
// }

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(InvTest, BasicAssertions) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(-10, 10);
  alloc::OwningArena<> alloc;
  const std::ptrdiff_t num_iters = 1000;
  for (std::ptrdiff_t dim = 1; dim < 5; ++dim) {
    auto s0 = alloc.scope();
    MutSquarePtrMatrix<std::int64_t> B{
      square_matrix<std::int64_t>(&alloc, dim)};
    for (std::ptrdiff_t i = 0; i < num_iters; ++i) {
      while (true) {
        for (std::ptrdiff_t n = 0; n < dim * dim; ++n)
          B.data()[n] = distrib(gen);
        if (NormalForm::rank(alloc, B) == dim) break;
      }
      auto s1 = alloc.scope();
      // Da * B^{-1} = Binv0
      // Da = Binv0 * B
      MutSquarePtrMatrix<std::int64_t> Da{
        square_matrix<std::int64_t>(&alloc, dim)};
      Da << B;
      auto Binv0 = NormalForm::inv(&alloc, Da);
      MutSquarePtrMatrix<std::int64_t> Bc{
        square_matrix<std::int64_t>(&alloc, dim)};
      Bc << B;
      auto [Binv1, s] = NormalForm::scaledInv(&alloc, Bc);
      EXPECT_TRUE(Da.isDiagonal());
      EXPECT_EQ((Binv0 * B), Da);
      Da.diag() << s;
      if (B * Binv1 != Da) {
        utils::print("\nB = ");
        B.print();
        utils::print("\nDa = ");
        Da.print();
        utils::print("\nBinv0 = ");
        Binv0.print();
        utils::print("\nBinv1 = ");
        Binv1.print();
        utils::print("\ns = ");
        utils::print(s);
        utils::print('\n');
      }
      EXPECT_EQ(B * Binv1, Da);
    }
  }
}

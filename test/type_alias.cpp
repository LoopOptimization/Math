#include <gtest/gtest.h>
#ifndef USE_MODULE
#include "Math/Array.cxx"
#include "Math/AxisTypes.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/MatrixDimensions.cxx"
#include <iostream>
#else
import Array;
import ManagedArray;
import std;
#endif

using namespace math;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
TEST(TypeAliasChecks, BasicAssertions) {
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
    std::cout << "A=\n" << A43 << "\n";
    EXPECT_EQ((A43[0, 0]), 2);
    EXPECT_EQ((A43[1, 0]), 3);
    EXPECT_EQ((A43[2, 0]), 6);
    EXPECT_EQ((A43[3, 0]), 2);
    EXPECT_EQ((A43[0, 1]), 5);
    EXPECT_EQ((A43[1, 1]), 6);
    EXPECT_EQ((A43[2, 1]), 1);
    EXPECT_EQ((A43[3, 1]), 6);
    EXPECT_EQ((A43[0, 2]), 8);
    EXPECT_EQ((A43[1, 2]), 3);
    EXPECT_EQ((A43[2, 2]), 1);
    EXPECT_EQ((A43[3, 2]), 1);
  }
}

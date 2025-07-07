import boost.ut;
import Array;
import ArrayParse;
import CorePrint;
import ManagedArray;
import MatDim;
import StaticArray;
import std;

using namespace boost::ut;
using namespace ::math;
using utils::operator""_mat;

auto autoConvert(::math::PtrMatrix<std::int64_t> A) -> std::int64_t {
  std::int64_t s = 0;
  for (std::ptrdiff_t m = 0; m < A.numRow(); ++m)
    for (std::ptrdiff_t n = 0; n < A.numCol(); ++n) s += A[m, n];
  return s;
}

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
void testBasicAssertions() {
  IntMatrix<> A{"[0 3 -2 1; 3 -1 -2 -2; 2 0 -3 0]"_mat};
  utils::print("A = \n");
  A.print();
  utils::print('\n');
  expect((A[0, 0]) == 0);
  expect((A[0, 1]) == 3);
  expect((A[0, 2]) == -2);
  expect((A[0, 3]) == 1);
  expect((A[1, 0]) == 3);
  expect((A[1, 1]) == -1);
  expect((A[1, 2]) == -2);
  expect((A[1, 3]) == -2);
  expect((A[2, 0]) == 2);
  expect((A[2, 1]) == 0);
  expect((A[2, 2]) == -3);
  expect((A[2, 3]) == 0);
#ifndef POLYMATHNOEXPLICITSIMDARRAY
  static_assert(std::same_as<::math::StaticDims<std::int64_t, 2, 3, false>,
                             ::math::StridedDims<2, 3, 4>>);
#else
  static_assert(std::same_as<::math::StaticDims<std::int64_t, 2, 3, false>,
                             ::math::DenseDims<2, 3>>);
#endif
  expect(autoConvert("[1 2 3; 4 5 6]"_mat) == 21);
}
void testBasicAssertions2() {
  IntMatrix<> A = "[-1 0 1 0 0; 0 -1 1 0 0; 0 0 -1 1 0; 0 0 -1 0 1]"_mat;
  expect((A[0, 0]) == -1);
}

int main() {
  "StringParse BasicAssertions"_test = [] {
    testBasicAssertions();
  };
  "StringParse2 BasicAssertions"_test = [] {
    testBasicAssertions2();
  };
  return 0;
}

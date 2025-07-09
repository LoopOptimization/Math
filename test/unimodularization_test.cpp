import boost.ut;
import ArrayParse;
import CorePrint;
import ManagedArray;
import std;
import Unimodularization;

using namespace boost::ut;
using namespace ::math;
using utils::operator""_mat;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
void testBasicAssertions() {
  /*
  IntMatrix<> VE{"[0 1; 1 0; 0 1; 1 0]"_mat};
  utils::print("VE="); VE.print(); utils::print('\n');
  auto VB = unimodularize(VE);
  ASSERT_TRUE(VB.has_value());
  utils::print("VB="); VB->print(); utils::print('\n');

  IntMatrix<> A23{"[9 5; -5 -2; 1 0]"_mat};
  auto B = unimodularize(A23);
  ASSERT_TRUE(B.has_value());
  utils::print("B="); B->print(); utils::print('\n');
  // EXPECT_EQ(j, length(bsc));
  // EXPECT_EQ(j, length(bs));
*/
  IntMatrix<> A13{"[6; -5; 15]"_mat};
  utils::print("A13=");
  A13.print();
  utils::print('\n');
  expect((A13[0, 0]) == 6);
  expect((A13[1, 0]) == -5);
  expect((A13[2, 0]) == 15);
  auto test6_10_15 = unimodularize(A13); //, 1, 93, 1001);
  expect(test6_10_15.has_value());
  A13[0, 0] = 102;
  A13[1, 0] = 190;
  A13[2, 0] = 345;
  auto test102_190_345 = unimodularize(A13); //, 1, 93, 1001);
  expect(test102_190_345.has_value());
}

int main() {
  "UnimodularizationTest BasicAssertions"_test = [] { testBasicAssertions(); };
  return 0;
}

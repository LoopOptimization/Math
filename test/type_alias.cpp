import boost.ut;
import Array;
import CorePrint;
import ManagedArray;
import std;

using namespace boost::ut;
using namespace ::math;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
void testBasicAssertions() {
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
    expect((A43[0, 0]) == 2);
    expect((A43[1, 0]) == 3);
    expect((A43[2, 0]) == 6);
    expect((A43[3, 0]) == 2);
    expect((A43[0, 1]) == 5);
    expect((A43[1, 1]) == 6);
    expect((A43[2, 1]) == 1);
    expect((A43[3, 1]) == 6);
    expect((A43[0, 2]) == 8);
    expect((A43[1, 2]) == 3);
    expect((A43[2, 2]) == 1);
    expect((A43[3, 2]) == 1);
  }
}

int main() {
  "TypeAliasChecks BasicAssertions"_test = [] { testBasicAssertions(); };
  return 0;
}

import boost.ut;

import BaseUtils;
import std;

using namespace boost::ut;

void testBasicAssertions() {

  for (int i = 0; i++ < 32;) {
    double di = i;
    for (int j = 0; j++ < i;) {
      double dj = j;
      auto [x, y] = ::math::lower_bound_factor(di, dj);
      expect(di == x * y);
      expect(x <= dj);
      expect(std::round(x) == x);
      expect(std::round(y) == y);
    }
  }
}

int main() {
  "FactorLowerBound BasicAssertions"_test = [] { testBasicAssertions(); };
  return 0;
}

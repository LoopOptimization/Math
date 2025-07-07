import boost.ut;
import std;
import Tuple;

using namespace boost::ut;
using containers::Tuple, containers::tie, containers::Add, containers::Pair;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
void testBasicAssertions() {
  double x = 2.3;
  double y = 4.5;
  long z = -5;
  unsigned long w = 15;
  Tuple t{._0 = w, ._1 = x, ._2 = y, ._3 = z};
  {
    auto [a, b, c, d] = t;
    static_assert(std::same_as<decltype(a), unsigned long>);
    static_assert(std::same_as<decltype(x), double>);
    static_assert(std::same_as<decltype(y), double>);
    static_assert(std::same_as<decltype(z), long>);
    expect(a == w);
    expect(b == x);
    expect(c == y);
    expect(d == z);
  }
  Tuple t1 = t.map([](auto arg) { return 3 * arg; });
  {
    auto [a, b, c, d] = t1;
    static_assert(std::same_as<decltype(a), unsigned long>);
    static_assert(std::same_as<decltype(x), double>);
    static_assert(std::same_as<decltype(y), double>);
    static_assert(std::same_as<decltype(z), long>);
    expect(a == 3 * w);
    expect(b == 3 * x);
    expect(c == 3 * y);
    expect(d == 3 * z);
  }
  t1.apply([](auto &arg) { arg *= 2; });
  {
    auto [a, b, c, d] = t1;
    static_assert(std::same_as<decltype(a), unsigned long>);
    static_assert(std::same_as<decltype(x), double>);
    static_assert(std::same_as<decltype(y), double>);
    static_assert(std::same_as<decltype(z), long>);
    expect(a == 6 * w);
    expect(b == 6 * x);
    expect(c == 6 * y);
    expect(d == 6 * z);
    tie(w, x, y, z) = t1;
    expect(a == w);
    expect(b == x);
    expect(c == y);
    expect(d == z);
  }

  // auto tr = tie(Add(x), y);
  double xtest = x + 4.7;
  tie(Add(x), y) = Tuple(4.7, 5.0);
  expect(x == xtest);
  expect(y == 5.0);
  tie(Add(x), y) = Pair(4.7, 8.0);
  expect(x == xtest + 4.7);
  expect(y == 8.0);
  auto t16 = Tuple(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
  expect(t16.template get<0>() == 0);
  expect(t16.template get<1>() == 1);
  expect(t16.template get<2>() == 2);
  expect(t16.template get<3>() == 3);
  expect(t16.template get<4>() == 4);
  expect(t16.template get<5>() == 5);
  expect(t16.template get<6>() == 6);
  expect(t16.template get<7>() == 7);
  expect(t16.template get<8>() == 8);
  expect(t16.template get<9>() == 9);
  expect(t16.template get<10>() == 10);
  expect(t16.template get<11>() == 11);
  expect(t16.template get<12>() == 12);
  expect(t16.template get<13>() == 13);
  expect(t16.template get<14>() == 14);
  expect(t16.template get<15>() == 15);
}

int main() {
  "TupleTest BasicAssertions"_test = [] {
    testBasicAssertions();
  };
  return 0;
}

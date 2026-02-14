import Testing;
import Optional;
import std;

using namespace testing;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
void testBasicAssertions() {
  utils::Optional<std::ptrdiff_t> x{3}, y,
    z{std::numeric_limits<std::ptrdiff_t>::min()};
  expect(x.hasValue());
  expect(*x == 3);
  expect(!y.hasValue());
  expect(!z.hasValue());
  x = 14;
  expect(x.hasValue());
  expect(*x == 14);
  y = 33;
  expect(y.hasValue());
  expect(*y == 33);
  y = 0;
  expect(y.hasValue());
  expect(!*y);
  y = {};
  expect(!y.hasValue());

  std::ptrdiff_t a = 42, b = 11, c = 8;
  utils::Optional<std::ptrdiff_t *> p{&a}, q;
  static_assert(sizeof(utils::Optional<std::ptrdiff_t *>) ==
                sizeof(std::ptrdiff_t *));
  expect(p.hasValue());
  expect(!q.hasValue());
  **p += 10;
  expect(**p == 52);
  q = &b;
  expect(q.hasValue());
  p = &c;
  **p += 18;
  expect(a == 52);
  expect(**q == 11);
  expect(c == 26);
  p = {};
  expect(!p.hasValue());
}

int main() {
  "Optional BasicAssertions"_test = [] { testBasicAssertions(); };
  return 0;
}

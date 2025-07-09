import boost.ut;

using namespace boost::ut;
import GCD;
import SIMD;
import std;
import Tuple;

template <std::ptrdiff_t W>
constexpr auto vecget(simd::Vec<W, std::int64_t> v, std::ptrdiff_t i)
  -> std::int64_t {
  return v[i];
}
template <>
constexpr auto vecget<1>(std::int64_t v, std::ptrdiff_t) -> std::int64_t {
  return v;
}
template <std::ptrdiff_t W>
constexpr void vecset(simd::Vec<W, std::int64_t> &v, std::int64_t newval,
                      std::ptrdiff_t i) {
  v[i] = newval;
}
template <>
constexpr void vecset<1>(std::int64_t &v, std::int64_t newval, std::ptrdiff_t) {
  v = newval;
}

template <std::ptrdiff_t W>
inline auto fillGCD(std::mt19937 &gen)
  -> containers::Tuple<simd::Vec<W, std::int64_t>, simd::Vec<W, std::int64_t>,
                       std::array<std::int64_t, std::size_t(W)>, std::int64_t> {
  std::array<std::int64_t, std::size_t(W)> z;
  simd::Vec<W, std::int64_t> a, b;
  std::uniform_int_distribution<> distrib(std::numeric_limits<int>::min(),
                                          std::numeric_limits<int>::max());
  std::int64_t rg = 0;
  for (std::ptrdiff_t j = 0; j < W; ++j) {
    std::int64_t w = distrib(gen), x = distrib(gen), y = std::gcd(w, x);
    expect(y == ::math::gcd(w, x));
    z[j] = y;
    vecset<W>(a, w, j);
    vecset<W>(b, x, j);
    rg = std::gcd(rg, y);
  }
  return {a, b, z, rg};
}

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
int main() {
  "SIMDGCDTest BasicAssertions"_test = [] {
    constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 20000; ++i) {
      auto [a, b, z, rg] = fillGCD<W>(gen);
      simd::Vec<W, std::int64_t> g = ::math::gcd<W>(a, b);
      for (std::ptrdiff_t j = 0; j < W; ++j) expect(vecget<W>(g, j) == z[j]);
      expect(rg == ::math::gcdreduce<W>(g));
    }
  };
  return 0;
}

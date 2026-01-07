import boost.ut;

using namespace boost::ut;
import GCD;
import SIMD;
import std;
import Tuple;

namespace {
template <std::ptrdiff_t W>
constexpr auto vecget(simd::Vec<W, std::int64_t> v, std::ptrdiff_t i)
  -> std::int64_t {
  return v[i];
}
template <std::ptrdiff_t W>
constexpr void vecset(simd::Vec<W, std::int64_t> &v, std::int64_t newval,
                      std::ptrdiff_t i) {
  v[i] = newval;
}

template <std::ptrdiff_t W>
inline auto fillGCD(std::mt19937 &gen)
  -> containers::Tuple<simd::Vec<W, std::int64_t>, simd::Vec<W, std::int64_t>,
                       std::array<std::int64_t, std::size_t(W)>, std::int64_t> {
  std::array<std::int64_t, std::size_t(W)> z;
  simd::Vec<W, std::int64_t> a, b;
  std::uniform_int_distribution<> distrib(std::numeric_limits<int>::min() + 1,
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

template <typename T, std::ptrdiff_t W>
inline auto fillGCDUnroll(std::mt19937 &gen,
                          T lb = std::numeric_limits<T>::min() + 1,
                          T ub = std::numeric_limits<T>::max())
  -> containers::Tuple<std::array<T, std::size_t(W)>,
                       std::array<T, std::size_t(W)>,
                       std::array<T, std::size_t(W)>, T> {
  std::array<T, std::size_t(W)> a, b, z;
  std::uniform_int_distribution<T> distrib(lb, ub);
  T rg = 0;
  for (std::ptrdiff_t j = 0; j < W; ++j) {
    T w = T(distrib(gen)), x = T(distrib(gen));
    T y = std::gcd(w, x);
    expect(y == ::math::gcd(w, x));
    a[j] = w;
    b[j] = x;
    z[j] = y;
    rg = std::gcd(rg, y);
  }
  return {a, b, z, rg};
}

template <typename T, std::ptrdiff_t W>
inline auto makeSVec(const std::array<T, std::size_t(W)> &arr)
  -> simd::SVec<T, W> {
  simd::SVec<T, W> result = simd::SVec<T, W>::range(T(0));
  for (std::ptrdiff_t i = 0; i < W; ++i) {
    // Use extract_value to read, but we need to construct differently
    // since extract_value returns by value, not reference
  }
  // Use a different approach: construct from array data directly
  // by accessing the underlying structure
  static_assert(sizeof(result) >= sizeof(arr), "SVec too small");
  std::memcpy(&result, arr.data(), sizeof(arr));
  return result;
}
} // namespace

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
int main() {
  "SIMDGCDTest BasicAssertions"_test = [] -> void {
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

  "UnrollGCDTest Int32"_test = [] -> void {
    constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z, rg] = fillGCDUnroll<std::int32_t, W>(gen);
      auto a = makeSVec<std::int32_t, W>(aarr);
      auto b = makeSVec<std::int32_t, W>(barr);
      auto g = ::math::gcd(a, b);
      for (std::ptrdiff_t j = 0; j < W; ++j) expect(g.extract_value(j) == z[j]);
      expect(rg == ::math::gcdreduce(g));
    }
  };

  "UnrollGCDTest Int64"_test = [] -> void {
    constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z, rg] = fillGCDUnroll<std::int64_t, W>(gen);
      auto a = makeSVec<std::int64_t, W>(aarr);
      auto b = makeSVec<std::int64_t, W>(barr);
      auto g = ::math::gcd(a, b);
      for (std::ptrdiff_t j = 0; j < W; ++j) expect(g.extract_value(j) == z[j]);
      expect(rg == ::math::gcdreduce(g));
    }
  };

  // TODO: Re-enable float/double tests after fixing SVec construction
  "UnrollGCDTest Float"_test = [] -> void {
    static constexpr std::int32_t bound = 1
                                          << std::numeric_limits<float>::digits;
    static constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z, rg_unused] =
        fillGCDUnroll<std::int32_t, W>(gen, -bound, bound);
      // Convert to float for testing
      std::array<float, W> aarrf, barrf;
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        aarrf[j] = float(aarr[j]);
        barrf[j] = float(barr[j]);
      }
      auto af = makeSVec<float, W>(aarrf);
      auto bf = makeSVec<float, W>(barrf);
      auto g = ::math::gcd(af, bf);
      for (std::ptrdiff_t j = 0; j < W; ++j)
        expect(std::int32_t(g.extract_value(j)) == z[j]);
    }
  };

  "UnrollGCDTest Double"_test = [] -> void {
    static constexpr std::int64_t bound =
      1Z << std::numeric_limits<double>::digits;
    static constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z, rg_unused] =
        fillGCDUnroll<std::int64_t, W>(gen, -bound, bound);
      // Convert to double for testing
      std::array<double, W> aarrd, barrd;
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        aarrd[j] = double(aarr[j]);
        barrd[j] = double(barr[j]);
      }
      auto ad = makeSVec<double, W>(aarrd);
      auto bd = makeSVec<double, W>(barrd);
      auto g = ::math::gcd(ad, bd);
      for (std::ptrdiff_t j = 0; j < W; ++j)
        expect(std::int64_t(g.extract_value(j)) == z[j]);
    }
  };

  "UnrollGCDTest EdgeCases"_test = [] -> void {
    // Test with zero values
    auto a =
      makeSVec<std::int32_t, 4>(std::array<std::int32_t, 4>{0, 12, -15, 0});
    auto b =
      makeSVec<std::int32_t, 4>(std::array<std::int32_t, 4>{18, 0, 25, 0});
    auto g = ::math::gcd(a, b);
    expect(g.extract_value(0) == 18);
    expect(g.extract_value(1) == 12);
    expect(g.extract_value(2) == 5);
    expect(g.extract_value(3) == 0);

    // Test negative values
    auto c =
      makeSVec<std::int64_t, 4>(std::array<std::int64_t, 4>{-12, -18, 24, -30});
    auto d =
      makeSVec<std::int64_t, 4>(std::array<std::int64_t, 4>{18, -24, -36, 45});
    auto h = ::math::gcd(c, d);
    expect(h.extract_value(0) == 6);
    expect(h.extract_value(1) == 6);
    expect(h.extract_value(2) == 12);
    expect(h.extract_value(3) == 15);
  };

  "UnrollLCMTest Int32"_test = [] -> void {
    constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z_unused, rg_unused] =
        fillGCDUnroll<std::int32_t, W>(gen, -32768, 32768);
      auto a = makeSVec<std::int32_t, W>(aarr);
      auto b = makeSVec<std::int32_t, W>(barr);
      auto lcm_result = ::math::lcm(a, b);
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        std::int32_t expected = std::lcm(aarr[j], barr[j]);
        expect(fatal(eq(lcm_result.extract_value(j), expected)));
      }
    }
  };

  "UnrollLCMTest Int64"_test = [] -> void {
    constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z_unused, rg_unused] =
        fillGCDUnroll<std::int64_t, W>(gen, -65535, 65535);
      auto a = makeSVec<std::int64_t, W>(aarr);
      auto b = makeSVec<std::int64_t, W>(barr);
      auto lcm_result = ::math::lcm(a, b);
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        std::int64_t expected = std::lcm(aarr[j], barr[j]);
        expect(fatal(eq(lcm_result.extract_value(j), expected)));
      }
    }
  };

  "UnrollLCMTest Float"_test = [] -> void {
    static constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z_unused, rg_unused] =
        fillGCDUnroll<std::int32_t, W>(gen, -4096, 4096);
      // Convert to float for testing
      std::array<float, W> aarrf, barrf;
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        aarrf[j] = float(aarr[j]);
        barrf[j] = float(barr[j]);
      }
      auto af = makeSVec<float, W>(aarrf);
      auto bf = makeSVec<float, W>(barrf);
      auto lcm_result = ::math::lcm(af, bf);
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        std::int32_t expected = std::lcm(aarr[j], barr[j]);
        expect(fatal(eq(std::int32_t(lcm_result.extract_value(j)), expected)));
      }
    }
  };

  "UnrollLCMTest Double"_test = [] -> void {
    static constexpr std::ptrdiff_t W = 16;
    std::random_device rd;
    std::mt19937 gen(rd());
    for (std::ptrdiff_t i = 0; i < 5000; ++i) {
      auto [aarr, barr, z_unused, rg_unused] =
        fillGCDUnroll<std::int64_t, W>(gen, -65536, 65536);
      // Convert to double for testing
      std::array<double, W> aarrd, barrd;
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        aarrd[j] = double(aarr[j]);
        barrd[j] = double(barr[j]);
      }
      auto ad = makeSVec<double, W>(aarrd);
      auto bd = makeSVec<double, W>(barrd);
      auto lcm_result = ::math::lcm(ad, bd);
      for (std::ptrdiff_t j = 0; j < W; ++j) {
        std::int64_t expected = std::lcm(aarr[j], barr[j]);
        expect(fatal(eq(std::int64_t(lcm_result.extract_value(j)), expected)));
      }
    }
  };

  "UnrollLCMTest EdgeCases"_test = [] -> void {
    // Test with ones and same values
    auto a =
      makeSVec<std::int32_t, 4>(std::array<std::int32_t, 4>{1, 12, 15, 7});
    auto b =
      makeSVec<std::int32_t, 4>(std::array<std::int32_t, 4>{18, 1, 15, 21});
    auto l = ::math::lcm(a, b);
    expect(l.extract_value(0) == 18); // lcm(1, 18) = 18
    expect(l.extract_value(1) == 12); // lcm(12, 1) = 12
    expect(l.extract_value(2) == 15); // lcm(15, 15) = 15
    expect(l.extract_value(3) == 21); // lcm(7, 21) = 21

    // Test negative values
    auto c =
      makeSVec<std::int64_t, 4>(std::array<std::int64_t, 4>{-12, -18, 24, -30});
    auto d =
      makeSVec<std::int64_t, 4>(std::array<std::int64_t, 4>{18, -24, -36, 45});
    auto m = ::math::lcm(c, d);
    expect(m.extract_value(0) == 36); // lcm(-12, 18) = 36
    expect(m.extract_value(1) == 72); // lcm(-18, -24) = 72
    expect(m.extract_value(2) == 72); // lcm(24, -36) = 72
    expect(m.extract_value(3) == 90); // lcm(-30, 45) = 90

    // Test with small coprime numbers
    auto e =
      makeSVec<std::int32_t, 4>(std::array<std::int32_t, 4>{3, 5, 7, 11});
    auto f =
      makeSVec<std::int32_t, 4>(std::array<std::int32_t, 4>{5, 7, 11, 13});
    auto n = ::math::lcm(e, f);
    expect(n.extract_value(0) == 15);  // lcm(3, 5) = 15
    expect(n.extract_value(1) == 35);  // lcm(5, 7) = 35
    expect(n.extract_value(2) == 77);  // lcm(7, 11) = 77
    expect(n.extract_value(3) == 143); // lcm(11, 13) = 143
  };

  return 0;
}

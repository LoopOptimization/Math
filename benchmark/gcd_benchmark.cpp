
import Nanobench;
import GCD;
import SIMD;
import std;

namespace {

constexpr std::ptrdiff_t W = 16;
constexpr std::size_t X = 32768 / (2 * sizeof(simd::SVec<std::int64_t, W>));
// Bound for exact integer representation in double
static constexpr std::int64_t bound = 1Z << std::numeric_limits<double>::digits;

// Shared test data - same values used for int64 and double benchmarks
std::array<simd::SVec<std::int64_t, W>, X> a_vecs, b_vecs;

// Native SIMD width for Vec benchmarks (dgcdx)
constexpr std::ptrdiff_t VW = simd::Width<std::int64_t>;
constexpr std::size_t VX = 32768 / (2 * sizeof(simd::Vec<VW, std::int64_t>));
std::array<simd::Vec<VW, std::int64_t>, VX> a_native_vecs, b_native_vecs;

// Double Vec arrays for dgcdx double benchmark
constexpr std::ptrdiff_t VWD = simd::Width<double>;
constexpr std::size_t VXD = 32768 / (2 * sizeof(simd::Vec<VWD, double>));
std::array<simd::Vec<VWD, double>, VXD> a_native_dvecs, b_native_dvecs;

} // namespace

void init_gcd_data() {
  std::mt19937_64 rng{42};
  std::uniform_int_distribution<std::int64_t> dist(-bound, bound);
  for (std::size_t i = 0; i < X; ++i) {
    std::array<std::int64_t, std::size_t(W)> a, b;
    for (std::ptrdiff_t j = 0; j < W; ++j) {
      a[j] = dist(rng);
      b[j] = dist(rng);
    }
    std::memcpy(&a_vecs[i], a.data(), sizeof(a));
    std::memcpy(&b_vecs[i], b.data(), sizeof(b));
  }
  // Initialize native Vec arrays for dgcdx benchmarks
  for (std::size_t i = 0; i < VX; ++i) {
    std::array<std::int64_t, std::size_t(VW)> a, b;
    for (std::ptrdiff_t j = 0; j < VW; ++j) {
      a[j] = dist(rng);
      b[j] = dist(rng);
    }
    std::memcpy(&a_native_vecs[i], a.data(), sizeof(a));
    std::memcpy(&b_native_vecs[i], b.data(), sizeof(b));
  }
  // Initialize double Vec arrays for dgcdx double benchmark
  for (std::size_t i = 0; i < VXD; ++i) {
    std::array<double, std::size_t(VWD)> a, b;
    for (std::ptrdiff_t j = 0; j < VWD; ++j) {
      a[j] = static_cast<double>(dist(rng));
      b[j] = static_cast<double>(dist(rng));
    }
    std::memcpy(&a_native_dvecs[i], a.data(), sizeof(a));
    std::memcpy(&b_native_dvecs[i], b.data(), sizeof(b));
  }
}

void BM_gcd_scalar_int64(Bench &bench) {
  std::mt19937_64 rng{42};
  std::uniform_int_distribution<std::int64_t> dist(-bound, bound);
  constexpr std::size_t N = 32768 / (2 * sizeof(std::int64_t));
  std::array<std::int64_t, N> a_scalars, b_scalars;
  for (std::size_t i = 0; i < N; ++i) {
    a_scalars[i] = dist(rng);
    b_scalars[i] = dist(rng);
  }

  bench.run("std::gcd<int64_t>", [&] {
    for (std::size_t i = 0; i < N; ++i) {
      auto g = std::gcd(a_scalars[i], b_scalars[i]);
      doNotOptimizeAway(g);
    }
  });

  bench.run("math::gcd<int64_t>", [&] {
    for (std::size_t i = 0; i < N; ++i) {
      auto g = ::math::gcd(a_scalars[i], b_scalars[i]);
      doNotOptimizeAway(g);
    }
  });
}

void BM_gcd_svec_int64(Bench &bench) {
  bench.run("math::gcd<SVec<int64_t,16>>", [&] {
    for (std::size_t i = 0; i < X; ++i) {
      auto g = ::math::gcd(a_vecs[i], b_vecs[i]);
      doNotOptimizeAway(g);
    }
  });
}

void BM_gcd_svec_double(Bench &bench) {
  bench.run("math::gcd<SVec<double,16>>", [&] {
    for (std::size_t i = 0; i < X; ++i) {
      // Cast int64 SVec to double SVec - same values for fair comparison
      auto a = static_cast<simd::SVec<double, W>>(a_vecs[i]);
      auto b = static_cast<simd::SVec<double, W>>(b_vecs[i]);
      auto g = ::math::gcd(a, b);
      doNotOptimizeAway(g);
    }
  });
}

void BM_gcdreduce_svec_int64(Bench &bench) {
  bench.run("math::gcdreduce<SVec<int64_t,16>>", [&] {
    for (std::size_t i = 0; i < X; ++i) {
      auto g = ::math::gcd(a_vecs[i], b_vecs[i]);
      auto r = ::math::gcdreduce(g);
      doNotOptimizeAway(r);
    }
  });
}

void BM_dgcdx_scalar_int64(Bench &bench) {
  std::mt19937_64 rng{42};
  std::uniform_int_distribution<std::int64_t> dist(-bound, bound);
  constexpr std::size_t N = 32768 / (2 * sizeof(std::int64_t));
  std::array<std::int64_t, N> a_scalars, b_scalars;
  for (std::size_t i = 0; i < N; ++i) {
    a_scalars[i] = dist(rng);
    b_scalars[i] = dist(rng);
  }

  bench.run("math::dgcdx<int64_t>", [&] {
    for (std::size_t i = 0; i < N; ++i) {
      auto result = ::math::dgcdx(a_scalars[i], b_scalars[i]);
      doNotOptimizeAway(result);
    }
  });
}

void BM_dgcdx_vec_int64(Bench &bench) {
  bench.run("math::dgcdx<Vec<int64_t>>", [&] {
    for (std::size_t i = 0; i < VX; ++i) {
      auto result = ::math::dgcdx(a_native_vecs[i], b_native_vecs[i]);
      doNotOptimizeAway(result);
    }
  });
}

void BM_dgcdx_vec_double(Bench &bench) {
  bench.run("math::dgcdx<Vec<double>>", [&] {
    for (std::size_t i = 0; i < VXD; ++i) {
      auto result = ::math::dgcdx(a_native_dvecs[i], b_native_dvecs[i]);
      doNotOptimizeAway(result);
    }
  });
}


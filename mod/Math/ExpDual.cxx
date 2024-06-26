#pragma once
#include "Math/Dual.hpp"
#include "Math/Exp.hpp"
#include <cstddef>
#include <limits>
namespace poly::math {

template <int l = 8> constexpr auto smax(auto x, auto y, auto z) {
  double m = std::max(std::max(value(x), value(y)), value(z));
  static constexpr double f = l, i = 1 / f;
  return m + i * log(exp(f * (x - m)) + exp(f * (y - m)) + exp(f * (z - m)));
}
template <int l = 8> constexpr auto smax(auto w, auto x, auto y, auto z) {
  double m =
    std::max(std::max(value(w), value(y)), std::max(value(x), value(z)));
  static constexpr double f = l, i = 1 / f;
  return m + i * log(exp(f * (w - m)) + exp(f * (x - m)) + exp(f * (y - m)) +
                     exp(f * (z - m)));
}
template <int l = 8, typename T, ptrdiff_t N>
constexpr auto smax(SVector<T, N> x) -> T {
  static_assert(!std::is_integral_v<T>);
  static constexpr double f = l, i = 1 / f;
  double m = -std::numeric_limits<double>::infinity();
  for (ptrdiff_t n = 0; n < N; ++n) m = std::max(m, value(x[n]));
  T a{};
  for (ptrdiff_t n = 0; n < N; ++n) a += exp(f * (x[n] - m));
  return m + i * log(a);
}

} // namespace poly::math

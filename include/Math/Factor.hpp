#pragma once

#include <array>
#include <cmath>
namespace poly::math {

constexpr auto lower_bound_factor(double N, double x) -> std::array<double, 2> {
  double y = std::ceil(N / x);
  while (x * y != N) {
    x = std::floor(N / y);
    y = std::ceil(N / x);
  }
  return {x, y};
}

} // namespace poly::math

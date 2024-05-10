#pragma once
#include "Math/Exp.hpp"
#include "SIMD/Intrin.hpp"
#include "SIMD/Unroll.hpp"
#include <cstddef>

namespace poly::math {
template <ptrdiff_t W>
[[gnu::always_inline]] constexpr auto
expm1b_kernel(std::integral_constant<int, 2>,
              simd::Vec<W, double> x) -> simd::Vec<W, double> {
  return x * (((0.009618130135925114 * x + 0.055504115022757844) * x +
               0.2402265069590989) *
              x * 0.6931471805599393);
}
template <ptrdiff_t W>
[[gnu::always_inline]] constexpr auto
expm1b_kernel(std::integral_constant<int, 3>,
              simd::Vec<W, double> x) -> simd::Vec<W, double> {
  return x * (((0.04166666762124105 * x + 0.1666666704849642) * x +
               0.49999999999999983) *
                x +
              0.9999999999999998);
}
template <ptrdiff_t W>
[[gnu::always_inline]] constexpr auto
expm1b_kernel(std::integral_constant<int, 10>,
              simd::Vec<W, double> x) -> simd::Vec<W, double> {
  return x * ((((0.5393833837413015 * x + 1.1712561359457612) * x +
                2.0346785922926713) *
                 x +
               2.6509490552382577) *
                x +
              2.302585092994046);
}
template <int B, ptrdiff_t W>
constexpr auto exp_impl(simd::Vec<W, double> x) -> simd::Vec<W, double> {
  constexpr std::integral_constant<int, B> base{};
  // #if __FAST_MATH__
  using V = simd::Vec<W, double>;
  auto maxmask = simd::cmp::ge<W, double>(
    x, simd::vbroadcast<W, double>(max_exp(0.0, base)));
  V maxval = simd::vbroadcast<W, double>(std::numeric_limits<double>::max());
  auto minmask = simd::cmp::le<W, double>(
    x, simd::vbroadcast<W, double>(subnormal_exp(0.0, base)));
  V alt = simd::select<double>(maxmask, maxval, V{});
  auto altmask = maxmask | minmask;
  // #else
  // if (x >= max_exp(x, base)) return std::numeric_limits<double>::infinity();
  // #endif
  V floatN = x * simd::vbroadcast<W, double>(LogBo256INV(base)) +
             simd::vbroadcast<W, double>(magic_round_const<double>);
  auto N = std::bit_cast<simd::Vec<W, int64_t>>(floatN);
  floatN -= simd::vbroadcast<W, double>(magic_round_const<double>);

  V r = floatN * simd::vbroadcast<W, double>(LogBo256U(base)) + x;
  r = floatN * simd::vbroadcast<W, double>(LogBo256L(base)) + r;
  V jU = simd::gather(J_TABLE, simd::mask::None<W>{},
                      N & simd::vbroadcast<W, int64_t>(255));
  V small = jU * expm1b_kernel<W>(base, r) + jU;
  simd::Vec<W, int64_t> twopk = std::bit_cast<simd::Vec<W, int64_t>>(
    (std::bit_cast<simd::Vec<W, int64_t>>(N) >> 8) << 52);
  V z = std::bit_cast<V>(twopk + std::bit_cast<simd::Vec<W, int64_t>>(small));
  return simd::select<double>(altmask, alt, z);
}

template <ptrdiff_t W>
constexpr auto exp(simd::Vec<W, double> x) -> simd::Vec<W, double> {
  return exp_impl<3, W>(x);
}
template <ptrdiff_t W>
constexpr auto exp2(simd::Vec<W, double> x) -> simd::Vec<W, double> {
  return exp_impl<2, W>(x);
}
template <ptrdiff_t W>
constexpr auto exp10(simd::Vec<W, double> x) -> simd::Vec<W, double> {
  return exp_impl<10, W>(x);
}
template <int B, ptrdiff_t R, ptrdiff_t U, ptrdiff_t W>
constexpr auto
exp_impl(simd::Unroll<R, U, W, double> x) -> simd::Unroll<R, U, W, double> {
  if constexpr (R * U == 1) {
    return {exp_impl<B, W>(x.vec)};
  } else {
    simd::Unroll<R, U, W, double> ret;
    for (ptrdiff_t i = 0; i < R * U; ++i)
      ret.data[i] = exp_impl<B, W>(x.data[i]);
    return ret;
  }
}
template <ptrdiff_t R, ptrdiff_t U, ptrdiff_t W>
constexpr auto
exp(simd::Unroll<R, U, W, double> x) -> simd::Unroll<R, U, W, double> {
  return exp_impl<3>(x);
}
template <ptrdiff_t R, ptrdiff_t U, ptrdiff_t W>
constexpr auto
exp2(simd::Unroll<R, U, W, double> x) -> simd::Unroll<R, U, W, double> {
  return exp_impl<2>(x);
}
template <ptrdiff_t R, ptrdiff_t U, ptrdiff_t W>
constexpr auto
exp10(simd::Unroll<R, U, W, double> x) -> simd::Unroll<R, U, W, double> {
  return exp_impl<10>(x);
}

template <ptrdiff_t W>
constexpr auto sigmoid(simd::Vec<W, double> x) -> simd::Vec<W, double> {
  return 1.0 / (1.0 + exp<W>(-x));
}
} // namespace poly::math

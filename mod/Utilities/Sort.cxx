module Sort;

import SIMD;
import std;

namespace {

template <typename T> using V8 = simd::SVec<T, 8>;
template <typename T> using V16 = simd::SVec<T, 16>;
template <typename T> using V32 = simd::SVec<T, 32>;

template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] constexpr void minmax(simd::Unroll<1, N, W, T> &x,
                                             simd::Unroll<1, N, W, T> &y) {
  simd::Unroll<1, N, W, T> z = x;
  x = min(x, y);
  y = max(z, y);
}

template <std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse2x(simd::Vec<W, T> x) -> simd::Vec<W, T> {
  static_assert(W >= 2 && W <= 16);
  if constexpr (W == 2) {
    return __builtin_shufflevector(x, x, 1, 0);
  } else if constexpr (W == 4) {
    return __builtin_shufflevector(x, x, 1, 0, 3, 2);
  } else if constexpr (W == 8) {
    return __builtin_shufflevector(x, x, 1, 0, 3, 2, 5, 4, 7, 6);
  } else {
    return __builtin_shufflevector(x, x, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10,
                                   13, 12, 15, 14);
  }
}
template <std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse4x(simd::Vec<W, T> x) -> simd::Vec<W, T> {
  static_assert(W >= 4 && W <= 16);
  if constexpr (W == 4) {
    return __builtin_shufflevector(x, x, 3, 2, 1, 0);
  } else if constexpr (W == 8) {
    return __builtin_shufflevector(x, x, 3, 2, 1, 0, 7, 6, 5, 4);
  } else {
    return __builtin_shufflevector(x, x, 3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8,
                                   15, 14, 13, 12);
  }
}
template <std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse8x(simd::Vec<W, T> x) -> simd::Vec<W, T> {
  static_assert(W == 8 || W == 16);
  if constexpr (W == 8) {
    return __builtin_shufflevector(x, x, 7, 6, 5, 4, 3, 2, 1, 0);
  } else {
    return __builtin_shufflevector(x, x, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12,
                                   11, 10, 9, 8);
  }
}

template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse2(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(W >= 2);
  if constexpr (N != 1) {
    simd::Unroll<1, N, W, T> ret;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t i = 0; i < N; ++i) ret.data_[i] = reverse2x(x.data_[i]);
    return ret;
  } else return {reverse2x(x.vec_)};
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse4(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W >= 4 && W >= 2);
  if constexpr (N != 1) {
    simd::Unroll<1, N, W, T> ret;
    if constexpr (W == 2) {
#pragma clang loop unroll(full)
      for (std::ptrdiff_t i = 0; i < (N >> 1); ++i) {
        ret.data_[2 * i] = reverse2x(x.data_[2 * i + 1]);
        ret.data_[2 * i + 1] = reverse2x(x.data_[2 * i]);
      }
    } else
#pragma clang loop unroll(full)
      for (std::ptrdiff_t i = 0; i < N; ++i)
        ret.data_[i] = reverse4x(x.data_[i]);
    return ret;
  } else return {reverse4x(x.vec_)};
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse8(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W >= 8 && W >= 2);
  if constexpr (N != 1) {
    simd::Unroll<1, N, W, T> ret;
    if constexpr (W == 2) {
#pragma clang loop unroll(full)
      for (std::ptrdiff_t i = 0; i < (N >> 2); ++i) {
        ret.data_[4 * i] = reverse2x(x.data_[4 * i + 3]);
        ret.data_[4 * i + 1] = reverse2x(x.data_[4 * i + 2]);
        ret.data_[4 * i + 2] = reverse2x(x.data_[4 * i + 1]);
        ret.data_[4 * i + 3] = reverse2x(x.data_[4 * i]);
      }
    } else if constexpr (W == 4) {
#pragma clang loop unroll(full)
      for (std::ptrdiff_t i = 0; i < (N >> 1); ++i) {
        ret.data_[2 * i] = reverse4x(x.data_[2 * i + 1]);
        ret.data_[2 * i + 1] = reverse4x(x.data_[2 * i]);
      }
    } else
#pragma clang loop unroll(full)
      for (std::ptrdiff_t i = 0; i < N; ++i)
        ret.data_[i] = reverse8x(x.data_[i]);
    return ret;
  } else return {reverse8x(x.vec_)};
}
template <int _0, int _1, int _2, int _3, std::ptrdiff_t N, std::ptrdiff_t W,
          typename T>
[[gnu::always_inline]] auto vshuf(simd::Unroll<1, N, W, T> x,
                                  simd::Unroll<1, N, W, T> y)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8 || N * W == 16);
  if constexpr (W == 16) {
    return {__builtin_shufflevector(x.vec_, y.vec_, _0, _1, 16 + _2, 16 + _3,
                                    4 + _0, 4 + _1, 20 + _2, 20 + _3, 8 + _0,
                                    8 + _1, 24 + _2, 24 + _3, 12 + _0, 12 + _1,
                                    28 + _2, 28 + _3)};
  } else if constexpr (W == 8) {
    if constexpr (N == 1)
      return {__builtin_shufflevector(x.vec_, y.vec_, _0, _1, 8 + _2, 8 + _3,
                                      4 + _0, 4 + _1, 12 + _2, 12 + _3)};
    else
      return {__builtin_shufflevector(x.data_[0], y.data_[0], _0, _1, 8 + _2,
                                      8 + _3, 4 + _0, 4 + _1, 12 + _2, 12 + _3),
              __builtin_shufflevector(x.data_[1], y.data_[1], _0, _1, 8 + _2,
                                      8 + _3, 4 + _0, 4 + _1, 12 + _2,
                                      12 + _3)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t i = 0; i < N; ++i)
      ret.data_[i] =
        __builtin_shufflevector(x.data_[i], y.data_[i], _0, _1, 4 + _2, 4 + _3);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t i = 0; i < (N >> 1); ++i) {
      ret.data_[2 * i] =
        __builtin_shufflevector(x.data_[2 * i], x.data_[2 * i + 1], _0, _1);
      ret.data_[2 * i + 1] =
        __builtin_shufflevector(y.data_[2 * i], y.data_[2 * i + 1], _2, _3);
    }
    return ret;
  }
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto interleave_halves(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8 || N * W == 16);
  if constexpr (W == 16) {
    return {__builtin_shufflevector(x.vec_, x.vec_, 0, 4, 1, 5, 2, 6, 3, 7, 8,
                                    12, 9, 13, 10, 14, 11, 15)};
  } else if constexpr (W == 8) {
    if constexpr (N == 1) {
      return {__builtin_shufflevector(x.vec_, x.vec_, 0, 4, 1, 5, 2, 6, 3, 7)};
    } else {
      ::simd::Vec<W, T> a = x.data_[0], b = x.data_[1];
      return {__builtin_shufflevector(a, a, 0, 4, 1, 5, 2, 6, 3, 7),
              __builtin_shufflevector(b, b, 0, 4, 1, 5, 2, 6, 3, 7)};
    }
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t i = 0; i < (N >> 1); ++i) {
      simd::Vec<W, T> a = x.data_[2 * i], b = x.data_[2 * i + 1];
      ret.data_[2 * i] = __builtin_shufflevector(a, b, 0, 4, 1, 5);
      ret.data_[2 * i + 1] = __builtin_shufflevector(a, b, 2, 6, 3, 7);
    }
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t i = 0; i < (N >> 2); ++i) {
      simd::Vec<W, T> a = x.data_[4 * i], b = x.data_[4 * i + 1],
                      c = x.data_[4 * i + 2], d = x.data_[4 * i + 3];
      ret.data_[4 * i] = __builtin_shufflevector(a, c, 0, 2);
      ret.data_[4 * i + 1] = __builtin_shufflevector(a, c, 1, 3);
      ret.data_[4 * i + 2] = __builtin_shufflevector(b, d, 0, 2);
      ret.data_[4 * i + 3] = __builtin_shufflevector(b, d, 1, 3);
    }
    return ret;
  }
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto interleave(simd::Unroll<1, N, W, T> x,
                                       simd::Unroll<1, N, W, T> y)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8 || N * W == 16);
  if constexpr (W == 16) {
    return {__builtin_shufflevector(x.vec_, y.vec_, 0, 17, 2, 19, 4, 21, 6, 23,
                                    8, 25, 10, 27, 12, 29, 14, 31)};
  } else if constexpr (W == 8) {
    if constexpr (N == 1)
      return {
        __builtin_shufflevector(x.vec_, y.vec_, 0, 9, 2, 11, 4, 13, 6, 15)};
    else
      return {__builtin_shufflevector(x.data_[0], y.data_[0], 0, 9, 2, 11, 4,
                                      13, 6, 15),
              __builtin_shufflevector(x.data_[1], y.data_[1], 0, 9, 2, 11, 4,
                                      13, 6, 15)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t i = 0; i < N; ++i)
      ret.data_[i] =
        __builtin_shufflevector(x.data_[i], y.data_[i], 0, 5, 2, 7);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
#pragma clang loop unroll(full)
    for (std::ptrdiff_t i = 0; i < N; ++i)
      ret.data_[i] = __builtin_shufflevector(x.data_[i], y.data_[i], 0, 3);
    return ret;
  }
}
// Treats inputs as sets of `8`, sorting each set of `8` independently.
// If `xlo` and `xhi` are of len 8, it sorts them.
// [sorted...]
// If they're of length 16, it returns two sorted sets of 8.
// [ lower8_set0..., lower8_set1..., upper8_set0..., upper8_set1... ]
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
auto sort16(simd::Unroll<1, N, W, T> a, simd::Unroll<1, N, W, T> b) {
  using V = simd::Unroll<1, N, W, T>;
  minmax(a, b);

  b = reverse2(b);
  minmax(a, b);

  V c = a;
  a = vshuf<0, 2, 0, 2>(a, b);
  b = vshuf<1, 3, 1, 3>(c, b);
  minmax(a, b);

  b = reverse4(b);
  minmax(a, b);

  c = a;
  a = vshuf<0, 1, 0, 1>(a, b);
  b = vshuf<2, 3, 2, 3>(c, b);
  minmax(a, b);

  c = a;
  a = vshuf<0, 2, 1, 3>(a, b);
  b = vshuf<1, 3, 0, 2>(c, b);
  minmax(a, b);

  b = reverse8(b);
  minmax(a, b);

  c = a;
  a = vshuf<0, 2, 1, 3>(a, b);
  b = vshuf<1, 3, 0, 2>(c, b);
  minmax(a, b);

  c = a;
  a = vshuf<0, 2, 1, 3>(a, b);
  b = vshuf<1, 3, 0, 2>(c, b);
  minmax(a, b);

  a = interleave_halves(a);
  b = interleave_halves(b);

  c = a;
  a = vshuf<0, 2, 0, 2>(a, b);
  b = vshuf<1, 3, 1, 3>(c, b);
  minmax(a, b);

  V b2 = reverse2(b);
  V b1 = reverse2(a);

  a = interleave(a, b2);
  b = interleave(b1, b);

  return a.cat(b);
}

template <typename T> auto sort32(V32<T> x) -> V32<T> {
  auto [xlo, xhi] = x.split();

  // Sort each half of 16 elements
  auto [l01, h01] = sort16(xlo, xhi).split();
  auto [l0, l1] = l01.split();
  auto [h0, h1] = h01.split();
  // V16<T> lo_sorted = sort16(xlo);
  // V16<T> hi_sorted = sort16(xhi);

  // Bitonic merge: reverse the second half and merge
  V16<T> a = l0.cat(h0);
  V16<T> b = l1.cat(h1).reverse();

  // After compare-exchange at distance 16, all elements in a <= all in b
  minmax(a, b);

  // Now sort each 16-element half independently
  // a.split()
  // b.split()
  auto [a0, a1] = a.split();
  auto [b0, b1] = b.split();
  V32<T> ab = sort16(a0.cat(b0), a1.cat(b1));
  // a = sort16(a);
  // b = sort16(b);
  auto [ab0, ab1] = ab.split();
  auto [x0, y0] = ab0.split();
  auto [x1, y1] = ab1.split();
  return x0.cat(x1).cat(y0.cat(y1));

  // return a.cat(b);
}
} // namespace

namespace utils {

auto sort(math::SVector<double, 16> x) -> math::SVector<double, 16> {
  auto [xlo, xhi] = x.simd().split();
  return {sort16(xlo, xhi)};
}
auto sort(math::SVector<float, 16> x) -> math::SVector<float, 16> {
  auto [xlo, xhi] = x.simd().split();
  return {sort16(xlo, xhi)};
}

auto sort(math::SVector<double, 32> x) -> math::SVector<double, 32> {
  return {sort32(x.simd())};
}
auto sort(math::SVector<float, 32> x) -> math::SVector<float, 32> {
  return {sort32(x.simd())};
}

} // namespace utils

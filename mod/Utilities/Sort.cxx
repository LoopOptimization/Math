module Sort;

import SIMD;
import std;

namespace {

template <typename T> using V8 = simd::SVec<T, 8>;

template <typename T>
[[gnu::always_inline]] constexpr void minmax(V8<T> &x, V8<T> &y) {
  V8<T> z = x;
  x = min(x, y);
  y = max(z, y);
}

template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto swap_el_pairs(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8);
  if constexpr (W == 8) {
    return {__builtin_shufflevector(x.vec_, x.vec_, 2, 3, 0, 1, 6, 7, 4, 5)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1];
    ret.data_[0] = __builtin_shufflevector(a, a, 2, 3, 0, 1);
    ret.data_[1] = __builtin_shufflevector(b, b, 2, 3, 0, 1);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1], c = x.data_[2],
                    d = x.data_[3];
    ret.data_[0] = b;
    ret.data_[1] = a;
    ret.data_[2] = d;
    ret.data_[3] = c;
    return ret;
  }
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse2(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8);
  if constexpr (W == 8) {
    return {__builtin_shufflevector(x.vec_, x.vec_, 1, 0, 3, 2, 5, 4, 7, 6)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1];
    ret.data_[0] = __builtin_shufflevector(a, a, 1, 0, 3, 2);
    ret.data_[1] = __builtin_shufflevector(b, b, 1, 0, 3, 2);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1], c = x.data_[2],
                    d = x.data_[3];
    ret.data_[0] = __builtin_shufflevector(a, a, 1, 0);
    ret.data_[1] = __builtin_shufflevector(b, b, 1, 0);
    ret.data_[2] = __builtin_shufflevector(c, c, 1, 0);
    ret.data_[3] = __builtin_shufflevector(d, d, 1, 0);
    return ret;
  }
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse4(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8);
  if constexpr (W == 8) {
    return {__builtin_shufflevector(x.vec_, x.vec_, 3, 2, 1, 0, 7, 6, 5, 4)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1];
    ret.data_[0] = __builtin_shufflevector(a, a, 3, 2, 1, 0);
    ret.data_[1] = __builtin_shufflevector(b, b, 3, 2, 1, 0);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1], c = x.data_[2],
                    d = x.data_[3];
    ret.data_[0] = __builtin_shufflevector(b, b, 1, 0);
    ret.data_[1] = __builtin_shufflevector(a, a, 1, 0);
    ret.data_[2] = __builtin_shufflevector(d, d, 1, 0);
    ret.data_[3] = __builtin_shufflevector(c, c, 1, 0);
    return ret;
  }
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto reverse8(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8);
  if constexpr (W == 8) {
    return {__builtin_shufflevector(x.vec_, x.vec_, 7, 6, 5, 4, 3, 2, 1, 0)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1];
    ret.data_[0] = __builtin_shufflevector(b, b, 3, 2, 1, 0);
    ret.data_[1] = __builtin_shufflevector(a, a, 3, 2, 1, 0);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1], c = x.data_[2],
                    d = x.data_[3];
    ret.data_[0] = __builtin_shufflevector(d, d, 1, 0);
    ret.data_[1] = __builtin_shufflevector(c, c, 1, 0);
    ret.data_[2] = __builtin_shufflevector(b, b, 1, 0);
    ret.data_[3] = __builtin_shufflevector(a, a, 1, 0);
    return ret;
  }
}
template <int _0, int _1, int _2, int _3, std::ptrdiff_t N, std::ptrdiff_t W,
          typename T>
[[gnu::always_inline]] auto vshuf(simd::Unroll<1, N, W, T> x,
                                  simd::Unroll<1, N, W, T> y)
  -> simd::Unroll<1, N, W, T> {
  if constexpr (W == 8) {
    return {__builtin_shufflevector(x.vec_, y.vec_, _0, _1, 8 + _2, 8 + _3,
                                    4 + _0, 4 + _1, 12 + _2, 12 + _3)};
  } else if constexpr (W == 4) {
    simd::Vec<W, T> xa = x.data_[0], xb = x.data_[1], ya = y.data_[0],
                    yb = y.data_[1];
    simd::Unroll<1, N, W, T> ret;
    ret.data_[0] = __builtin_shufflevector(xa, ya, _0, _1, 4 + _2, 4 + _3);
    ret.data_[1] = __builtin_shufflevector(xb, yb, _0, _1, 4 + _2, 4 + _3);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Vec<W, T> xa = x.data_[0], xb = x.data_[1], xc = x.data_[2],
                    xd = x.data_[3], ya = y.data_[0], yb = y.data_[1],
                    yc = y.data_[2], yd = y.data_[3];
    simd::Unroll<1, N, W, T> ret;
    ret.data_[0] = __builtin_shufflevector(xa, xb, _0, _1);
    ret.data_[1] = __builtin_shufflevector(ya, yb, _2, _3);
    ret.data_[2] = __builtin_shufflevector(xc, xd, _0, _1);
    ret.data_[3] = __builtin_shufflevector(yc, yd, _2, _3);
    return ret;
  }
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto interleave_halves(simd::Unroll<1, N, W, T> x)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8);
  if constexpr (W == 8) {
    return {__builtin_shufflevector(x.vec_, x.vec_, 0, 4, 1, 5, 2, 6, 3, 7)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1];
    ret.data_[0] = __builtin_shufflevector(a, b, 0, 4, 1, 5);
    ret.data_[1] = __builtin_shufflevector(a, b, 2, 6, 3, 7);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> a = x.data_[0], b = x.data_[1], c = x.data_[2],
                    d = x.data_[3];
    ret.data_[0] = __builtin_shufflevector(a, c, 0, 2);
    ret.data_[1] = __builtin_shufflevector(a, c, 1, 3);
    ret.data_[2] = __builtin_shufflevector(b, d, 0, 2);
    ret.data_[3] = __builtin_shufflevector(b, d, 1, 3);
    return ret;
  }
}
template <std::ptrdiff_t N, std::ptrdiff_t W, typename T>
[[gnu::always_inline]] auto interleave(simd::Unroll<1, N, W, T> x,
                                       simd::Unroll<1, N, W, T> y)
  -> simd::Unroll<1, N, W, T> {
  static_assert(N * W == 8);
  if constexpr (W == 8) {
    return {__builtin_shufflevector(x.vec_, y.vec_, 0, 9, 2, 11, 4, 13, 6, 15)};
  } else if constexpr (W == 4) {
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> xa = x.data_[0], xb = x.data_[1], ya = y.data_[0],
                    yb = y.data_[1];
    ret.data_[0] = __builtin_shufflevector(xa, ya, 0, 5, 2, 7);
    ret.data_[1] = __builtin_shufflevector(xb, yb, 0, 5, 2, 7);
    return ret;
  } else {
    static_assert(W == 2);
    simd::Unroll<1, N, W, T> ret;
    simd::Vec<W, T> xa = x.data_[0], xb = x.data_[1], xc = x.data_[2],
                    xd = x.data_[3], ya = y.data_[0], yb = y.data_[1],
                    yc = y.data_[2], yd = y.data_[3];
    ret.data_[0] = __builtin_shufflevector(xa, ya, 0, 3);
    ret.data_[1] = __builtin_shufflevector(xb, yb, 0, 3);
    ret.data_[2] = __builtin_shufflevector(xc, yc, 0, 3);
    ret.data_[3] = __builtin_shufflevector(xd, yd, 0, 3);
    return ret;
  }
}
template <typename T>
auto sortImpl(math::SVector<T, 16> x) -> math::SVector<T, 16> {
  auto [xlo, xhi] = x.split();
  V8<T> a = xlo.simd(), b = xhi.simd();
  minmax(a, b);

  b = reverse2(b);
  minmax(a, b);

  V8<T> c = a;
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

  V8<T> b2 = reverse2(b);
  V8<T> b1 = reverse2(a);

  a = interleave(a, b2);
  b = interleave(b1, b);

  using SV8 = math::SVector<T, 8>;
  return SV8(a).cat(SV8(b));
}
} // namespace

namespace utils {

auto sort(math::SVector<double, 16> x) -> math::SVector<double, 16> {
  return sortImpl(x);
}
auto sort(math::SVector<float, 16> x) -> math::SVector<float, 16> {
  return sortImpl(x);
}

} // namespace utils

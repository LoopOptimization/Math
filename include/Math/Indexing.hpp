#pragma once
#include "Math/AxisTypes.hpp"
#include "Math/Iterators.hpp"
#include "Math/MatrixDimensions.hpp"
#include "SIMD/Indexing.hpp"
#include <cstddef>

namespace poly::math {

[[maybe_unused]] static inline constexpr struct Begin {
  friend inline auto operator<<(std::ostream &os, Begin) -> std::ostream & {
    return os << 0;
  }
} begin;
[[maybe_unused]] static inline constexpr struct End {
  friend inline auto operator<<(std::ostream &os, End) -> std::ostream & {
    return os << "end";
  }
} end;
/// TODO: remove `OffsetBegin`
/// We probably won't support non-zero-based indexing
struct OffsetBegin {
  [[no_unique_address]] ptrdiff_t offset;
  friend inline auto operator<<(std::ostream &os,
                                OffsetBegin r) -> std::ostream & {
    return os << r.offset;
  }
};

constexpr auto operator+(ptrdiff_t x, Begin) -> OffsetBegin {
  return OffsetBegin{x};
}
constexpr auto operator+(Begin, ptrdiff_t x) -> OffsetBegin {
  return OffsetBegin{x};
}
constexpr auto operator+(ptrdiff_t x, OffsetBegin y) -> OffsetBegin {
  return OffsetBegin{x + y.offset};
}
constexpr auto operator+(OffsetBegin y, ptrdiff_t x) -> OffsetBegin {
  return OffsetBegin{ptrdiff_t(x) + y.offset};
}
[[maybe_unused]] static constexpr inline struct OffsetEnd {
  [[no_unique_address]] ptrdiff_t offset;
  friend inline auto operator<<(std::ostream &os,
                                OffsetEnd r) -> std::ostream & {
    return os << "end - " << r.offset;
  }
} last{1};
constexpr auto operator-(End, ptrdiff_t x) -> OffsetEnd { return OffsetEnd{x}; }
constexpr auto operator-(OffsetEnd y, ptrdiff_t x) -> OffsetEnd {
  return OffsetEnd{y.offset + x};
}
constexpr auto operator+(OffsetEnd y, ptrdiff_t x) -> OffsetEnd {
  return OffsetEnd{y.offset - x};
}

// Union type
template <typename T>
concept ScalarRelativeIndex =
  std::same_as<T, End> || std::same_as<T, Begin> ||
  std::same_as<T, OffsetBegin> || std::same_as<T, OffsetEnd>;

template <typename T>
concept ScalarIndex =
  std::convertible_to<T, ptrdiff_t> || ScalarRelativeIndex<T>;

[[maybe_unused]] static constexpr inline struct Colon {
  [[nodiscard]] inline constexpr auto operator()(auto B, auto E) const {
    return Range{standardizeRangeBound(B), standardizeRangeBound(E)};
  }
} _;

constexpr auto canonicalize(ptrdiff_t e, ptrdiff_t) -> ptrdiff_t { return e; }
constexpr auto canonicalize(Begin, ptrdiff_t) -> ptrdiff_t { return 0; }
constexpr auto canonicalize(OffsetBegin b, ptrdiff_t) -> ptrdiff_t {
  return b.offset;
}
constexpr auto canonicalize(End, ptrdiff_t M) -> ptrdiff_t { return M; }
constexpr auto canonicalize(OffsetEnd e, ptrdiff_t M) -> ptrdiff_t {
  return M - e.offset;
}
template <typename B, typename E>
constexpr auto canonicalizeRange(Range<B, E> r,
                                 ptrdiff_t M) -> Range<ptrdiff_t, ptrdiff_t> {
  return Range<ptrdiff_t, ptrdiff_t>{canonicalize(r.b, M),
                                     canonicalize(r.e, M)};
}
constexpr auto canonicalizeRange(Colon,
                                 ptrdiff_t M) -> Range<ptrdiff_t, ptrdiff_t> {
  return Range<ptrdiff_t, ptrdiff_t>{0, M};
}

static_assert(ScalarIndex<OffsetEnd>);

template <typename T>
concept AbstractSlice = requires(T t, ptrdiff_t M) {
  { canonicalizeRange(t, M) } -> std::same_as<Range<ptrdiff_t, ptrdiff_t>>;
};
static_assert(AbstractSlice<Range<ptrdiff_t, ptrdiff_t>>);
static_assert(AbstractSlice<Colon>);

[[nodiscard]] inline constexpr auto calcOffset(Length<> len,
                                               ptrdiff_t i) -> ptrdiff_t {
  invariant(i >= 0z);
  invariant(i < len);
  return i;
}
[[nodiscard]] inline constexpr auto calcOffset(Length<1>,
                                               ptrdiff_t i) -> ptrdiff_t {
  invariant(i == 0z);
  return 0z;
}
// FIXME: probably not needed
[[nodiscard]] inline constexpr auto
calcOffset(std::integral_constant<ptrdiff_t, 1>, ptrdiff_t i) -> ptrdiff_t {
  invariant(i == 0z);
  return 0z;
}
[[nodiscard]] inline constexpr auto calcOffset(Length<>, Begin) -> ptrdiff_t {
  return 0z;
}
[[nodiscard]] inline constexpr auto calcOffset(Length<> len,
                                               OffsetBegin i) -> ptrdiff_t {
  return calcOffset(len, i.offset);
}
[[nodiscard]] inline constexpr auto calcOffset(Length<> len,
                                               OffsetEnd i) -> ptrdiff_t {
  invariant(i.offset <= len);
  return ptrdiff_t(len) - i.offset;
}
[[nodiscard]] inline constexpr auto calcRangeOffset(Length<> len,
                                                    ptrdiff_t i) -> ptrdiff_t {
  invariant(i <= len);
  return i;
}
[[nodiscard]] inline constexpr auto calcRangeOffset(Length<>,
                                                    Begin) -> ptrdiff_t {
  return 0z;
}
[[nodiscard]] inline constexpr auto
calcRangeOffset(Length<> len, OffsetBegin i) -> ptrdiff_t {
  return calcRangeOffset(len, i.offset);
}
[[nodiscard]] inline constexpr auto calcRangeOffset(Length<> len,
                                                    OffsetEnd i) -> ptrdiff_t {
  invariant(i.offset <= len);
  return ptrdiff_t(len) - i.offset;
}
// note that we don't check i.b < len because we want to allow
// empty ranges, and r.b <= r.e <= len is checked in calcNewDim.
template <class B, class E>
constexpr auto calcOffset(Length<> len, Range<B, E> i) -> ptrdiff_t {
  return calcRangeOffset(len, i.b);
}
[[nodiscard, gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(Length<>, Colon) -> ptrdiff_t {
  return 0z;
}

[[nodiscard, gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(SquareDims<>, ptrdiff_t i) -> ptrdiff_t {
  return i;
}
[[nodiscard, gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(DenseDims<>, ptrdiff_t i) -> ptrdiff_t {
  return i;
}

struct StridedRange {
  Length<> len;
  RowStride<> stride;
  [[gnu::artificial, gnu::always_inline]] explicit inline constexpr
  operator ptrdiff_t() const {
    return ptrdiff_t(len);
  }
  friend inline auto operator<<(std::ostream &os,
                                StridedRange x) -> std::ostream & {
    return os << "Length: " << ptrdiff_t(x.len)
              << " (stride: " << ptrdiff_t(x.stride) << ")";
  }
};

template <ptrdiff_t U, ptrdiff_t W, typename M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(Length<> len, simd::index::Unroll<U, W, M> i) {
  if constexpr (std::same_as<M, simd::mask::None<W>>)
    invariant((i.index + U * W - 1) < len);
  else invariant(i.index + (U - 1) * W + i.mask.lastUnmasked() - 1 < len);
  return i.index;
}

template <class I>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(StridedRange d, I i) -> ptrdiff_t {
  return ptrdiff_t(d.stride) * calcOffset(d.len, i);
}

template <class R, class C>
[[nodiscard, gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(StridedDims<> d, R r, C c) -> ptrdiff_t {
  return ptrdiff_t(stride(d)) * calcOffset(length(ptrdiff_t(Row<>(d))), r) +
         calcOffset(length(ptrdiff_t(Col<>(d))), c);
}

// constexpr auto is_integral_const(auto) -> bool { return false; }
// template <typename T, T V>
// constexpr auto is_integral_const(std::integral_constant<T, V>) -> bool {
//   return true;
// }
[[nodiscard, gnu::artificial, gnu::always_inline]] inline constexpr auto
row(StridedRange r) -> Row<> {
  return row(ptrdiff_t(r.len));
}
[[nodiscard, gnu::artificial, gnu::always_inline]] inline constexpr auto
col(StridedRange) -> Col<1> {
  return {};
}
[[nodiscard, gnu::artificial, gnu::always_inline]] inline constexpr auto
stride(StridedRange r) -> RowStride<> {
  return r.stride;
}

template <typename T>
concept StaticInt =
  std::is_same_v<T, std::integral_constant<typename T::value_type, T::value>>;

template <typename T>
concept DenseLayout =
  std::integral<T> || std::is_convertible_v<T, DenseDims<>> || StaticInt<T>;

static_assert(StaticInt<std::integral_constant<ptrdiff_t, 3>>);
static_assert(!StaticInt<int64_t>);

template <ptrdiff_t R, ptrdiff_t C, ptrdiff_t X>
constexpr auto IsStridedColVectorDim(StridedDims<R, C, X>) -> bool {
  return C == 1;
}
constexpr auto IsStridedColVectorDim(auto) -> bool { return false; }

template <typename D> consteval auto IsStridedColVectorDim() -> bool {
  return IsStridedColVectorDim(std::declval<D>());
}

template <typename T>
concept IsOne =
  std::same_as<std::remove_cvref_t<T>, std::integral_constant<ptrdiff_t, 1>>;

template <typename T>
concept RowVectorDimension = requires(T t) {
  { Length(t) } -> std::same_as<T>;
};
static_assert(RowVectorDimension<Length<3>>);
static_assert(RowVectorDimension<Length<>>);
static_assert(!RowVectorDimension<ptrdiff_t>);

template <typename T>
concept StaticLength = RowVectorDimension<T> && !std::same_as<T, Length<>>;

template <typename D>
concept ColVectorSMatDimension =
  std::same_as<decltype(Col(std::declval<D>())), Col<1>>;
template <typename D>
concept ColVectorDimension =
  std::same_as<D, StridedRange> || ColVectorSMatDimension<D>;

// constexpr auto row(RowVectorDimension auto) -> Row<1> { return {}; }
constexpr auto row(ColVectorSMatDimension auto s) { return Row(s); }
constexpr auto col(ColVectorSMatDimension auto) -> Col<1> { return {}; }
constexpr auto stride(ColVectorSMatDimension auto) -> RowStride<1> {
  return {};
}
template <class I>
constexpr auto calcOffset(ColVectorSMatDimension auto d, I i) -> ptrdiff_t {
  return unwrapStride(stride(d)) * calcOffset(unwrapRow(Row(d)), i);
}

template <typename D>
concept VectorDimension = RowVectorDimension<D> || ColVectorDimension<D>;
// Concept for aligning array dimensions with indices.
template <class I, class D>
concept Index =
  (VectorDimension<D> &&
   (ScalarIndex<I> || AbstractSlice<I> || simd::index::issimd<I>)) ||
  (DenseLayout<D> && (ScalarIndex<I> || simd::index::issimd<I>)) ||
  (MatrixDimension<D> && requires(I i) {
    { i.rowIdx };
    { i.colIdx };
  });
static_assert(Index<CartesianIndex<ptrdiff_t, ptrdiff_t>, DenseDims<>>);
struct Empty {};

[[gnu::artificial, gnu::always_inline]] inline constexpr auto
calcNewDim(VectorDimension auto, ScalarIndex auto) -> Empty {
  return {};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
calcNewDim(SquareDims<>, ptrdiff_t) -> Empty {
  return {};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
calcNewDim(DenseDims<>, ptrdiff_t) -> Empty {
  return {};
}
// constexpr auto calcNewDim(auto, ptrdiff_t) -> Empty { return {}; }

constexpr auto calcNewDim(Length<> len,
                          Range<ptrdiff_t, ptrdiff_t> r) -> Length<> {
  invariant(r.e <= len);
  invariant(r.b <= r.e);
  return length(ptrdiff_t(r.e - r.b));
}
template <class B, class E>
constexpr auto calcNewDim(Length<> len, Range<B, E> r) -> Length<> {
  return calcNewDim(len, canonicalizeRange(r, ptrdiff_t(len)));
}
constexpr auto calcNewDim(StridedRange len,
                          Range<ptrdiff_t, ptrdiff_t> r) -> StridedRange {
  return StridedRange{calcNewDim(len.len, r), len.stride};
}
template <class B, class E>
constexpr auto calcNewDim(StridedRange len, Range<B, E> r) -> StridedRange {
  return StridedRange{calcNewDim(len.len, r), len.stride};
}
template <ScalarIndex R, ScalarIndex C>
constexpr auto calcNewDim(StridedDims<>, R, C) -> Empty {
  return {};
}
constexpr auto calcNewDim(Length<> len, Colon) -> Length<> { return len; };
constexpr auto calcNewDim(StaticInt auto len, Colon) { return len; };
constexpr auto calcNewDim(StridedRange len, Colon) -> StridedRange {
  return len;
};

template <AbstractSlice B, ScalarIndex C>
constexpr auto calcNewDim(StridedDims<> d, B b, C) -> StridedRange {
  Length<> rowDims = calcNewDim(length(ptrdiff_t(Row(d))), b);
  return StridedRange{rowDims, RowStride(d)};
}

template <ScalarIndex R, AbstractSlice C>
constexpr auto calcNewDim(StridedDims<> d, R, C c) {
  return calcNewDim(length(ptrdiff_t(Col(d))), c);
}

template <AbstractSlice B, AbstractSlice C>
constexpr auto calcNewDim(StridedDims<> d, B r, C c) {
  auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
  auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
  return StridedDims(row(rowDims), col(colDims), RowStride(d));
}
template <ptrdiff_t NR, ptrdiff_t NC, AbstractSlice B, AbstractSlice C>
constexpr auto calcNewDim(DenseDims<NR, NC> d, B r, C c) {
  if constexpr ((NR >= 0) && std::same_as<B, Colon>) {
    if constexpr ((NC >= 0) && std::same_as<C, Colon>) {
      return DenseDims(Row<NR>{}, Col<NC>{});
    } else {
      auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
      return StridedDims(Row<NR>{}, col(colDims), rowStride(unwrapCol(Col(d))));
    }
  } else if constexpr ((NC >= 0) && std::same_as<C, Colon>) {
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return DenseDims(row(rowDims), Col<NC>{});
  } else {
    auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return StridedDims(row(rowDims), col(colDims), RowStride(d));
  }
}
template <ptrdiff_t NR, ptrdiff_t NC, ptrdiff_t X, AbstractSlice B,
          AbstractSlice C>
constexpr auto calcNewDim(StridedDims<NR, NC, X> d, B r, C c) {
  if constexpr ((NR >= 0) && std::same_as<B, Colon>) {
    if constexpr ((NC >= 0) && std::same_as<C, Colon>) {
      return StridedDims(Row<NR>{}, Col<NC>{}, RowStride(d));
    } else {
      auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
      return StridedDims(Row<NR>{}, col(colDims), RowStride(d));
    }
  } else if constexpr ((NC >= 0) && std::same_as<C, Colon>) {
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return StridedDims(row(rowDims), Col<NC>{}, RowStride(d));
  } else {
    auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return StridedDims(row(rowDims), col(colDims), RowStride(d));
  }
}
template <AbstractSlice B>
constexpr auto calcNewDim(DenseDims<> d, B r, Colon) {
  auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
  return DenseDims(row(rowDims), Col(d));
}
template <AbstractSlice B>
constexpr auto calcNewDim(SquareDims<> d, B r, Colon) {
  auto rowDims = ptrdiff_t(calcNewDim(length(ptrdiff_t(Row(d))), r));
  return DenseDims(row(rowDims), col(d));
}
template <ptrdiff_t R, ptrdiff_t C, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, simd::index::Unroll<R>,
                          simd::index::Unroll<C, W, M> c) {
  return simd::index::UnrollDims<R, C, W, M>{c.mask, RowStride(d)};
}
template <ptrdiff_t R, ptrdiff_t C, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, simd::index::Unroll<C, W, M> r,
                          simd::index::Unroll<R>) {
  return simd::index::UnrollDims<R, C, W, M, true>{r.mask, RowStride(d)};
}

template <ptrdiff_t C, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, ptrdiff_t,
                          simd::index::Unroll<C, W, M> c) {
  return simd::index::UnrollDims<1, C, W, M>{c.mask, RowStride(d)};
}
template <ptrdiff_t R, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, simd::index::Unroll<R, W, M> r,
                          ptrdiff_t) {
  if constexpr (W == 1)
    return simd::index::UnrollDims<R, 1, 1, M>{r.mask, RowStride(d)};
  else return simd::index::UnrollDims<1, R, W, M, true>{r.mask, RowStride(d)};
}

template <ptrdiff_t U, ptrdiff_t W, typename M>
constexpr auto calcNewDim(ptrdiff_t, simd::index::Unroll<U, W, M> i) {
  return simd::index::UnrollDims<1, U, W, M, false, 1>{i.mask, RowStride<1>{}};
}

template <ptrdiff_t U, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedRange x, simd::index::Unroll<U, W, M> i) {
  if constexpr (W == 1)
    return simd::index::UnrollDims<U, 1, 1, M, false, -1>{i.mask, stride(x)};
  else return simd::index::UnrollDims<1, U, W, M, true, -1>{i.mask, stride(x)};
}
template <ptrdiff_t U, ptrdiff_t W, typename M>
constexpr auto calcNewDim(ColVectorSMatDimension auto x,
                          simd::index::Unroll<U, W, M> i) {
  if constexpr (W == 1)
    return simd::index::UnrollDims<U, 1, 1, M, false, -1>{i.mask, stride(x)};
  else return simd::index::UnrollDims<1, U, W, M, true, -1>{i.mask, stride(x)};
}

} // namespace poly::math

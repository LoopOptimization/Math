module;
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <type_traits>

export module Indexing;

import AxisTypes;
import CompressReference;
import Range;
import MatDim;
import SIMD;

export namespace math {

/// TODO: remove `OffsetBegin`
/// We probably won't support non-zero-based indexing
struct OffsetBegin {
  [[no_unique_address]] ptrdiff_t offset_;

private:
  friend auto operator<<(std::ostream &os, OffsetBegin r) -> std::ostream & {
    return os << r.offset_;
  }
  [[gnu::always_inline, gnu::artificial]] friend inline constexpr auto
  operator+(ptrdiff_t x, OffsetBegin y) -> OffsetBegin {
    return OffsetBegin{x + y.offset_};
  }
  [[gnu::always_inline, gnu::artificial]] friend inline constexpr auto
  operator+(OffsetBegin y, ptrdiff_t x) -> OffsetBegin {
    return OffsetBegin{x + y.offset_};
  }
};
[[maybe_unused]] inline constexpr struct Begin {
private:
  friend auto operator<<(std::ostream &os, Begin) -> std::ostream & {
    return os << 0;
  }
  [[gnu::always_inline, gnu::artificial]] friend inline constexpr auto
  operator+(ptrdiff_t x, Begin) -> OffsetBegin {
    return OffsetBegin{x};
  }
  [[gnu::always_inline, gnu::artificial]] friend inline constexpr auto
  operator+(Begin, ptrdiff_t x) -> OffsetBegin {
    return OffsetBegin{x};
  }
} begin;

[[maybe_unused]] constexpr inline struct OffsetEnd {
  [[no_unique_address]] ptrdiff_t offset_;

private:
  friend auto operator<<(std::ostream &os, OffsetEnd r) -> std::ostream & {
    return os << "end - " << r.offset_;
  }
  [[gnu::always_inline, gnu::artificial]] friend inline constexpr auto
  operator-(OffsetEnd y, ptrdiff_t x) -> OffsetEnd {
    return OffsetEnd{y.offset_ + x};
  }
  [[gnu::always_inline, gnu::artificial]] friend inline constexpr auto
  operator+(OffsetEnd y, ptrdiff_t x) -> OffsetEnd {
    return OffsetEnd{y.offset_ - x};
  }

} last{1};
[[maybe_unused]] inline constexpr struct End {
private:
  friend auto operator<<(std::ostream &os, End) -> std::ostream & {
    return os << "end";
  }
  [[gnu::always_inline, gnu::artificial]] friend inline constexpr auto
  operator-(End, ptrdiff_t x) -> OffsetEnd {
    return OffsetEnd{x};
  }
} end;

// Union type
template <typename T>
concept ScalarRelativeIndex =
  std::same_as<T, End> || std::same_as<T, Begin> ||
  std::same_as<T, OffsetBegin> || std::same_as<T, OffsetEnd>;

template <typename T>
concept ScalarIndex =
  std::convertible_to<T, ptrdiff_t> || ScalarRelativeIndex<T>;

[[maybe_unused]] constexpr inline struct Colon {
  [[nodiscard, gnu::always_inline]] constexpr auto operator()(auto B,
                                                              auto E) const {
    return Range{math::standardizeRangeBound(B),
                 math::standardizeRangeBound(E)};
  }
} _;

constexpr auto canonicalize(ptrdiff_t e, ptrdiff_t) -> ptrdiff_t { return e; }
constexpr auto canonicalize(Begin, ptrdiff_t) -> ptrdiff_t { return 0; }
constexpr auto canonicalize(OffsetBegin b, ptrdiff_t) -> ptrdiff_t {
  return b.offset_;
}
constexpr auto canonicalize(End, ptrdiff_t M) -> ptrdiff_t { return M; }
constexpr auto canonicalize(OffsetEnd e, ptrdiff_t M) -> ptrdiff_t {
  return M - e.offset_;
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

[[nodiscard, gnu::always_inline]] constexpr auto
calcOffset(Length<> len, ptrdiff_t i) -> ptrdiff_t {
  invariant(i >= 0z);
  invariant(i < len);
  return i;
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcOffset(Length<1>, ptrdiff_t i) -> ptrdiff_t {
  invariant(i == 0z);
  return 0z;
}
// FIXME: probably not needed
[[nodiscard, gnu::always_inline]] constexpr auto
calcOffset(std::integral_constant<ptrdiff_t, 1>, ptrdiff_t i) -> ptrdiff_t {
  invariant(i == 0z);
  return 0z;
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcOffset(Length<>, Begin) -> ptrdiff_t {
  return 0z;
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcOffset(Length<> len, OffsetBegin i) -> ptrdiff_t {
  return calcOffset(len, i.offset_);
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcOffset(Length<> len, OffsetEnd i) -> ptrdiff_t {
  invariant(i.offset_ <= len);
  return ptrdiff_t(len) - i.offset_;
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcRangeOffset(Length<> len, ptrdiff_t i) -> ptrdiff_t {
  invariant(i <= len);
  return i;
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcRangeOffset(Length<>, Begin) -> ptrdiff_t {
  return 0z;
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcRangeOffset(Length<> len, OffsetBegin i) -> ptrdiff_t {
  return calcRangeOffset(len, i.offset_);
}
[[nodiscard, gnu::always_inline]] constexpr auto
calcRangeOffset(Length<> len, OffsetEnd i) -> ptrdiff_t {
  invariant(i.offset_ <= len);
  return ptrdiff_t(len) - i.offset_;
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

template <ptrdiff_t L = -1, ptrdiff_t X = -1> struct StridedRange {
  Length<L> len_;
  RowStride<X> stride_;
  [[gnu::artificial, gnu::always_inline]] explicit inline constexpr
  operator ptrdiff_t() const {
    return ptrdiff_t(len_);
  }
  constexpr explicit operator Row<L>() const { return asrow(len_); }
  constexpr explicit operator Col<1>() const { return {}; }
  constexpr explicit operator RowStride<X>() const { return stride_; }

private:
  [[nodiscard, gnu::artificial,
    gnu::always_inline]] friend inline constexpr auto
  row(StridedRange<> r) -> Row<> {
    return row(ptrdiff_t(r.len_));
  }
  [[nodiscard, gnu::artificial,
    gnu::always_inline]] friend inline constexpr auto
  col(StridedRange<>) -> Col<1> {
    return {};
  }
  [[nodiscard, gnu::artificial,
    gnu::always_inline]] friend inline constexpr auto
  stride(StridedRange<> r) -> RowStride<> {
    return r.stride_;
  }

  friend inline auto operator<<(std::ostream &os,
                                StridedRange x) -> std::ostream & {
    return os << "Length: " << ptrdiff_t(x.len_)
              << " (stride: " << ptrdiff_t(x.stride_) << ")";
  }
};
template <ptrdiff_t L, ptrdiff_t X> Row(StridedRange<L, X>) -> Row<L>;
template <ptrdiff_t L, ptrdiff_t X> Col(StridedRange<L, X>) -> Col<1z>;
template <ptrdiff_t L, ptrdiff_t X>
RowStride(StridedRange<L, X>) -> RowStride<X>;

static_assert(
  std::same_as<Col<1>, decltype(Col(std::declval<StridedRange<>>()))>);
static_assert(ColVectorDimension<StridedRange<>>);

template <ptrdiff_t U, ptrdiff_t W, typename M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(Length<> len, simd::index::Unroll<U, W, M> i) {
  if constexpr (std::same_as<M, simd::mask::None<W>>)
    invariant((i.index_ + U * W - 1) < len);
  else invariant(i.index_ + (U - 1) * W + i.mask_.lastUnmasked() - 1 < len);
  return i.index_;
}

template <class I>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
calcOffset(StridedRange<> d, I i) -> ptrdiff_t {
  return ptrdiff_t(d.stride_) * calcOffset(d.len_, i);
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
concept DenseLayout =
  RowVectorDimension<T> || std::is_convertible_v<T, DenseDims<>>;

template <typename T>
concept StaticLength = RowVectorDimension<T> && !std::same_as<T, Length<>>;

// constexpr auto row(RowVectorDimension auto) -> Row<1> { return {}; }
// constexpr auto row(ColVectorDimension auto s) { return Row(s); }
// constexpr auto col(ColVectorDimension auto) -> Col<1> { return {}; }
// constexpr auto stride(ColVectorDimension auto) -> RowStride<1> { return {}; }
template <class I>
constexpr auto calcOffset(ColVectorDimension auto d, I i) -> ptrdiff_t {
  return unwrapStride(stride(d)) * calcOffset(length(unwrapRow(Row(d))), i);
}

// Concept for aligning array dimensions with indices.
template <class I, class D>
concept Index =
  (VectorDimension<D> &&
   (ScalarIndex<I> || AbstractSlice<I> || simd::index::issimd<I>)) ||
  (DenseLayout<D> && (ScalarIndex<I> || simd::index::issimd<I>)) ||
  (MatrixDimension<D> && requires(I i) {
    { i.row_idx_ };
    { i.col_idx_ };
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
constexpr auto calcNewDim(StridedRange<> len,
                          Range<ptrdiff_t, ptrdiff_t> r) -> StridedRange<> {
  return StridedRange<>{calcNewDim(len.len_, r), len.stride_};
}
template <class B, class E>
constexpr auto calcNewDim(StridedRange<> len, Range<B, E> r) -> StridedRange<> {
  return StridedRange<>{calcNewDim(len.len_, r), len.stride_};
}
template <ScalarIndex R, ScalarIndex C>
constexpr auto calcNewDim(StridedDims<>, R, C) -> Empty {
  return {};
}
constexpr auto calcNewDim(Length<> len, Colon) -> Length<> { return len; };
constexpr auto calcNewDim(StaticInt auto len, Colon) { return len; };
constexpr auto calcNewDim(StridedRange<> len, Colon) -> StridedRange<> {
  return len;
};

template <AbstractSlice B, ScalarIndex C>
constexpr auto calcNewDim(StridedDims<> d, B b, C) -> StridedRange<> {
  Length<> rowDims = calcNewDim(length(ptrdiff_t(Row(d))), b);
  return StridedRange<>{rowDims, RowStride(d)};
}

template <ScalarIndex R, AbstractSlice C>
constexpr auto calcNewDim(StridedDims<> d, R, C c) {
  return calcNewDim(length(ptrdiff_t(Col(d))), c);
}

template <AbstractSlice B, AbstractSlice C>
constexpr auto calcNewDim(StridedDims<> d, B r, C c) {
  auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
  auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
  return StridedDims(asrow(rowDims), ascol(colDims), RowStride(d));
}
template <ptrdiff_t NR, ptrdiff_t NC, AbstractSlice B, AbstractSlice C>
constexpr auto calcNewDim(DenseDims<NR, NC> d, B r, C c) {
  if constexpr ((NR >= 0) && std::same_as<B, Colon>) {
    if constexpr ((NC >= 0) && std::same_as<C, Colon>) {
      return DenseDims(Row<NR>{}, Col<NC>{});
    } else {
      auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
      return StridedDims(Row<NR>{}, ascol(colDims),
                         stride(unwrapCol(Col(d))));
    }
  } else if constexpr ((NC >= 0) && std::same_as<C, Colon>) {
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return DenseDims(asrow(rowDims), Col<NC>{});
  } else {
    auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return StridedDims(asrow(rowDims), ascol(colDims), RowStride(d));
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
      return StridedDims(Row<NR>{}, ascol(colDims), RowStride(d));
    }
  } else if constexpr ((NC >= 0) && std::same_as<C, Colon>) {
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return StridedDims(asrow(rowDims), Col<NC>{}, RowStride(d));
  } else {
    auto colDims = calcNewDim(length(ptrdiff_t(Col(d))), c);
    auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
    return StridedDims(asrow(rowDims), ascol(colDims), RowStride(d));
  }
}
template <AbstractSlice B>
constexpr auto calcNewDim(DenseDims<> d, B r, Colon) {
  auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
  return DenseDims(asrow(rowDims), Col(d));
}
template <AbstractSlice B>
constexpr auto calcNewDim(SquareDims<> d, B r, Colon) {
  auto rowDims = calcNewDim(length(ptrdiff_t(Row(d))), r);
  return DenseDims(asrow(rowDims), col(d));
}
template <ptrdiff_t R, ptrdiff_t C, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, simd::index::Unroll<R>,
                          simd::index::Unroll<C, W, M> c) {
  return simd::index::UnrollDims<R, C, W, M>{c.mask_, RowStride(d)};
}
template <ptrdiff_t R, ptrdiff_t C, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, simd::index::Unroll<C, W, M> r,
                          simd::index::Unroll<R>) {
  return simd::index::UnrollDims<R, C, W, M, true>{r.mask_, RowStride(d)};
}

template <ptrdiff_t C, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, ptrdiff_t,
                          simd::index::Unroll<C, W, M> c) {
  return simd::index::UnrollDims<1, C, W, M>{c.mask_, RowStride(d)};
}
template <ptrdiff_t R, ptrdiff_t W, typename M>
constexpr auto calcNewDim(StridedDims<> d, simd::index::Unroll<R, W, M> r,
                          ptrdiff_t) {
  if constexpr (W == 1)
    return simd::index::UnrollDims<R, 1, 1, M>{r.mask_, RowStride(d)};
  else return simd::index::UnrollDims<1, R, W, M, true>{r.mask_, RowStride(d)};
}

template <ptrdiff_t U, ptrdiff_t W, typename M>
constexpr auto calcNewDim(Length<>, simd::index::Unroll<U, W, M> i) {
  return simd::index::UnrollDims<1, U, W, M, false, 1>{i.mask_, RowStride<1>{}};
}

template <ptrdiff_t U, ptrdiff_t W, typename M>
constexpr auto calcNewDim(ColVectorDimension auto x,
                          simd::index::Unroll<U, W, M> i) {
  if constexpr (W == 1)
    return simd::index::UnrollDims<U, 1, 1, M, false, -1>{i.mask_,
                                                          stride(x)};
  else
    return simd::index::UnrollDims<1, U, W, M, true, -1>{i.mask_, stride(x)};
}

} // namespace math

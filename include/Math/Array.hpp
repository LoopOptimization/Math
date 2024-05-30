#pragma once

#include "Alloc/Arena.hpp"
#include "Alloc/Mallocator.hpp"
#include "Containers/Pair.hpp"
#include "Containers/Storage.hpp"
#include "Math/ArrayOps.hpp"
#include "Math/AxisTypes.hpp"
#include "Math/Indexing.hpp"
#include "Math/Iterators.hpp"
#include "Math/Matrix.hpp"
#include "Math/MatrixDimensions.hpp"
#include "Math/Rational.hpp"
#include "SIMD/Intrin.hpp"
#include "SIMD/UnrollIndex.hpp"
#include "Utilities/ArrayPrint.hpp"
#include "Utilities/Assign.hpp"
#include "Utilities/LoopMacros.hpp"
#include "Utilities/Optional.hpp"
#include "Utilities/Reference.hpp"
#include "Utilities/TypeCompression.hpp"
#include "Utilities/TypePromotion.hpp"
#include "Utilities/Valid.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <ranges>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <version>

// https://llvm.org/doxygen/Compiler_8h_source.html#l00307
#ifndef POLY_MATH_HAS_CPP_ATTRIBUTE
#if defined(__cplusplus) && defined(__has_cpp_attribute)
#define POLY_MATH_HAS_CPP_ATTRIBUTE(x) __has_cpp_attribute(x)
#else
#define POLY_MATH_HAS_CPP_ATTRIBUTE(x) 0
#endif
#endif
#if POLY_MATH_HAS_CPP_ATTRIBUTE(gsl::Owner)
#define POLY_MATH_GSL_OWNER [[gsl::Owner]]
#else
#define POLY_MATH_GSL_OWNER
#endif
/// POLY_MATH_GSL_POINTER - Apply this to non-owning classes like
/// StringRef to enable lifetime warnings.
#if POLY_MATH_HAS_CPP_ATTRIBUTE(gsl::Pointer)
#define POLY_MATH_GSL_POINTER [[gsl::Pointer]]
#else
#define POLY_MATH_GSL_POINTER
#endif

namespace poly::math {
using axis::unwrapRow, axis::unwrapCol;

static_assert(Dimension<Length<>>);
static_assert(Dimension<DenseDims<3, 2>>);
static_assert(
  std::same_as<Col<1>, decltype(Col(std::declval<DenseDims<3, 1>>()))>);
static_assert(Dimension<DenseDims<3, 1>>);
static_assert(VectorDimension<StridedRange>);
static_assert(VectorDimension<DenseDims<3, 1>>);
static_assert(!MatrixDimension<DenseDims<3, 1>>);
#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)
template <typename T>
using DefaultAlloc = std::allocator<utils::compressed_t<T>>;
#else
template <typename T>
using DefaultAlloc = alloc::Mallocator<utils::compressed_t<T>>;
#endif

template <class T, Dimension S,
          ptrdiff_t N =
            containers::PreAllocStorage<utils::compressed_t<T>, S>(),
          class A = DefaultAlloc<T>>
struct ManagedArray;

using utils::Valid, utils::Optional;

template <class T, Dimension S, bool Compress = utils::Compressible<T>>
struct Array;

template <typename T>
inline auto printMatrix(std::ostream &os,
                        Array<T, StridedDims<>> A) -> std::ostream &;

// Cases we need to consider:
// 1. Slice-indexing
// 2.a. `ptrdiff_t` indexing, not compressed
// 2.b. `ptrdiff_t` indexing, compressed
// 3.a.i. Vector indexing, contig, no mask
// 3.a.ii. Vector indexing, contig, mask
// 3.b.i. Vector indexing, discontig, no mask
// 3.b.ii. Vector indexing, discontig, mask
// all of the above for `T*` and `const T*`
template <typename T, typename P, typename S, typename I>
[[gnu::flatten, gnu::always_inline]] constexpr auto
index(P *ptr, S shape, I i) noexcept -> decltype(auto) {
  auto offset = calcOffset(shape, i);
  auto new_dim = calcNewDim(shape, i);
  invariant(ptr != nullptr);
  using D = decltype(new_dim);
  if constexpr (simd::index::issimd<D>) return simd::ref(ptr + offset, new_dim);
  else {
    constexpr bool compress = !std::same_as<T, std::remove_const_t<P>>;
    if constexpr (!std::same_as<D, Empty>)
      if constexpr (std::is_const_v<P>)
        return Array<T, D, compress>{ptr + offset, new_dim};
      else return MutArray<T, D, compress>{ptr + offset, new_dim};
    else if constexpr (!compress) return ptr[offset];
    else if constexpr (std::is_const_v<P>) return T::decompress(ptr + offset);
    else return utils::Reference<T>{ptr + offset};
  }
}
// for (row/col)vectors, we drop the row/col, essentially broadcasting
template <typename T, typename P, typename S, typename R, typename C>
[[gnu::flatten, gnu::always_inline]] constexpr auto
index(P *ptr, S shape, R wr, C wc) noexcept -> decltype(auto) {
  if constexpr (MatrixDimension<S>) {
    auto r = unwrapRow(wr);
    auto c = unwrapCol(wc);
    auto offset = calcOffset(shape, r, c);
    auto new_dim = calcNewDim(shape, r, c);
    using D = decltype(new_dim);
    if constexpr (simd::index::issimd<D>) {
      return simd::ref(ptr + offset, new_dim);
    } else {
      constexpr bool compress = !std::same_as<T, std::remove_const_t<P>>;
      if constexpr (!std::same_as<D, Empty>)
        if constexpr (std::is_const_v<P>)
          return Array<T, D, compress>{ptr + offset, new_dim};
        else return MutArray<T, D, compress>{ptr + offset, new_dim};
      else if constexpr (!compress) return ptr[offset];
      else if constexpr (std::is_const_v<P>) return T::decompress(ptr + offset);
      else return utils::Reference<T>{ptr + offset};
    }
  } else if constexpr (ColVectorDimension<S>)
    return index<T>(ptr, shape, unwrapRow(wr));
  else // if constexpr (RowVectorDimension<S>)
    return index<T>(ptr, shape, unwrapCol(wc));
  // else if constexpr (std::integral<R>)
  //   return index<T>(ptr, unwrapRow(Row(shape)), unwrapRow(wr));
  // else
  //   return simd::transpose(index<T>(ptr, unwrapRow(Row(shape)),
  //   unwrapRow(wr)));
}
static_assert(
  std::same_as<StridedDims<3>,
               decltype(calcNewDim(std::declval<DenseDims<3>>(), _, _(1, 5)))>);
static_assert(
  std::same_as<StridedDims<3>, decltype(calcNewDim(
                                 std::declval<StridedDims<3>>(), _, _(1, 5)))>);

template <typename T, bool Column = false> struct SliceIterator {
  using dim_type = std::conditional_t<Column, StridedRange, Length<>>;
  using value_type =
    std::conditional_t<std::is_const_v<T>,
                       Array<std::remove_cvref_t<T>, dim_type>,
                       MutArray<std::remove_reference_t<T>, dim_type>>;
  using storage_type =
    std::conditional_t<std::is_const_v<T>, const utils::compressed_t<T>,
                       utils::compressed_t<T>>;
  storage_type *data_;
  Length<> len_;
  RowStride<> row_stride_;
  ptrdiff_t idx_;
  // constexpr auto operator=(const SliceIterator &) -> SliceIterator & =
  // default;
  // constexpr auto operator*() -> value_type;
  constexpr auto operator*() const -> value_type;
  constexpr auto operator++() -> SliceIterator & {
    idx_++;
    return *this;
  }
  constexpr auto operator++(int) -> SliceIterator {
    SliceIterator ret{*this};
    ++(*this);
    return ret;
  }
  constexpr auto operator--() -> SliceIterator & {
    idx_--;
    return *this;
  }
  constexpr auto operator--(int) -> SliceIterator {
    SliceIterator ret{*this};
    --(*this);
    return ret;
  }
  friend constexpr auto operator-(SliceIterator a,
                                  SliceIterator b) -> ptrdiff_t {
    return a.idx_ - b.idx_;
  }
  friend constexpr auto operator+(SliceIterator a,
                                  ptrdiff_t i) -> SliceIterator<T, Column> {
    return {a.data_, a.len_, a.row_stride_, a.idx_ + i};
  }
  friend constexpr auto operator==(SliceIterator a, SliceIterator b) -> bool {
    return a.idx_ == b.idx_;
  }
  friend constexpr auto operator<=>(SliceIterator a,
                                    SliceIterator b) -> std::strong_ordering {
    return a.idx_ <=> b.idx_;
  }
  friend constexpr auto operator==(SliceIterator a, Row<> r) -> bool
  requires(!Column)
  {
    return a.idx_ == r;
  }
  friend constexpr auto operator<=>(SliceIterator a,
                                    Row<> r) -> std::strong_ordering
  requires(!Column)
  {
    return a.idx_ <=> r;
  }
  friend constexpr auto operator==(SliceIterator a, Col<> r) -> bool
  requires(Column)
  {
    return a.idx_ == r;
  }
  friend constexpr auto operator<=>(SliceIterator a,
                                    Col<> r) -> std::strong_ordering
  requires(Column)
  {
    return a.idx_ <=> r;
  }
};

template <typename T, bool Column = false> struct SliceRange {
  T *data_;
  Length<> len_;
  RowStride<> row_stride_;
  ptrdiff_t stop_;
  [[nodiscard]] constexpr auto begin() const -> SliceIterator<T, Column> {
    return {data_, len_, row_stride_, 0};
  }
  [[nodiscard]] constexpr auto end() const {
    if constexpr (Column) return col(stop_);
    else return row(stop_);
  }
};
/// Constant Array
template <class T, Dimension S, bool Compress>
struct POLY_MATH_GSL_POINTER Array {
  static_assert(!std::is_const_v<T>, "T shouldn't be const");

  using storage_type = std::conditional_t<Compress, utils::compressed_t<T>, T>;
  using value_type = T;
  using reference = T &;
  using const_reference = const T &;
  using size_type = ptrdiff_t;
  using difference_type = int;
  using iterator = storage_type *;
  using const_iterator = const storage_type *;
  using pointer = storage_type *;
  using const_pointer = const storage_type *;
  using concrete = std::true_type;
  static constexpr bool trivial =
    std::is_trivially_default_constructible_v<T> &&
    std::is_trivially_destructible_v<T>;
  static constexpr bool isdense = DenseLayout<S>;
  static constexpr bool flatstride = isdense || std::same_as<S, StridedRange>;
  // static_assert(flatstride != std::same_as<S, StridedDims<>>);

  explicit constexpr Array() = default;
  constexpr Array(const Array &) = default;
  constexpr Array(Array &&) noexcept = default;
  constexpr auto operator=(const Array &) -> Array & = default;
  constexpr auto operator=(Array &&) noexcept -> Array & = default;
  constexpr Array(const storage_type *p, S s) : ptr(p), sz(s) {}
  constexpr Array(Valid<const storage_type> p, S s) : ptr(p), sz(s) {}
  template <ptrdiff_t R, ptrdiff_t C>
  constexpr Array(const storage_type *p, Row<R> r, Col<C> c)
    : ptr(p), sz(S{r, c}) {}
  template <ptrdiff_t R, ptrdiff_t C>
  constexpr Array(Valid<const storage_type> p, Row<R> r, Col<C> c)
    : ptr(p), sz(dimension<S>(r, c)) {}
  template <std::convertible_to<S> V>
  constexpr Array(Array<T, V> a) : ptr(a.data()), sz(a.dim()) {}
  template <size_t N>
  constexpr Array(const std::array<T, N> &a) : ptr(a.data()), sz(length(N)) {}
  [[nodiscard]] constexpr auto data() const noexcept -> const storage_type * {
    invariant(ptr != nullptr || ptrdiff_t(sz) == 0);
    return ptr;
  }
  [[nodiscard]] constexpr auto wrappedPtr() noexcept -> Valid<T> { return ptr; }

  [[nodiscard]] constexpr auto
  begin() const noexcept -> StridedIterator<const T>
  requires(std::is_same_v<S, StridedRange>)
  {
    const storage_type *p = ptr;
    return StridedIterator{p, sz.stride_};
  }
  [[nodiscard]] constexpr auto begin() const noexcept
    -> const storage_type *requires(isdense) { return ptr; }

  [[nodiscard]] constexpr auto end() const noexcept
  requires(flatstride)
  {
    return begin() + ptrdiff_t(sz);
  }
  [[nodiscard]] constexpr auto rbegin() const noexcept
  requires(flatstride)
  {
    return std::reverse_iterator(end());
  }
  [[nodiscard]] constexpr auto rend() const noexcept
  requires(flatstride)
  {
    return std::reverse_iterator(begin());
  }
  [[nodiscard]] constexpr auto front() const noexcept -> const T & {
    return *ptr;
  }
  [[nodiscard]] constexpr auto back() const noexcept -> const T & {
    if constexpr (flatstride) return *(end() - 1);
    else return ptr[sride(sz) * ptrdiff_t(row(sz)) - 1];
  }
  // indexing has two components:
  // 1. offsetting the pointer
  // 2. calculating new dim
  // static constexpr auto slice(Valid<T>, Index<S> auto i){
  //   auto
  // }
  [[gnu::flatten, gnu::always_inline]] constexpr auto
  operator[](Index<S> auto i) const noexcept -> decltype(auto) {
    return index<T>(ptr, sz, i);
  }
  // for (row/col)vectors, we drop the row/col, essentially broadcasting
  template <class R, class C>
  [[gnu::flatten, gnu::always_inline]] constexpr auto
  operator[](R r, C c) const noexcept -> decltype(auto) {
    return index<T>(ptr, sz, r, c);
  }
  [[nodiscard]] constexpr auto minRowCol() const -> ptrdiff_t {
    return std::min(ptrdiff_t(numRow()), ptrdiff_t(numCol()));
  }

  [[nodiscard]] constexpr auto diag() const noexcept {
    StridedRange r{minRowCol(), ptrdiff_t(RowStride(sz)) + 1};
    invariant(ptr != nullptr);
    return Array<T, StridedRange>{ptr, r};
  }
  [[nodiscard]] constexpr auto antiDiag() const noexcept {
    StridedRange r{minRowCol(), ptrdiff_t(RowStride(sz)) - 1};
    invariant(ptr != nullptr);
    return Array<T, StridedRange>{ptr + ptrdiff_t(Col(sz)) - 1, r};
  }
  [[nodiscard]] constexpr auto isSquare() const noexcept -> bool {
    return ptrdiff_t(Row(sz)) == ptrdiff_t(Col(sz));
  }
  [[nodiscard]] constexpr auto checkSquare() const -> Optional<ptrdiff_t> {
    ptrdiff_t N = ptrdiff_t(numRow());
    if (N != ptrdiff_t(numCol())) return {};
    return N;
  }

  [[nodiscard]] constexpr auto numRow() const noexcept {
    if constexpr (RowVectorDimension<S>) return Row<1>{};
    else return row(sz);
  }
  [[nodiscard]] constexpr auto numCol() const noexcept { return col(sz); }
  [[nodiscard]] constexpr auto rowStride() const noexcept {
    if constexpr (std::same_as<S, Length<>>) return RowStride<1>{};
    else return stride(sz);
  }
  [[nodiscard]] constexpr auto empty() const -> bool {
    if constexpr (StaticLength<S>) return S::staticint() == 0;
    else return sz == S{};
  }
  [[nodiscard]] constexpr auto size() const noexcept {
    if constexpr (StaticLength<S>) return S::staticint();
    else return ptrdiff_t(sz);
  }
  [[nodiscard]] constexpr auto dim() const noexcept -> S { return sz; }
  constexpr void clear()
  requires(std::same_as<S, Length<>>)
  {
    sz = S{};
  }
  [[nodiscard]] constexpr auto t() const { return Transpose{*this}; }
  [[nodiscard]] constexpr auto isExchangeMatrix() const -> bool
  requires(MatrixDimension<S>)
  {
    ptrdiff_t N = ptrdiff_t(numRow());
    if (N != ptrdiff_t(numCol())) return false;
    for (ptrdiff_t i = 0; i < N; ++i) {
      for (ptrdiff_t j = 0; j < N; ++j)
        if ((*this)[i, j] != (i + j == N - 1)) return false;
    }
  }
  [[nodiscard]] constexpr auto isDiagonal() const -> bool
  requires(MatrixDimension<S>)
  {
    for (ptrdiff_t r = 0; r < numRow(); ++r)
      for (ptrdiff_t c = 0; c < numCol(); ++c)
        if (r != c && (*this)[r, c] != 0) return false;
    return true;
  }
  [[nodiscard]] constexpr auto view() const noexcept -> Array<T, S> {
    invariant(ptr != nullptr);
    return Array<T, S>{ptr, this->sz};
  }

  [[nodiscard]] constexpr auto
  deleteCol(ptrdiff_t c) const -> ManagedArray<T, S> {
    static_assert(MatrixDimension<S>);
    auto new_dim = dim().similar(numRow() - 1);
    ManagedArray<T, decltype(new_dim)> A(new_dim);
    for (ptrdiff_t m = 0; m < numRow(); ++m) {
      A[m, _(0, c)] = (*this)[m, _(0, c)];
      A[m, _(c, math::end)] = (*this)[m, _(c + 1, math::end)];
    }
    return A;
  }
  [[nodiscard]] constexpr auto
  operator==(const Array &other) const noexcept -> bool {
    if constexpr (MatrixDimension<S> && !DenseLayout<S>) {
      ptrdiff_t M = ptrdiff_t(other.numRow());
      if ((numRow() != M) || (numCol() != other.numCol())) return false;
      // may not be dense, iterate over rows
      for (ptrdiff_t i = 0; i < M; ++i)
        if ((*this)[i, _] != other[i, _]) return false;
    } else {
      constexpr ptrdiff_t W = simd::Width<T>;
      ptrdiff_t N = size();
      if (N != other.size()) return false;
      if constexpr (W == 1) {
        for (ptrdiff_t i = 0; i < N; ++i)
          if ((*this)[i] != other[i]) return false;
      } else {
        for (ptrdiff_t i = 0;; i += W) {
          auto u{simd::index::unrollmask<1, W>(N, i)};
          if (!u) break;
          if (simd::cmp::ne<W, T>((*this)[u], other[u])) return false;
        }
      }
    }
    return true;
  }
  // FIXME: strided should skip over elements
  [[nodiscard]] constexpr auto norm2() const noexcept -> value_type {
    static_assert(DenseLayout<S>);
    value_type ret{0};
    for (auto x : *this) ret += x * x;
    return ret;
    // return std::transform_reduce(begin(), end(), begin(), 0.0);
  }
  // FIXME: strided should skips over elements
  [[nodiscard]] constexpr auto sum() const noexcept -> value_type {
    static_assert(VectorDimension<S> || DenseLayout<S>);
    value_type ret{0};
    for (auto x : *this) ret += x;
    return ret;
    // return std::reduce(begin(), end());
  }
  // interpret a bigger object as smaller
  template <typename U> [[nodiscard]] auto reinterpret() const {
    static_assert(sizeof(storage_type) % sizeof(U) == 0);
    static_assert(std::same_as<U, double>);
    if constexpr (std::same_as<U, T>) return *this;
    else {
      auto r = unwrapRow(numRow());
      auto c = unwrapCol(numCol());
      constexpr auto ratio = sizeof(storage_type) / sizeof(U);
#ifdef __cpp_lib_start_lifetime_as
      U *p = std::start_lifetime_as_array(reinterpret_cast<const U *>(data()),
                                          ratio * ptrdiff_t(rowStride()) * r);
#else
      const U *p = std::launder(reinterpret_cast<const U *>(data()));
#endif
      if constexpr (IsOne<decltype(r)>) {
        if constexpr (StaticInt<decltype(c)>)
          return Array<U, std::integral_constant<ptrdiff_t, c * ratio>>{p, {}};
        else return Array<U, ptrdiff_t>{p, c * ratio};
      } else if constexpr (DenseLayout<S>) {
        return Array<U, DenseDims<>>(p, DenseDims(row(r), col(c * ratio)));
      } else {
        ptrdiff_t stride = ptrdiff_t(rowStride()) * ratio;
        if constexpr (IsOne<decltype(c)>) {
          constexpr auto sr = std::integral_constant<ptrdiff_t, ratio>{};
          return Array<U, StridedDims<-1, ratio, -1>>(
            p, StridedDims(row(r), col(sr), rowStride(stride)));
        } else
          return Array<U, StridedDims<>>(
            p, StridedDims(row(r), col(c * ratio), rowStride(stride)));
      }
    }
  }
  friend void PrintTo(const Array &x, ::std::ostream *os) { *os << x; }
  [[nodiscard]] friend constexpr auto
  operator<=>(Array x, Array y) -> std::strong_ordering {
    ptrdiff_t M = x.size();
    ptrdiff_t N = y.size();
    for (ptrdiff_t i = 0, L = std::min(M, N); i < L; ++i)
      if (auto cmp = x[i] <=> y[i]; cmp != 0) return cmp;
    return M <=> N;
  };
  [[nodiscard]] constexpr auto eachRow() const -> SliceRange<const T, false>
  requires(MatrixDimension<S>)
  {
    return {data(), length(ptrdiff_t(Col(this->sz))), RowStride(this->sz),
            ptrdiff_t(Row(this->sz))};
  }
  [[nodiscard]] constexpr auto eachCol() const -> SliceRange<const T, true>
  requires(MatrixDimension<S>)
  {
    return {data(), length(ptrdiff_t(Row(this->sz))), RowStride(this->sz),
            ptrdiff_t(Col(this->sz))};
  }
  friend auto operator<<(std::ostream &os, Array x) -> std::ostream &
  requires(utils::Printable<T>)
  {
    if constexpr (MatrixDimension<S>)
      return printMatrix(os, Array<T, StridedDims<>>{x});
    else return utils::printVector(os, x.begin(), x.end());
  }
  [[nodiscard]] constexpr auto split(ptrdiff_t at) const
    -> containers::Pair<Array<T, Length<>>, Array<T, Length<>>>
  requires(VectorDimension<S>)
  {
    invariant(at <= size());
    return {(*this)[_(0, at)], (*this)[_(at, math::end)]};
  }
  [[nodiscard]] constexpr auto
  popFront() const -> containers::Pair<T, Array<T, Length<>>>
  requires(VectorDimension<S>)
  {
    invariant(0 < size());
    return {ptr[0], (*this)[_(1, math::end)]};
  }
#ifndef NDEBUG
  [[gnu::used]] void dump() const {
    // we can't combine `gnu::used` with `requires(utils::Printable<T>)`
    // requires(utils::Printable<T>)
    if constexpr (utils::Printable<T>)
      std::cout << "Size: " << ptrdiff_t(sz) << " " << *this << "\n";
  }
  // [[gnu::used]] void dump(const char *filename) const {
  //   if constexpr (std::integral<T>) {
  //     std::FILE *f = std::fopen(filename, "w");
  //     if (f == nullptr) return;
  //     (void)std::fprintf(f, "C= [");
  //     if constexpr (MatrixDimension<S>) {
  //       for (ptrdiff_t i = 0; i < Row(sz); ++i) {
  //         if (i) (void)std::fprintf(f, "\n");
  //         (void)std::fprintf(f, "%ld", long((*this)[i, 0]));
  //         for (ptrdiff_t j = 1; j < Col(sz); ++j)
  //           (void)std::fprintf(f, " %ld", long((*this)[i, j]));
  //       }
  //     } else {
  //       (void)std::fprintf(f, "%ld", long((*this)[0]));
  //       for (ptrdiff_t i = 1; (i < ptrdiff_t(sz)); ++i)
  //         (void)std::fprintf(f, ", %ld", long((*this)[i]));
  //     }
  //     (void)std::fprintf(f, "]");
  //     (void)std::fclose(f);
  //   }
  // }
#endif
protected:
  // NOLINTNEXTLINE(misc-non-private-member-variables-in-classes)
  const storage_type *ptr;
  // NOLINTNEXTLINE(misc-non-private-member-variables-in-classes)
  [[no_unique_address]] S sz;
};

static_assert(
  std::is_trivially_default_constructible_v<Array<int64_t, DenseDims<>>>);

template <class T, Dimension S, bool Compress>
struct POLY_MATH_GSL_POINTER MutArray
  : Array<T, S, Compress>,
    ArrayOps<T, S, MutArray<T, S, Compress>> {
  using BaseT = Array<T, S, Compress>;
  // using BaseT::BaseT;
  using BaseT::operator[], BaseT::data, BaseT::begin, BaseT::end, BaseT::rbegin,
    BaseT::rend, BaseT::front, BaseT::back;
  using storage_type = typename BaseT::storage_type;

  explicit constexpr MutArray() = default;
  constexpr MutArray(const MutArray &) = default;
  constexpr MutArray(MutArray &&) noexcept = default;
  // constexpr auto operator=(const MutArray &) -> MutArray & = delete;
  constexpr auto operator=(const MutArray &) -> MutArray & = default;
  constexpr auto operator=(MutArray &&) noexcept -> MutArray & = default;
  constexpr MutArray(storage_type *p, S s) : BaseT(p, s) {}

  constexpr void truncate(S nz) {
    S oz = this->sz;
    this->sz = nz;
    if constexpr (std::same_as<S, Length<>>) {
      invariant(ptrdiff_t(nz) <= ptrdiff_t(oz));
    } else if constexpr (std::convertible_to<S, DenseDims<>>) {
      auto new_x = ptrdiff_t{RowStride(nz)}, old_x = ptrdiff_t{RowStride(oz)},
           new_n = ptrdiff_t{Col(nz)}, old_n = ptrdiff_t{Col(oz)},
           new_m = ptrdiff_t{Row(nz)}, old_m = ptrdiff_t{Row(oz)};
      invariant(new_m <= old_m);
      invariant(new_n <= old_n);
      invariant(new_x <= old_x);
      ptrdiff_t cols_to_copy = new_n, rows_to_copy = new_m;
      // we only need to copy if memory shifts position
      bool copy_cols = ((cols_to_copy > 0) && (new_x != old_x));
      // if we're in place, we have 1 less row to copy
      if (((--rows_to_copy) > 0) && (copy_cols)) {
        // truncation, we need to copy rows to increase stride
        storage_type *src = data(), *dst = src;
        do {
          src += old_x;
          dst += new_x;
          std::copy_n(src, cols_to_copy, dst);
        } while (--rows_to_copy);
      }
    } else {
      static_assert(MatrixDimension<S>, "Can only resize 1 or 2d containers.");
      invariant(nz.row() <= oz.row());
      invariant(nz.col() <= oz.col());
    }
  }

  constexpr void truncate(Row<> r) {
    if constexpr (std::same_as<S, Length<>>) {
      return truncate(S(r));
    } else if constexpr (std::convertible_to<S, DenseDims<>>) {
      static_assert(!std::convertible_to<S, SquareDims<>>,
                    "if truncating a row, matrix must be strided or dense.");
      invariant(r <= Row(this->sz));
      DenseDims new_sz = this->sz;
      truncate(new_sz.set(r));
    } else {
      static_assert(std::convertible_to<S, StridedDims<>>);
      invariant(r <= Row(this->sz));
      this->sz.set(r);
    }
  }
  constexpr void truncate(Col<> c) {
    if constexpr (std::same_as<S, Length<>>) {
      return truncate(S(c));
    } else if constexpr (std::is_same_v<S, DenseDims<>>) {
      static_assert(!std::convertible_to<S, SquareDims<>>,
                    "if truncating a col, matrix must be strided or dense.");
      invariant(c <= Col(this->sz));
      DenseDims new_sz = this->sz;
      truncate(new_sz.set(c));
    } else {
      static_assert(std::convertible_to<S, StridedDims<>>);
      invariant(c <= Col(this->sz));
      this->sz.set(c);
    }
  }

  template <class... Args>
  constexpr MutArray(Args &&...args)
    : Array<T, S, Compress>(std::forward<Args>(args)...) {}

  template <std::convertible_to<T> U, std::convertible_to<S> V>
  constexpr MutArray(Array<U, V> a) : Array<T, S>(a) {}
  template <size_t N>
  constexpr MutArray(std::array<T, N> &a) : Array<T, S>(a.data(), length(N)) {}
  [[nodiscard]] constexpr auto data() noexcept -> storage_type * {
    invariant(this->ptr != nullptr || ptrdiff_t(this->sz) == 0);
    return const_cast<storage_type *>(this->ptr);
  }
  [[nodiscard]] constexpr auto wrappedPtr() noexcept -> Valid<T> {
    return data();
  }

  [[nodiscard]] constexpr auto begin() noexcept -> StridedIterator<T>
  requires(std::is_same_v<S, StridedRange>)
  {
    return StridedIterator{const_cast<storage_type *>(this->ptr),
                           this->sz.stride_};
  }
  [[nodiscard]] constexpr auto
  begin() noexcept -> storage_type *requires(!std::is_same_v<S, StridedRange>) {
    return const_cast<storage_type *>(this->ptr);
  }

  [[nodiscard]] constexpr auto end() noexcept {
    return begin() + ptrdiff_t(this->sz);
  }
  [[nodiscard]] constexpr auto rbegin() noexcept {
    return std::reverse_iterator(end());
  }
  [[nodiscard]] constexpr auto rend() noexcept {
    return std::reverse_iterator(begin());
  }
  constexpr auto front() noexcept -> T & { return *begin(); }
  constexpr auto back() noexcept -> T & { return *(end() - 1); }
  [[gnu::flatten, gnu::always_inline]] constexpr auto
  operator[](Index<S> auto i) noexcept -> decltype(auto) {
    return index<T>(data(), this->sz, i);
  }
  template <class R, class C>
  [[gnu::flatten, gnu::always_inline]] constexpr auto
  operator[](R r, C c) noexcept -> decltype(auto) {
    return index<T>(data(), this->sz, r, c);
  }
  constexpr void fill(T value) {
    std::fill_n(this->data(), ptrdiff_t(this->dim()), value);
  }
  [[nodiscard]] constexpr auto diag() noexcept {
    Length l =
      length(std::min(ptrdiff_t(Row(this->sz)), ptrdiff_t(Col(this->sz))));
    RowStride rs = RowStride(this->sz);
    StridedRange r{l, ++rs};
    return MutArray<T, StridedRange>{data(), r};
  }
  [[nodiscard]] constexpr auto antiDiag() noexcept {
    Col<> c = Col(this->sz);
    Length l =
      length(ptrdiff_t(std::min(ptrdiff_t(Row(this->sz)), ptrdiff_t(c))));
    RowStride rs = RowStride(this->sz);
    StridedRange r{l, --rs};
    return MutArray<T, StridedRange>{data() + ptrdiff_t(c) - 1, r};
  }
  constexpr void erase(ptrdiff_t i)
  requires(std::same_as<S, Length<>>)
  {
    auto old_len = ptrdiff_t(this->sz--);
    if (i < this->sz)
      std::copy(this->data() + i + 1, data() + old_len, this->data() + i);
  }
  constexpr void erase(Row<> r) {
    if constexpr (std::same_as<S, Length<>>) {
      return erase(S(r));
    } else if constexpr (std::convertible_to<S, DenseDims<>>) {
      static_assert(!std::convertible_to<S, SquareDims<>>,
                    "if erasing a row, matrix must be strided or dense.");
      auto col = ptrdiff_t{Col(this->sz)},
           new_row = ptrdiff_t{Row(this->sz)} - 1;
      this->sz.set(row(new_row));
      if ((col == 0) || (r == new_row)) return;
      storage_type *dst = data() + ptrdiff_t(r) * col;
      std::copy_n(dst + col, (new_row - ptrdiff_t(r)) * col, dst);
    } else {
      static_assert(std::convertible_to<S, StridedDims<>>);
      auto stride = ptrdiff_t{RowStride(this->sz)},
           col = ptrdiff_t{Col(this->sz)},
           new_row = ptrdiff_t{Row(this->sz)} - 1;
      this->sz.set(row(new_row));
      if ((col == 0) || (r == new_row)) return;
      invariant(col <= stride);
      if ((col + (512 / (sizeof(T)))) <= stride) {
        storage_type *dst = data() + ptrdiff_t(r) * stride;
        for (auto m = ptrdiff_t(r); m < new_row; ++m) {
          storage_type *src = dst + stride;
          std::copy_n(src, col, dst);
          dst = src;
        }
      } else {
        storage_type *dst = data() + ptrdiff_t(r) * stride;
        std::copy_n(dst + stride, (new_row - ptrdiff_t(r)) * stride, dst);
      }
    }
  }
  constexpr void erase(Col<> c) {
    if constexpr (std::same_as<S, Length<>>) {
      return erase(S(c));
    } else if constexpr (std::convertible_to<S, DenseDims<>>) {
      static_assert(!std::convertible_to<S, SquareDims<>>,
                    "if erasing a col, matrix must be strided or dense.");
      auto new_col = ptrdiff_t{Col(this->sz)}, old_col = new_col--,
           row = ptrdiff_t{Row(this->sz)};
      this->sz.set(col(new_col));
      ptrdiff_t cols_to_copy = new_col - ptrdiff_t(c);
      if ((cols_to_copy == 0) || (row == 0)) return;
      // we only need to copy if memory shifts position
      for (ptrdiff_t m = 0; m < row; ++m) {
        storage_type *dst = data() + m * new_col + ptrdiff_t(c);
        storage_type *src = data() + m * old_col + ptrdiff_t(c) + 1;
        std::copy_n(src, cols_to_copy, dst);
      }
    } else {
      static_assert(std::convertible_to<S, StridedDims<>>);
      auto stride = ptrdiff_t{RowStride(this->sz)},
           new_col = ptrdiff_t{Col(this->sz)} - 1,
           row = ptrdiff_t{Row(this->sz)};
      this->sz.set(col(new_col));
      ptrdiff_t cols_to_copy = new_col - ptrdiff_t(c);
      if ((cols_to_copy == 0) || (row == 0)) return;
      // we only need to copy if memory shifts position
      for (ptrdiff_t m = 0; m < row; ++m) {
        storage_type *dst = data() + m * stride + ptrdiff_t(c);
        std::copy_n(dst + 1, cols_to_copy, dst);
      }
    }
  }
  constexpr void moveLast(Col<> j) {
    static_assert(MatrixDimension<S>);
    if (j == this->numCol()) return;
    Col Nd = this->numCol() - 1;
    for (ptrdiff_t m = 0; m < this->numRow(); ++m) {
      auto x = (*this)[m, ptrdiff_t(j)];
      for (auto n = ptrdiff_t(j); n < Nd;) {
        ptrdiff_t o = n++;
        (*this)[m, o] = (*this)[m, n];
      }
      (*this)[m, Nd] = x;
    }
  }
  constexpr auto eachRow() -> SliceRange<T, false>
  requires(MatrixDimension<S>)
  {
    return {data(), length(ptrdiff_t(Col(this->sz))), RowStride(this->sz),
            ptrdiff_t(Row(this->sz))};
  }
  constexpr auto eachCol() -> SliceRange<T, true>
  requires(MatrixDimension<S>)
  {
    return {data(), length(ptrdiff_t(Row(this->sz))), RowStride(this->sz),
            ptrdiff_t(Col(this->sz))};
  }
  template <typename U> [[nodiscard]] auto reinterpret() {
    static_assert(sizeof(storage_type) % sizeof(U) == 0);
    static_assert(std::same_as<U, double>);
    if constexpr (std::same_as<U, T>) return *this;
    else {
      auto r = unwrapRow(this->numRow());
      auto c = unwrapCol(this->numCol());
      constexpr size_t ratio = sizeof(storage_type) / sizeof(U);
#ifdef __cpp_lib_start_lifetime_as
      U *p = std::start_lifetime_as_array(reinterpret_cast<const U *>(data()),
                                          ratio * ptrdiff_t(rowStride()) * r);
#else
      U *p = std::launder(reinterpret_cast<U *>(data()));
#endif
      if constexpr (IsOne<decltype(r)>) {
        if constexpr (StaticInt<decltype(c)>)
          return MutArray<U, std::integral_constant<ptrdiff_t, c * ratio>>(
            p, std::integral_constant<ptrdiff_t, c * ratio>{});
        else return MutArray<U, ptrdiff_t>{p, c * ratio};
      } else if constexpr (DenseLayout<S>) {
        return MutArray<U, DenseDims<>>(p, DenseDims(row(r), col(c * ratio)));
      } else {
        ptrdiff_t stride = ptrdiff_t(this->rowStride()) * ratio;
        if constexpr (IsOne<decltype(c)>) {
          constexpr auto sr = std::integral_constant<ptrdiff_t, ratio>{};
          return MutArray<U, StridedDims<-1, ratio, -1>>(
            p, StridedDims(row(r), col(sr), rowStride(stride)));
        } else
          return MutArray<U, StridedDims<>>(
            p, StridedDims(row(r), col(c * ratio), rowStride(stride)));
      }
    }
  }
};

template <typename T, typename S>
Array(T *, S) -> Array<utils::decompressed_t<T>, S>;
template <typename T, typename S>
MutArray(T *, S) -> MutArray<utils::decompressed_t<T>, S>;

template <typename T, typename S> MutArray(MutArray<T, S>) -> MutArray<T, S>;

static_assert(std::convertible_to<Array<int64_t, SquareDims<>>,
                                  Array<int64_t, DenseDims<>>>);
static_assert(std::convertible_to<Array<int64_t, DenseDims<8, 8>>,
                                  Array<int64_t, DenseDims<>>>);
static_assert(std::convertible_to<Array<int64_t, SquareDims<>>,
                                  Array<int64_t, StridedDims<>>>);
static_assert(std::convertible_to<Array<int64_t, DenseDims<>>,
                                  Array<int64_t, StridedDims<>>>);
static_assert(std::convertible_to<MutArray<int64_t, SquareDims<>>,
                                  Array<int64_t, DenseDims<>>>);
static_assert(std::convertible_to<MutArray<int64_t, SquareDims<>>,
                                  Array<int64_t, StridedDims<>>>);
static_assert(std::convertible_to<MutArray<int64_t, DenseDims<>>,
                                  Array<int64_t, StridedDims<>>>);
static_assert(std::convertible_to<MutArray<int64_t, SquareDims<>>,
                                  MutArray<int64_t, DenseDims<>>>);
static_assert(std::convertible_to<MutArray<int64_t, SquareDims<>>,
                                  MutArray<int64_t, StridedDims<>>>);
static_assert(std::convertible_to<MutArray<int64_t, DenseDims<>>,
                                  MutArray<int64_t, StridedDims<>>>);
static_assert(AbstractVector<Array<int64_t, Length<>>>);
static_assert(!AbstractVector<Array<int64_t, StridedDims<>>>);
static_assert(AbstractMatrix<Array<int64_t, StridedDims<>>>);
static_assert(RowVector<Array<int64_t, Length<>>>);
static_assert(ColVector<Transpose<Array<int64_t, Length<>>>>);
static_assert(ColVector<Array<int64_t, StridedRange>>);
static_assert(RowVector<Transpose<Array<int64_t, StridedRange>>>);

// template <typename T, bool Column>
// inline constexpr auto SliceIterator<T, Column>::operator*()
//   -> SliceIterator<T, Column>::value_type {
//   if constexpr (Column) return {data + idx, StridedRange{len, rowStride}};
//   else return {data + rowStride * idx, len};
// }
template <typename T, bool Column>
constexpr auto SliceIterator<T, Column>::operator*() const
  -> SliceIterator<T, Column>::value_type {
  if constexpr (Column) return {data_ + idx_, StridedRange{len_, row_stride_}};
  else return {data_ + ptrdiff_t(row_stride_) * idx_, len_};
}
// template <typename T, bool Column>
// inline constexpr auto operator*(SliceIterator<T, Column> it)
//   -> SliceIterator<T, Column>::value_type {
//   if constexpr (Column)
//     return {it.data + it.idx, StridedRange{it.len, it.rowStride}};
//   else return {it.data + it.rowStride * it.idx, it.len};
// }

static_assert(std::weakly_incrementable<SliceIterator<int64_t, false>>);
static_assert(std::forward_iterator<SliceIterator<int64_t, false>>);
static_assert(std::ranges::forward_range<SliceRange<int64_t, false>>);
static_assert(std::ranges::range<SliceRange<int64_t, false>>);
using ITEST = std::iter_rvalue_reference_t<SliceIterator<int64_t, true>>;
static_assert(std::is_same_v<ITEST, MutArray<int64_t, StridedRange, false>>);

/// Non-owning view of a managed array, capable of resizing,
/// but not of re-allocating in case the capacity is exceeded.
template <class T, Dimension S>
struct POLY_MATH_GSL_POINTER ResizeableView : MutArray<T, S> {
  using BaseT = MutArray<T, S>;
  using U = containers::default_capacity_type_t<S>;
  using storage_type = typename BaseT::storage_type;
  constexpr ResizeableView() noexcept : BaseT(nullptr, S{}), capacity_(U{}) {}
  constexpr ResizeableView(storage_type *p, S s, U c) noexcept
    : BaseT(p, s), capacity_(c) {}
  constexpr ResizeableView(alloc::Arena<> *a, S s, U c) noexcept
    : ResizeableView{a->template allocate<storage_type>(ptrdiff_t(c)), s, c} {}

  [[nodiscard]] constexpr auto isFull() const -> bool {
    return ptrdiff_t(this->sz) == capacity_;
  }

  template <class... Args>
  constexpr auto emplace_back_within_capacity(Args &&...args) -> decltype(auto)
  requires(std::same_as<S, Length<>>)
  {
    auto s = ptrdiff_t(this->sz), c = ptrdiff_t(this->capacity_);
    invariant(s < c);
    T &ret =
      *std::construct_at<T>(this->data() + s, std::forward<Args>(args)...);
    this->sz = length(s + 1z);
    return ret;
  }
  /// Allocates extra space if needed
  /// Has a different name to make sure we avoid ambiguities.
  template <class... Args>
  constexpr auto emplace_backa(alloc::Arena<> *alloc,
                               Args &&...args) -> decltype(auto)
  requires(std::same_as<S, Length<>>)
  {
    if (isFull()) reserve(alloc, (ptrdiff_t(capacity_) + 1z) * 2z);
    return *std::construct_at(this->data() + ptrdiff_t(this->sz++),
                              std::forward<Args>(args)...);
  }
  constexpr void push_back_within_capacity(T value)
  requires(std::same_as<S, Length<>>)
  {
    auto s = ptrdiff_t(this->sz), c = ptrdiff_t(this->capacity_);
    invariant(s < c);
    std::construct_at<T>(this->data() + s, std::move(value));
    this->sz = length(s + 1z);
  }
  constexpr void push_back(alloc::Arena<> *alloc, T value)
  requires(std::same_as<S, Length<>>)
  {
    if (isFull()) reserve(alloc, (ptrdiff_t(capacity_) + 1z) * 2z);
    std::construct_at<T>(this->data() + ptrdiff_t(this->sz++),
                         std::move(value));
  }
  constexpr void pop_back()
  requires(std::same_as<S, Length<>>)
  {
    invariant(this->sz > 0);
    if constexpr (std::is_trivially_destructible_v<T>) --this->sz;
    else std::destroy_at(this->data() + ptrdiff_t(--this->sz));
  }
  constexpr auto pop_back_val() -> T
  requires(std::same_as<S, Length<>>)
  {
    invariant(this->sz > 0);
    return std::move(this->data()[ptrdiff_t(--this->sz)]);
  }
  constexpr void resize(ptrdiff_t nz)
  requires(std::same_as<S, Length<>>)
  {
    resize(length(nz));
  }
  // behavior
  // if S is StridedDims, then we copy data.
  // If the new dims are larger in rows or cols, we fill with 0.
  // If the new dims are smaller in rows or cols, we truncate.
  // New memory outside of dims (i.e., stride larger), we leave uninitialized.
  //
  constexpr void resize(S nz) {
    S oz = this->sz;
    this->sz = nz;
    if constexpr (std::same_as<S, Length<>>) {
      auto ozs = ptrdiff_t(oz), nzs = ptrdiff_t(nz);
      invariant(nzs <= capacity_);
      if constexpr (!std::is_trivially_destructible_v<T>) {
        if (nzs < ozs) std::destroy_n(this->data() + nz, oz - nz);
        else if (nzs > ozs)
          for (ptrdiff_t i = ozs; i < nzs; ++i)
            std::construct_at(this->data() + i);
      } else if (nz > oz)
        std::fill(this->data() + ozs, this->data() + nzs, T{});
    } else {
      static_assert(std::is_trivially_destructible_v<T>,
                    "Resizing matrices holding non-is_trivially_destructible_v "
                    "objects is not yet supported.");
      static_assert(MatrixDimension<S>, "Can only resize 1 or 2d containers.");
      auto new_x = ptrdiff_t{RowStride(nz)}, old_x = ptrdiff_t{RowStride(oz)},
           new_n = ptrdiff_t{Col(nz)}, old_n = ptrdiff_t{Col(oz)},
           new_m = ptrdiff_t{Row(nz)}, old_m = ptrdiff_t{Row(oz)};
      invariant(U(nz) <= capacity_);
      U len = U(nz);
      T *npt = this->data();
      // we can copy forward so long as the new stride is smaller
      // so that the start of the dst range is outside of the src range
      // we can also safely forward copy if we allocated a new ptr
      bool forward_copy = (new_x <= old_x);
      ptrdiff_t cols_to_copy = std::min(old_n, new_n);
      // we only need to copy if memory shifts position
      bool copy_cols = ((cols_to_copy > 0) && (new_x != old_x));
      // if we're in place, we have 1 less row to copy
      ptrdiff_t rows_to_copy = std::min(old_m, new_m) - 1;
      ptrdiff_t fill_count = new_n - cols_to_copy;
      if ((rows_to_copy) && (copy_cols || fill_count)) {
        if (forward_copy) {
          // truncation, we need to copy rows to increase stride
          T *src = this->data() + old_x;
          T *dst = npt + new_x;
          do {
            if (copy_cols) std::copy_n(src, cols_to_copy, dst);
            if (fill_count) std::fill_n(dst + cols_to_copy, fill_count, T{});
            src += old_x;
            dst += new_x;
          } while (--rows_to_copy);
        } else /* [[unlikely]] */ {
          // backwards copy, only needed when we increasing stride but not
          // reallocating, which should be comparatively uncommon.
          // Should probably benchmark or determine actual frequency
          // before adding `[[unlikely]]`.
          T *src = this->data() + (rows_to_copy + 1) * old_x;
          T *dst = npt + (rows_to_copy + 1) * new_x;
          do {
            src -= old_x;
            dst -= new_x;
            if (cols_to_copy)
              std::copy_backward(src, src + cols_to_copy, dst + cols_to_copy);
            if (fill_count) std::fill_n(dst + cols_to_copy, fill_count, T{});
          } while (--rows_to_copy);
        }
      }
      // zero init remaining rows
      for (ptrdiff_t m = old_m; m < new_m; ++m)
        std::fill_n(npt + m * new_x, new_n, T{});
    }
  }

  constexpr void resize(Row<> r) {
    if constexpr (std::same_as<S, Length<>>) {
      return resize(S(r));
    } else if constexpr (MatrixDimension<S>) {
      S nz = this->sz;
      return resize(nz.set(r));
    }
  }
  constexpr void resize(Col<> c) {
    if constexpr (std::same_as<S, Length<>>) {
      return resize(S(c));
    } else if constexpr (MatrixDimension<S>) {
      S nz = this->sz;
      return resize(nz.set(c));
    }
  }
  constexpr void resizeForOverwrite(S M) {
    invariant(ptrdiff_t(M) <= ptrdiff_t(this->sz));
    if constexpr (!std::is_trivially_destructible_v<T>) {
      ptrdiff_t nz = ptrdiff_t(M), oz = ptrdiff_t(this->sz);
      if (nz < oz) std::destroy_n(this->data() + nz, oz - nz);
      else if (nz > oz) // FIXME: user should initialize?
        for (ptrdiff_t i = oz; i < nz; ++i) std::construct_at(this->data() + i);
    }
    this->sz = M;
  }
  constexpr void resizeForOverwrite(ptrdiff_t M)
  requires(std::same_as<S, Length<>>)
  {
    resizeForOverwrite(length(M));
  }
  constexpr void resizeForOverwrite(Row<> r) {
    if constexpr (std::same_as<S, Length<>>) {
      return resizeForOverwrite(S(r));
    } else if constexpr (MatrixDimension<S>) {
      S nz = this->sz;
      return resizeForOverwrite(nz.set(r));
    }
  }
  constexpr void resizeForOverwrite(Col<> c) {
    if constexpr (std::same_as<S, Length<>>) {
      return resizeForOverwrite(S(c));
    } else if constexpr (MatrixDimension<S>) {
      S nz = this->sz;
      return resizeForOverwrite(nz.set(c));
    }
  }
  [[nodiscard]] constexpr auto getCapacity() const -> U { return capacity_; }

  // set size and 0.
  constexpr void setSize(Row<> r, Col<> c) {
    resizeForOverwrite({r, c});
    this->fill(0);
  }
  constexpr void resize(Row<> MM, Col<> NN) { resize(DenseDims{MM, NN}); }
  constexpr void resizeForOverwrite(Row<> M, Col<> N, RowStride<> X) {
    invariant(X >= N);
    if constexpr (std::is_same_v<S, StridedDims<>>)
      resizeForOverwrite(S{M, N, X});
    else if constexpr (std::is_same_v<S, SquareDims<>>) {
      invariant(ptrdiff_t(M) == ptrdiff_t(N));
      resizeForOverwrite(S{M});
    } else {
      static_assert(std::is_same_v<S, DenseDims<>>);
      resizeForOverwrite(S{M, N});
    }
  }
  constexpr void resizeForOverwrite(Row<> M, Col<> N) {
    if constexpr (std::is_same_v<S, StridedDims<>>)
      resizeForOverwrite(S{M, N, {ptrdiff_t(N)}});
    else if constexpr (std::is_same_v<S, SquareDims<>>) {
      invariant(ptrdiff_t(M) == ptrdiff_t(N));
      resizeForOverwrite(S{M});
    } else resizeForOverwrite(S{M, N});
  }

  constexpr auto
  insert_within_capacity(T *p, T x) -> T *requires(std::same_as<S, Length<>>) {
    invariant(p >= this->data());
    T *e = this->data() + ptrdiff_t(this->sz);
    invariant(p <= e);
    invariant(this->sz < capacity_);
    if constexpr (BaseT::trivial) {
      if (p < e) std::copy_backward(p, e, e + 1);
      *p = std::move(x);
    } else if (p < e) {
      std::construct_at<T>(e, std::move(*(e - 1)));
      for (; --e != p;) *e = std::move(*(e - 1));
      *p = std::move(x);
    } else std::construct_at<T>(e, std::move(x));
    ++this->sz;
    return p;
  }

  template <size_t SlabSize, bool BumpUp>
  constexpr void reserve(alloc::Arena<SlabSize, BumpUp> *alloc,
                         ptrdiff_t newCapacity) {
    if (newCapacity <= capacity_) return;
    this->ptr = alloc->template reallocate<false, T>(
      const_cast<T *>(this->ptr), ptrdiff_t(capacity_), newCapacity,
      ptrdiff_t(this->sz));
    // T *oldPtr =
    //   std::exchange(this->data(), alloc->template
    //   allocate<T>(newCapacity));
    // std::copy_n(oldPtr, U(this->sz), this->data());
    // alloc->deallocate(oldPtr, capacity);
    capacity_ = capacity(newCapacity);
  }

protected:
  // NOLINTNEXTLINE(misc-non-private-member-variables-in-classes)
  [[no_unique_address]] U capacity_{0};
};

static_assert(std::is_copy_assignable_v<Array<void *, Length<>>>);
static_assert(std::is_copy_assignable_v<MutArray<void *, Length<>>>);
static_assert(std::is_trivially_copyable_v<MutArray<void *, Length<>>>);
static_assert(std::is_trivially_move_assignable_v<MutArray<void *, Length<>>>);
[[nodiscard]] constexpr auto newCapacity(ptrdiff_t c) -> ptrdiff_t {
  return c ? c + c : 4z;
}

template <class T, class S>
concept AbstractSimilar = (MatrixDimension<S> && AbstractMatrix<T>) ||
                          (VectorDimension<S> && AbstractVector<T>);

/// Stores memory, then pointer.
/// Thus struct's alignment determines initial alignment
/// of the stack memory.
/// Information related to size is then grouped next to the pointer.
///
/// The Intel compiler + OpenMP appears to memcpy data around,
/// or at least build ManagedArrays bypassing the constructors listed here.
/// This caused invalid frees, as the pointer still pointed to the old
/// stack memory.
template <class T, Dimension S, ptrdiff_t StackStorage, class A>
struct POLY_MATH_GSL_OWNER ManagedArray : ResizeableView<T, S> {
  // static_assert(std::is_trivially_destructible_v<T>);
  using BaseT = ResizeableView<T, S>;
  using U = containers::default_capacity_type_t<S>;
  using storage_type = typename BaseT::storage_type;
  // We're deliberately not initializing storage.
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wuninitialized"
#else
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
#endif
  constexpr ManagedArray(A a) noexcept
    : BaseT{memory_.data(), S{}, capacity(StackStorage)}, allocator_{a} {
#ifndef NDEBUG
    if (!StackStorage) return;
    if constexpr (std::numeric_limits<T>::has_signaling_NaN)
      std::fill_n(this->data(), StackStorage,
                  std::numeric_limits<T>::signaling_NaN());
    else if constexpr (std::numeric_limits<T>::is_specialized)
      std::fill_n(this->data(), StackStorage, std::numeric_limits<T>::min());
#endif
  }
  constexpr ManagedArray(S s, A a) noexcept
    : BaseT{memory_.data(), s, capacity(StackStorage)}, allocator_{a} {
    U len = U(capacity(ptrdiff_t(this->sz)));
    if (len > StackStorage) this->allocateAtLeast(len);
#ifndef NDEBUG
    if (!len) return;
    auto l = ptrdiff_t(len);
    if constexpr (std::numeric_limits<T>::has_signaling_NaN)
      std::fill_n(this->data(), l, std::numeric_limits<T>::signaling_NaN());
    else if constexpr (std::numeric_limits<T>::is_specialized)
      std::fill_n(this->data(), l, std::numeric_limits<T>::min());
#endif
  }
  constexpr ManagedArray(S s, T x, A a) noexcept
    : BaseT{memory_.data(), s, capacity(StackStorage)}, allocator_{a} {
    auto len = ptrdiff_t(this->sz);
    if (len > StackStorage) this->allocateAtLeast(capacity(len));
    if (len) std::fill_n(this->data(), len, x);
  }
  constexpr ManagedArray() noexcept : ManagedArray(A{}) {};
  constexpr ManagedArray(S s) noexcept : ManagedArray(s, A{}) {};
  constexpr ManagedArray(ptrdiff_t s) noexcept
  requires(std::same_as<S, SquareDims<>>)
    : ManagedArray(SquareDims<>{row(s)}, A{}) {};
  constexpr ManagedArray(S s, T x) noexcept : ManagedArray(s, x, A{}) {};

  // constexpr ManagedArray(std::type_identity<T>) noexcept :
  // ManagedArray(A{}){};
  constexpr ManagedArray(std::type_identity<T>, S s) noexcept
    : ManagedArray(s, A{}) {};
  // constexpr ManagedArray(std::type_identity<T>, A a) noexcept :
  // ManagedArray(a){}; constexpr ManagedArray(std::type_identity<T>, S s, A
  // a) noexcept
  //   : ManagedArray(s, a){};
  constexpr ManagedArray(T x) noexcept
  requires(std::same_as<S, Length<>>)
    : BaseT{memory_.data(), S{}, capacity(StackStorage)}, allocator_(A{}) {
    if constexpr (StackStorage == 0) this->growUndef(1);
    this->push_back_within_capacity(std::move(x));
  }

  template <class D>
  constexpr ManagedArray(const ManagedArray<T, D, StackStorage, A> &b) noexcept
    : BaseT{memory_.data(), S(b.dim()), capacity(StackStorage)},
      allocator_(b.get_allocator()) {
    auto len = ptrdiff_t(this->sz);
    this->growUndef(len);
    std::copy_n(b.data(), len, this->data());
  }
  template <std::convertible_to<T> Y, class D, class AY>
  constexpr ManagedArray(const ManagedArray<Y, D, StackStorage, AY> &b) noexcept
    : BaseT{memory_.data(), S{}, capacity(StackStorage)},
      allocator_(b.get_allocator()) {
    S d = b.dim();
    auto len = ptrdiff_t(d);
    this->growUndef(len);
    this->sz = d;
    (*this) << b;
  }
  template <std::convertible_to<T> Y, size_t M>
  constexpr ManagedArray(std::array<Y, M> il) noexcept
    : BaseT{memory_.data(), S{}, capacity(StackStorage)} {
    auto len = ptrdiff_t(M);
    this->growUndef(len);
    std::copy_n(il.begin(), len, this->data());
    this->sz = math::length(ptrdiff_t(M));
  }
  template <std::convertible_to<T> Y, class D, class AY>
  constexpr ManagedArray(const ManagedArray<Y, D, StackStorage, AY> &b,
                         S s) noexcept
    : BaseT{memory_.data(), S(s), capacity(StackStorage)},
      allocator_(b.get_allocator()) {
    auto len = ptrdiff_t(this->sz);
    invariant(len == U(b.size()));
    this->growUndef(len);
    T *p = this->data();
    std::copy_n(b.data(), len, this->data());
  }
  constexpr ManagedArray(const ManagedArray &b) noexcept
    : BaseT{memory_.data(), S(b.dim()), capacity(StackStorage)},
      allocator_(b.get_allocator()) {
    auto len = ptrdiff_t(this->sz);
    this->growUndef(len);
    std::copy_n(b.data(), len, this->data());
  }
  constexpr ManagedArray(const Array<T, S> &b) noexcept
    : BaseT{memory_.data(), S(b.dim()), capacity(StackStorage)} {
    auto len = ptrdiff_t(this->sz);
    this->growUndef(len);
    std::copy_n(b.data(), len, this->data());
  }
  template <AbstractSimilar<S> V>
  constexpr ManagedArray(const V &b) noexcept
    : BaseT{memory_.data(), S(shape(b)), capacity(StackStorage)} {
    this->growUndef(ptrdiff_t(this->sz));
    (*this) << b;
  }
  template <class D>
  constexpr ManagedArray(ManagedArray<T, D, StackStorage, A> &&b) noexcept
    : BaseT{memory_.data(), b.dim(), U(capacity(StackStorage))},
      allocator_(b.get_allocator()) {
    if (!b.isSmall()) { // steal
      this->ptr = b.data();
      this->capacity_ = b.getCapacity();
    } else std::copy_n(b.data(), ptrdiff_t(b.dim()), this->data());
    b.resetNoFree();
  }
  constexpr ManagedArray(ManagedArray &&b) noexcept
    : BaseT{memory_.data(), b.dim(), U(capacity(StackStorage))},
      allocator_(b.get_allocator()) {
    if constexpr (StackStorage) {
      if (!b.isSmall()) { // steal
        this->ptr = b.data();
        this->capacity_ = b.getCapacity();
      } else std::copy_n(b.data(), ptrdiff_t(b.dim()), this->data());
    } else {
      this->ptr = b.ptr;
      this->capacity_ = b.getCapacity();
    }
    b.resetNoFree();
  }
  template <class D>
  constexpr ManagedArray(ManagedArray<T, D, StackStorage, A> &&b, S s) noexcept
    : BaseT{memory_.data(), s, U(capacity(StackStorage))},
      allocator_(b.get_allocator()) {
    if (!b.isSmall()) { // steal
      this->ptr = b.data();
      this->capacity_ = b.getCapacity();
    } else std::copy_n(b.data(), ptrdiff_t(b.dim()), this->data());
    b.resetNoFree();
  }
  template <std::convertible_to<T> Y>
  constexpr ManagedArray(const SmallSparseMatrix<Y> &B)
    : BaseT{memory_.data(), B.dim(), capacity(StackStorage)} {
    auto len = ptrdiff_t(this->sz);
    this->growUndef(len);
    this->fill(0);
    ptrdiff_t k = 0;
    for (ptrdiff_t i = 0; i < this->numRow(); ++i) {
      uint32_t m = B.getRows()[i] & 0x00ffffff;
      ptrdiff_t j = 0;
      while (m) {
        uint32_t tz = std::countr_zero(m);
        m >>= tz + 1;
        j += tz;
        (*this)[i, j++] = T(B.getNonZeros()[k++]);
      }
    }
    invariant(k == B.getNonZeros().size());
  }
  constexpr ManagedArray(const ColVector auto &v)
  requires(MatrixDimension<S>)
    : BaseT{memory_.data(), S(shape(v)), U(capacity(StackStorage))} {
    this->growUndef(ptrdiff_t(this->sz));
    MutArray<T, decltype(v.dim())>(this->data(), v.dim()) << v;
    // (*this) << v;
  }
  constexpr ManagedArray(const RowVector auto &v)
  requires(MatrixDimension<S>)
    : BaseT{memory_.data(), S(CartesianIndex(1, v.size())),
            U(capacity(StackStorage))} {
    this->growUndef(ptrdiff_t(this->sz));
    (*this) << v;
  }
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#else
#pragma clang diagnostic pop
#endif

  template <class D>
  constexpr auto operator=(
    const ManagedArray<T, D, StackStorage, A> &b) noexcept -> ManagedArray &
  requires(!std::same_as<S, D>)
  {
    // this condition implies `this->data() == nullptr`
    if (this->data() == b.data()) return *this;
    S d = b.dim();
    auto len = ptrdiff_t(d);
    this->growUndef(len);
    std::copy_n(b.data(), len, this->data());
    this->sz = d;
    return *this;
  }
  template <class D>
  constexpr auto
  operator=(ManagedArray<T, D, StackStorage, A> &&b) noexcept -> ManagedArray &
  requires(!std::same_as<S, D>)
  {
    // this condition implies `this->data() == nullptr`
    if (this->data() == b.data()) return *this;
    // here, we commandeer `b`'s memory
    S d = b.dim();
    this->allocator_ = std::move(b.get_allocator());
    // if `b` is small, we need to copy memory
    // no need to shrink our capacity
    if (b.isSmall()) std::copy_n(b.data(), ptrdiff_t(d), this->data());
    else this->maybeDeallocate(b.data(), ptrdiff_t(b.getCapacity()));
    b.resetNoFree();
    this->sz = d;
    return *this;
  }
  constexpr auto operator=(const ManagedArray &b) noexcept -> ManagedArray & {
    if (this == &b) return *this;
    S d = b.dim();
    auto len = ptrdiff_t(d);
    this->growUndef(len);
    std::copy_n(b.data(), len, this->data());
    this->sz = d;
    return *this;
  }
  constexpr auto operator=(ManagedArray &&b) noexcept -> ManagedArray & {
    if (this == &b) return *this;
    // here, we commandeer `b`'s memory
    S d = b.dim();
    this->allocator_ = std::move(b.get_allocator());
    if (b.isSmall()) std::copy_n(b.data(), ptrdiff_t(d), this->data());
    else this->maybeDeallocate(b.data(), ptrdiff_t(b.getCapacity()));
    b.resetNoFree();
    this->sz = d;
    return *this;
  }
  constexpr void resetNoFree() {
    this->ptr = memory_.data();
    this->sz = S{};
    this->capacity_ = capacity(StackStorage);
  }
  constexpr ~ManagedArray() noexcept { this->maybeDeallocate(); }

  [[nodiscard]] static constexpr auto identity(ptrdiff_t M) -> ManagedArray {
    static_assert(MatrixDimension<S>);
    ManagedArray B(SquareDims<>{row(M)}, T{0});
    B.diag() << 1;
    return B;
  }
  [[nodiscard]] static constexpr auto identity(Row<> R) -> ManagedArray {
    static_assert(MatrixDimension<S>);
    return identity(ptrdiff_t(R));
  }
  [[nodiscard]] static constexpr auto identity(Col<> C) -> ManagedArray {
    static_assert(MatrixDimension<S>);
    return identity(ptrdiff_t(C));
  }

  constexpr void reserveForGrow1()
  requires(std::same_as<S, Length<>>)
  {
    auto s = ptrdiff_t(this->sz), c = ptrdiff_t(this->capacity_);
    if (s == c) [[unlikely]]
      reserve(length(newCapacity(c)));
  }

  template <class... Args>
  constexpr auto emplace_back(Args &&...args) -> decltype(auto)
  requires(std::same_as<S, Length<>>)
  {
    reserveForGrow1();
    return this->emplace_back_within_capacity(args...);
  }
  constexpr void push_back(T value)
  requires(std::same_as<S, Length<>>)
  {
    reserveForGrow1();
    this->push_back_within_capacity(std::move(value));
  }
  constexpr auto insert(T *p, T x) -> T *requires(std::same_as<S, Length<>>) {
    auto s = ptrdiff_t(this->sz), c = ptrdiff_t(this->capacity_);
    if (s == c) [[unlikely]] {
      ptrdiff_t d = p - this->data();
      reserve(length(newCapacity(c)));
      p = this->data() + d;
    }
    invariant(s == ptrdiff_t(this->sz));
    return this->insert_within_capacity(p, std::move(x));
  }
  // behavior
  // if S is StridedDims, then we copy data.
  // If the new dims are larger in rows or cols, we fill with 0.
  // If the new dims are smaller in rows or cols, we truncate.
  // New memory outside of dims (i.e., stride larger), we leave uninitialized.
  //
  constexpr void resize(S nz) {
    reallocForSize(nz);
    this->sz = nz;
  }
  constexpr void resize(Row<> r)
  requires(std::same_as<S, Length<>> || MatrixDimension<S>)
  {
    if constexpr (std::same_as<S, Length<>>) return resize(S(r));
    else return resize(auto{this->sz}.set(r));
  }
  constexpr void resize(Col<> c)
  requires(std::same_as<S, Length<>> || MatrixDimension<S>)
  {
    if constexpr (std::same_as<S, Length<>>) return resize(S(c));
    else if constexpr (MatrixDimension<S>) return resize(auto{this->sz}.set(c));
  }
  constexpr void resizeForOverwrite(S M) {
    auto nz = ptrdiff_t(M);
    if (nz > ptrdiff_t(this->sz)) growUndef(nz);
    if constexpr (!std::is_trivially_destructible_v<T>) {
      ptrdiff_t oz = ptrdiff_t(this->sz);
      if (nz < oz) std::destroy_n(this->data() + nz, oz - nz);
      else if (nz > oz) // TODO: user should be in charge of initialization?
        for (ptrdiff_t i = oz; i < nz; ++i) std::construct_at(this->data() + i);
    }
    this->sz = M;
  }
  constexpr void resizeForOverwrite(ptrdiff_t M)
  requires(std::same_as<S, Length<>>)
  {
    resizeForOverwrite(length(M));
  }
  constexpr void resize(ptrdiff_t M)
  requires(std::same_as<S, Length<>>)
  {
    resize(length(M));
  }
  constexpr void reserve(ptrdiff_t M)
  requires(std::same_as<S, Length<>>)
  {
    reserve(length(M));
  }
  constexpr void resizeForOverwrite(Row<> r) {
    if constexpr (std::same_as<S, Length<>>) {
      return resizeForOverwrite(S(r));
    } else if constexpr (MatrixDimension<S>) {
      S nz = this->sz;
      return resizeForOverwrite(nz.set(r));
    }
  }
  constexpr void resizeForOverwrite(Col<> c) {
    if constexpr (std::same_as<S, Length<>>) {
      return resizeForOverwrite(S(c));
    } else if constexpr (MatrixDimension<S>) {
      S nz = this->sz;
      return resizeForOverwrite(nz.set(c));
    }
  }
  static constexpr auto
  reserveCore(S nz, storage_type *op, ptrdiff_t old_len, U oc,
              bool was_allocated) -> containers::Pair<storage_type *, U> {
    auto new_capacity = ptrdiff_t(nz);
    invariant(new_capacity >= 0z);
    if (new_capacity <= oc) return {op, oc};
    // allocate new, copy, deallocate old
    auto [new_ptr, new_cap] = alloc::alloc_at_least(A{}, new_capacity);
    auto nc = ptrdiff_t(new_cap);
    invariant(nc >= new_capacity);
    if (old_len) {
      if constexpr (std::is_trivially_move_constructible_v<T> &&
                    std::is_trivially_destructible_v<T>)
        std::memcpy(new_ptr, op, old_len * sizeof(storage_type));
      else std::uninitialized_move_n(op, old_len, new_ptr);
    }
    maybeDeallocate(op, old_len, ptrdiff_t(oc), was_allocated);
    return {new_ptr, math::capacity(nc)};
  }
  constexpr void reserve(S nz) {
    auto [np, nc] = reserveCore(nz, this->data(), ptrdiff_t(this->sz),
                                this->capacity_, wasAllocated());
    this->ptr = np;
    this->capacity_ = nc;
  }
  [[nodiscard]] constexpr auto get_allocator() const noexcept -> A {
    return allocator_;
  }
  // set size and 0.
  constexpr void setSize(Row<> r, Col<> c) {
    resizeForOverwrite({r, c});
    this->fill(0);
  }
  constexpr void resize(Row<> MM, Col<> NN) { resize(DenseDims{MM, NN}); }
  constexpr void reserve(Row<> M, Col<> N) {
    if constexpr (std::is_same_v<S, StridedDims<>>)
      reserve(StridedDims{M, N, max(N, RowStride(this->dim()))});
    else if constexpr (std::is_same_v<S, SquareDims<>>)
      reserve(SquareDims{row(std::max(ptrdiff_t(M), ptrdiff_t(N)))});
    else reserve(DenseDims{M, N});
  }
  constexpr void reserve(Row<> M, RowStride<> X) {
    if constexpr (std::is_same_v<S, StridedDims<>>)
      reserve(S{M, col(ptrdiff_t(X)), X});
    else if constexpr (std::is_same_v<S, SquareDims<>>)
      reserve(SquareDims{row(std::max(ptrdiff_t(M), ptrdiff_t(X)))});
    else reserve(S{M, col(ptrdiff_t(X))});
  }
  constexpr void clearReserve(Row<> M, Col<> N) {
    this->clear();
    reserve(M, N);
  }
  constexpr void clearReserve(Row<> M, RowStride<> X) {
    this->clear();
    reserve(M, X);
  }
  constexpr void resizeForOverwrite(Row<> M, Col<> N, RowStride<> X) {
    invariant(X >= N);
    if constexpr (std::is_same_v<S, StridedDims<>>)
      resizeForOverwrite(S{M, N, X});
    else if constexpr (std::is_same_v<S, SquareDims<>>) {
      invariant(ptrdiff_t(M) == ptrdiff_t(N));
      resizeForOverwrite(S{M});
    } else resizeForOverwrite(S{M, N});
  }
  constexpr void resizeForOverwrite(Row<> M, Col<> N) {
    if constexpr (std::is_same_v<S, StridedDims<>>)
      resizeForOverwrite(S{M, N, {ptrdiff_t(N)}});
    else if constexpr (std::is_same_v<S, SquareDims<>>) {
      invariant(ptrdiff_t(M) == ptrdiff_t(N));
      resizeForOverwrite(S{M});
    } else resizeForOverwrite(S{M, N});
  }

private:
  constexpr void allocateAtLeast(U len) {
    auto l = size_t(ptrdiff_t(len));
    alloc::AllocResult<storage_type> res = alloc::alloc_at_least(allocator_, l);
    this->ptr = res.ptr;
    invariant(res.count >= l);
    this->capacity_ = capacity(res.count);
  }
  [[nodiscard]] constexpr auto isSmall() const -> bool {
    invariant(this->capacity_ >= StackStorage);
    return this->capacity_ == StackStorage;
  }
  [[nodiscard]] constexpr auto wasAllocated() const -> bool {
    return !isSmall();
  }
  // this method should only be called from the destructor
  // (and the implementation taking the new ptr and capacity)
  void maybeDeallocate() noexcept {
    maybeDeallocate(this->data(), ptrdiff_t(this->sz),
                    ptrdiff_t(this->capacity_), wasAllocated());
  }
  static void maybeDeallocate(storage_type *p, ptrdiff_t sz, ptrdiff_t cap,
                              bool was_allocated) noexcept {
    if constexpr (!std::is_trivially_destructible_v<storage_type>)
      std::destroy_n(p, sz);
    if (was_allocated) A{}.deallocate(p, cap);
  }
  // this method should be called whenever the buffer lives
  // NOTE: it is invalid to reassign `sz` before calling `maybeDeallocate`!
  void maybeDeallocate(storage_type *newPtr, ptrdiff_t newCapacity) noexcept {
    maybeDeallocate(this->data(), ptrdiff_t(this->sz),
                    ptrdiff_t(this->capacity_), wasAllocated());
    this->ptr = newPtr;
    this->capacity_ = capacity(newCapacity);
  }
  // grow, discarding old data
  void growUndef(ptrdiff_t M) {
    invariant(M >= 0);
    if (M <= this->capacity_) return;
    maybeDeallocate();
    // because this doesn't care about the old data,
    // we can allocate after freeing, which may be faster
    this->ptr = this->allocator_.allocate(M);
    this->capacity_ = capacity(M);
#ifndef NDEBUG
    if constexpr (std::numeric_limits<T>::has_signaling_NaN)
      std::fill_n(this->data(), M, std::numeric_limits<T>::signaling_NaN());
    else if constexpr (std::integral<T>)
      std::fill_n(this->data(), M, std::numeric_limits<T>::min());
#endif
  }
  constexpr void reallocForSize(S nz) {
    S oz = this->sz;
    if constexpr (std::same_as<S, Length<>>) {
      auto ozs = ptrdiff_t(oz), nzs = ptrdiff_t(nz);
      if (nz <= oz) {
        if constexpr (!std::is_trivially_destructible_v<T>)
          if (nz < oz) std::destroy_n(this->data() + nzs, ozs - nzs);
        return;
      }
      if (nz > this->capacity_) {
        auto ncu = size_t(nzs);
        auto [newPtr, newCap] = alloc::alloc_at_least(allocator_, ncu);
        invariant(newCap >= ncu);
        auto ncs = ptrdiff_t(newCap);
        invariant(ncs >= nzs);
        if (oz) std::copy_n(this->data(), ozs, newPtr);
        maybeDeallocate(newPtr, ncs);
      }
      if constexpr (!std::is_trivially_destructible_v<T>)
        for (ptrdiff_t i = ozs; i < nz; ++i)
          std::construct_at(this->data() + i);
      else std::fill(this->data() + ozs, this->data() + nzs, T{});
    } else {
      static_assert(std::is_trivially_destructible_v<T>,
                    "Resizing matrices holding non-is_trivially_destructible_v "
                    "objects is not yet supported.");
      static_assert(MatrixDimension<S>, "Can only resize 1 or 2d containers.");
      auto len = ptrdiff_t(nz);
      if (len == 0) return;
      auto new_x = ptrdiff_t{RowStride(nz)}, old_x = ptrdiff_t{RowStride(oz)},
           new_n = ptrdiff_t{Col(nz)}, old_n = ptrdiff_t{Col(oz)},
           new_m = ptrdiff_t{Row(nz)}, old_m = ptrdiff_t{Row(oz)};
      bool new_alloc = len > this->capacity_;
      bool in_place = !new_alloc;
      T *npt = this->data();
      if (new_alloc) {
        alloc::AllocResult<T> res = alloc::alloc_at_least(allocator_, len);
        npt = res.ptr;
        len = res.count;
      }
      // we can copy forward so long as the new stride is smaller
      // so that the start of the dst range is outside of the src range
      // we can also safely forward copy if we allocated a new ptr
      bool forward_copy = (new_x <= old_x) || new_alloc;
      ptrdiff_t cols_to_copy = std::min(old_n, new_n);
      // we only need to copy if memory shifts position
      bool copy_cols = new_alloc || ((cols_to_copy > 0) && (new_x != old_x));
      // if we're in place, we have 1 less row to copy
      ptrdiff_t rows_to_copy = std::min(old_m, new_m);
      ptrdiff_t fill_count = new_n - cols_to_copy;
      if ((rows_to_copy) && (copy_cols || fill_count)) {
        if (forward_copy) {
          // truncation, we need to copy rows to increase stride
          T *src = this->data();
          T *dst = npt;
          do {
            if (copy_cols && (!in_place)) std::copy_n(src, cols_to_copy, dst);
            if (fill_count) std::fill_n(dst + cols_to_copy, fill_count, T{});
            src += old_x;
            dst += new_x;
            in_place = false;
          } while (--rows_to_copy);
        } else /* [[unlikely]] */ {
          // backwards copy, only needed when we increasing stride but not
          // reallocating, which should be comparatively uncommon.
          // Should probably benchmark or determine actual frequency
          // before adding `[[unlikely]]`.
          invariant(in_place);
          T *src = this->data() + (rows_to_copy + in_place) * old_x;
          T *dst = npt + (rows_to_copy + in_place) * new_x;
          do {
            src -= old_x;
            dst -= new_x;
            if (cols_to_copy && (rows_to_copy > in_place))
              std::copy_backward(src, src + cols_to_copy, dst + cols_to_copy);
            if (fill_count) std::fill_n(dst + cols_to_copy, fill_count, T{});
          } while (--rows_to_copy);
        }
      }
      // zero init remaining rows
      for (ptrdiff_t m = old_m; m < new_m; ++m)
        std::fill_n(npt + m * new_x, new_n, T{});
      if (new_alloc) maybeDeallocate(npt, len);
    }
  }

  friend void PrintTo(const ManagedArray &x, ::std::ostream *os)
  requires(utils::Printable<T>)
  {
    *os << x;
  }

  [[no_unique_address]] A allocator_{};
  [[no_unique_address]] containers::Storage<storage_type, StackStorage> memory_;
};

static_assert(std::move_constructible<ManagedArray<intptr_t, Length<>>>);
static_assert(std::copyable<ManagedArray<intptr_t, Length<>>>);
// Check that `[[no_unique_address]]` is working.
// sizes should be:
// [ptr, dims, capacity, allocator_, array]
// 8 + 3*4 + 4 + 0 + 64*8 = 24 + 512 = 536
static_assert(sizeof(ManagedArray<int64_t, StridedDims<>, 64,
                                  alloc::Mallocator<int64_t>>) == 552);
// sizes should be:
// [ptr, dims, capacity, allocator_, array]
// 8 + 2*4 + 8 + 0 + 64*8 = 24 + 512 = 536
static_assert(
  sizeof(ManagedArray<int64_t, DenseDims<>, 64, alloc::Mallocator<int64_t>>) ==
  544);
// sizes should be:
// [ptr, dims, capacity, allocator_, array]
// 8 + 1*4 + 4 + 0 + 64*8 = 16 + 512 = 528
static_assert(
  sizeof(ManagedArray<int64_t, SquareDims<>, 64, alloc::Mallocator<int64_t>>) ==
  536);

template <class T, ptrdiff_t N = containers::PreAllocStorage<T, ptrdiff_t>()>
using Vector = ManagedArray<T, Length<>, N>;
template <class T> using PtrVector = Array<T, Length<>>;
template <class T> using MutPtrVector = MutArray<T, Length<>>;

static_assert(std::move_constructible<Vector<intptr_t>>);
static_assert(std::copy_constructible<Vector<intptr_t>>);
static_assert(std::copyable<Vector<intptr_t>>);
static_assert(AbstractVector<Array<int64_t, Length<>>>);
static_assert(AbstractVector<MutArray<int64_t, Length<>>>);
static_assert(AbstractVector<Vector<int64_t>>);
static_assert(!AbstractVector<int64_t>);
static_assert(!std::is_trivially_copyable_v<Vector<int64_t>>);
static_assert(!std::is_trivially_destructible_v<Vector<int64_t>>);

template <typename T> using StridedVector = Array<T, StridedRange>;
template <typename T> using MutStridedVector = MutArray<T, StridedRange>;

static_assert(AbstractVector<StridedVector<int64_t>>);
static_assert(AbstractVector<MutStridedVector<int64_t>>);
static_assert(std::is_trivially_copyable_v<StridedVector<int64_t>>);

template <class T, ptrdiff_t R = -1, ptrdiff_t C = -1, ptrdiff_t X = -1>
using PtrMatrix = Array<T, StridedDims<R, C, X>>;
template <class T, ptrdiff_t R = -1, ptrdiff_t C = -1, ptrdiff_t X = -1>
using MutPtrMatrix = MutArray<T, StridedDims<R, C, X>>;
template <class T,
          ptrdiff_t L = containers::PreAllocStorage<T, StridedDims<>>()>
using Matrix = ManagedArray<T, StridedDims<>, L>;
template <class T, ptrdiff_t R = -1, ptrdiff_t C = -1>
using DensePtrMatrix = Array<T, DenseDims<R, C>>;
template <class T, ptrdiff_t R = -1, ptrdiff_t C = -1>
using MutDensePtrMatrix = MutArray<T, DenseDims<R, C>>;
template <class T, ptrdiff_t L = containers::PreAllocStorage<T, DenseDims<>>()>
using DenseMatrix = ManagedArray<T, DenseDims<>, L>;
template <class T> using SquarePtrMatrix = Array<T, SquareDims<>>;
template <class T> using MutSquarePtrMatrix = MutArray<T, SquareDims<>>;
template <class T, ptrdiff_t L = containers::PreAllocStorage<T, SquareDims<>>()>
using SquareMatrix = ManagedArray<T, SquareDims<>, L>;

static_assert(sizeof(PtrMatrix<int64_t>) ==
              3 * sizeof(ptrdiff_t) + sizeof(int64_t *));
static_assert(sizeof(MutPtrMatrix<int64_t>) ==
              3 * sizeof(ptrdiff_t) + sizeof(int64_t *));
static_assert(sizeof(DensePtrMatrix<int64_t>) ==
              2 * sizeof(ptrdiff_t) + sizeof(int64_t *));
static_assert(sizeof(MutDensePtrMatrix<int64_t>) ==
              2 * sizeof(ptrdiff_t) + sizeof(int64_t *));
static_assert(sizeof(SquarePtrMatrix<int64_t>) ==
              sizeof(ptrdiff_t) + sizeof(int64_t *));
static_assert(sizeof(MutSquarePtrMatrix<int64_t>) ==
              sizeof(ptrdiff_t) + sizeof(int64_t *));
static_assert(std::is_trivially_copyable_v<PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> is not trivially copyable!");
static_assert(std::is_trivially_copyable_v<PtrVector<int64_t>>,
              "PtrVector<int64_t,0> is not trivially copyable!");
// static_assert(std::is_trivially_copyable_v<MutPtrMatrix<int64_t>>,
//               "MutPtrMatrix<int64_t> is not trivially copyable!");
static_assert(sizeof(ManagedArray<int32_t, DenseDims<3, 5>, 15>) ==
              sizeof(int32_t *) + 16 * sizeof(int32_t));
static_assert(sizeof(ManagedArray<int32_t, DenseDims<>, 15>) ==
              sizeof(int32_t *) + 3 * sizeof(ptrdiff_t) + 16 * sizeof(int32_t));

static_assert(!AbstractVector<PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractVector succeeded");
static_assert(!AbstractVector<MutPtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractVector succeeded");
static_assert(!AbstractVector<const PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractVector succeeded");

static_assert(AbstractMatrix<PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");
static_assert(
  std::same_as<decltype(PtrMatrix<int64_t>(nullptr, row(0),
                                           col(0))[ptrdiff_t(0), ptrdiff_t(0)]),
               const int64_t &>);
static_assert(
  std::same_as<std::remove_reference_t<decltype(MutPtrMatrix<int64_t>(
                 nullptr, row(0), col(0))[ptrdiff_t(0), ptrdiff_t(0)])>,
               int64_t>);

static_assert(AbstractMatrix<MutPtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");
static_assert(AbstractMatrix<const PtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");
static_assert(AbstractMatrix<const MutPtrMatrix<int64_t>>,
              "PtrMatrix<int64_t> isa AbstractMatrix failed");

static_assert(AbstractVector<MutPtrVector<int64_t>>,
              "PtrVector<int64_t> isa AbstractVector failed");
static_assert(AbstractVector<PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractVector failed");
static_assert(AbstractVector<const PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractVector failed");
static_assert(AbstractVector<const MutPtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractVector failed");

static_assert(AbstractVector<Vector<int64_t>>,
              "PtrVector<int64_t> isa AbstractVector failed");

static_assert(!AbstractMatrix<MutPtrVector<int64_t>>,
              "PtrVector<int64_t> isa AbstractMatrix succeeded");
static_assert(!AbstractMatrix<PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractMatrix succeeded");
static_assert(!AbstractMatrix<const PtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractMatrix succeeded");
static_assert(!AbstractMatrix<const MutPtrVector<int64_t>>,
              "PtrVector<const int64_t> isa AbstractMatrix succeeded");
static_assert(std::is_convertible_v<DenseMatrix<int64_t>, Matrix<int64_t>>);
static_assert(
  std::is_convertible_v<DenseMatrix<int64_t>, DensePtrMatrix<int64_t>>);
static_assert(std::is_convertible_v<DenseMatrix<int64_t>, PtrMatrix<int64_t>>);
static_assert(std::is_convertible_v<SquareMatrix<int64_t>, Matrix<int64_t>>);
static_assert(
  std::is_convertible_v<SquareMatrix<int64_t>, MutPtrMatrix<int64_t>>);

template <class T, class S>
ManagedArray(std::type_identity<T>, S s) -> ManagedArray<T, S>;

template <class S> using IntArray = Array<int64_t, S>;
template <VectorDimension S = ptrdiff_t>
using IntVector = ManagedArray<int64_t, S>;
template <MatrixDimension S = DenseDims<>>
using IntMatrix = ManagedArray<int64_t, S>;

static_assert(std::same_as<IntMatrix<>::value_type, int64_t>);
static_assert(AbstractMatrix<IntMatrix<>>);
static_assert(std::copyable<IntMatrix<>>);
static_assert(std::same_as<utils::eltype_t<Matrix<int64_t>>, int64_t>);

static_assert(std::convertible_to<Array<int64_t, SquareDims<>>,
                                  Array<int64_t, StridedDims<>>>);

static_assert(std::same_as<const int64_t &,
                           decltype(std::declval<PtrMatrix<int64_t>>()[0, 0])>);
static_assert(std::is_trivially_copyable_v<MutArray<int64_t, Length<>>>);
/// \brief Returns the maximum number of digits per column of a matrix.
constexpr auto getMaxDigits(PtrMatrix<Rational> A) -> Vector<ptrdiff_t> {
  ptrdiff_t M = ptrdiff_t(A.numRow());
  ptrdiff_t N = ptrdiff_t(A.numCol());
  Vector<ptrdiff_t> maxDigits{length(N), 0};
  invariant(ptrdiff_t(maxDigits.size()), N);
  // this is slow, because we count the digits of every element
  // we could optimize this by reducing the number of calls to countDigits
  for (ptrdiff_t i = 0; i < M; i++) {
    for (ptrdiff_t j = 0; j < N; j++) {
      ptrdiff_t c = countDigits(A[i, j]);
      maxDigits[j] = std::max(maxDigits[j], c);
    }
  }
  return maxDigits;
}

/// Returns the number of digits of the largest number in the matrix.
template <std::integral T>
constexpr auto getMaxDigits(PtrMatrix<T> A) -> Vector<T> {
  ptrdiff_t M = ptrdiff_t(A.numRow());
  ptrdiff_t N = ptrdiff_t(A.numCol());
  Vector<T> maxDigits{length(N), T{}};
  invariant(ptrdiff_t(maxDigits.size()), N);
  // first, we find the digits with the maximum value per column
  for (ptrdiff_t i = 0; i < M; i++) {
    for (ptrdiff_t j = 0; j < N; j++) {
      // negative numbers need one more digit
      // first, we find the maximum value per column,
      // dividing positive numbers by -10
      T Aij = A[i, j];
      if constexpr (std::signed_integral<T>)
        maxDigits[j] = std::min(maxDigits[j], Aij > 0 ? Aij / -10 : Aij);
      else maxDigits[j] = std::max(maxDigits[j], Aij);
    }
  }
  // then, we count the digits of the maximum value per column
  for (ptrdiff_t j = 0; j < maxDigits.size(); j++)
    maxDigits[j] = utils::countDigits(maxDigits[j]);
  return maxDigits;
}

template <typename T>
inline auto printMatrix(std::ostream &os, PtrMatrix<T> A) -> std::ostream & {
  // std::ostream &printMatrix(std::ostream &os, T const &A) {
  auto [M, N] = shape(A);
  if ((!M) || (!N)) return os << "[ ]";
  // first, we determine the number of digits needed per column
  auto maxDigits{getMaxDigits(A)};
  using U = decltype(countDigits(std::declval<T>()));
  for (ptrdiff_t i = 0; i < M; i++) {
    if (i) os << "  ";
    else os << "\n[ ";
    for (ptrdiff_t j = 0; j < N; j++) {
      auto Aij = A[i, j];
      for (U k = 0; k < U(maxDigits[j]) - countDigits(Aij); k++) os << " ";
      os << Aij;
      if (j != ptrdiff_t(N) - 1) os << " ";
      else if (i != ptrdiff_t(M) - 1) os << "\n";
    }
  }
  return os << " ]";
}
// We mirror `A` with a matrix of integers indicating sizes, and a vectors of
// chars. We fill the matrix with the number of digits of each element, and
// the vector with the characters of each element. We could use a vector of
// vectors of chars to avoid needing to copy memory on reallocation, but this
// would yield more complicated management. We should also generally be able
// to avoid allocations. We can use a Vector with a lot of initial capacity,
// and then resize based on a conservative estimate of the number of chars per
// elements.
inline auto printMatrix(std::ostream &os,
                        PtrMatrix<double> A) -> std::ostream & {
  // std::ostream &printMatrix(std::ostream &os, T const &A) {
  auto [M, N] = shape(A);
  if ((!M) || (!N)) return os << "[ ]";
  // first, we determine the number of digits needed per column
  Vector<char, 512> digits;
  digits.resizeForOverwrite(512);
  // we can't have more than 255 digits
  DenseMatrix<uint8_t> numDigits{DenseDims<>{row(M), col(N)}};
  char *ptr = digits.begin();
  char *pEnd = digits.end();
  for (ptrdiff_t m = 0; m < M; m++) {
    for (ptrdiff_t n = 0; n < N; n++) {
      auto Aij = A[m, n];
      while (true) {
        auto [p, ec] = std::to_chars(ptr, pEnd, Aij);
        if (ec == std::errc()) [[likely]] {
          numDigits[m, n] = std::distance(ptr, p);
          ptr = p;
          break;
        }
        // we need more space
        ptrdiff_t elemSoFar = m * ptrdiff_t(N) + n;
        ptrdiff_t charSoFar = std::distance(digits.begin(), ptr);
        // cld
        ptrdiff_t charPerElem = (charSoFar + elemSoFar - 1) / elemSoFar;
        ptrdiff_t newCapacity =
          (1 + charPerElem) * M * N; // +1 for good measure
        digits.resize(newCapacity);
        ptr = digits.begin() + charSoFar;
        pEnd = digits.end();
      }
    }
  }
  Vector<uint8_t> maxDigits;
  maxDigits.resizeForOverwrite(N);
  maxDigits << numDigits[0, _];
  for (ptrdiff_t m = 0; m < M; m++)
    for (ptrdiff_t n = 0; n < N; n++)
      maxDigits[n] = std::max(maxDigits[n], numDigits[m, n]);

  ptr = digits.begin();
  // we will allocate 512 bytes at a time
  for (ptrdiff_t i = 0; i < M; i++) {
    if (i) os << "  ";
    else os << "\n[ ";
    for (ptrdiff_t j = 0; j < N; j++) {
      ptrdiff_t nD = numDigits[i, j];
      for (ptrdiff_t k = 0; k < maxDigits[j] - nD; k++) os << " ";
      os << std::string_view(ptr, nD);
      if (j != ptrdiff_t(N) - 1) os << " ";
      else if (i != ptrdiff_t(M) - 1) os << "\n";
      ptr += nD;
    }
  }
  return os << " ]";
}

template <class T, class S, class P>
template <typename Op, typename RHS>
void ArrayOps<T, S, P>::vcopyTo(const RHS &B, Op op) {
  // static_assert(sizeof(utils::eltype_t<decltype(B)>) <= 8);
  MutArray<T, S, !std::same_as<T *, decltype(data_())>> self{Self()};
  // P &self{Self()};
  auto [M, N] = promote_shape(self, B);
  constexpr bool assign = std::same_as<Op, utils::CopyAssign>;
  using PT = utils::promote_eltype_t<P, RHS>;
#ifdef CASTTOSCALARIZE
  using E = math::scalarize_via_cast_t<
    std::remove_cvref_t<decltype(std::declval<P>().view())>>;
  if constexpr (!std::same_as<E, void> &&
                ((math::ScalarizeViaCastTo<E, decltype(B)>()) ||
                 (std::same_as<std::remove_cvref_t<decltype(B)>, double> &&
                  std::same_as<Op, std::multiplies<>>))) {
    auto d{reinterpret<E>(Self())};
    if constexpr (assign) d << reinterpret<E>(B);
    else d << op(d, reinterpret<E>(B));
#ifndef POLYMATHNOEXPLICITSIMDARRAY
  } else if constexpr (simd::SIMDSupported<PT>) {
#else
  } else if constexpr (AbstractVector<P>) {
#endif
#elifndef POLYMATHNOEXPLICITSIMDARRAY
  if constexpr (simd::SIMDSupported<PT>) {
#else
  if constexpr (AbstractVector<P>) {
#endif
#ifndef POLYMATHNOEXPLICITSIMDARRAY
    if constexpr (IsOne<decltype(M)>)
      vcopyToSIMD(self, B, N, utils::NoRowIndex{}, op);
    else if constexpr (IsOne<decltype(N)>)
      vcopyToSIMD(self, B, M, utils::NoRowIndex{}, op);
    else if constexpr (StaticInt<decltype(M)>) {
      constexpr std::array<ptrdiff_t, 2> UIR = unrollf<ptrdiff_t(M)>();
      constexpr ptrdiff_t U = UIR[0];
      if constexpr (U != 0)
        for (ptrdiff_t r = 0; r < (M - U + 1); r += U)
          vcopyToSIMD(self, B, N, simd::index::Unroll<U>{r}, op);
      constexpr ptrdiff_t R = UIR[1];
      if constexpr (R != 0)
        vcopyToSIMD(self, B, N, simd::index::Unroll<R>{M - R}, op);
    } else {
      ptrdiff_t r = 0;
      for (; r < (M - 3); r += 4)
        vcopyToSIMD(self, B, N, simd::index::Unroll<4>{r}, op);
      switch (M & 3) {
      case 0: return;
      case 1: return vcopyToSIMD(self, B, N, simd::index::Unroll<1>{r}, op);
      case 2: return vcopyToSIMD(self, B, N, simd::index::Unroll<2>{r}, op);
      default: return vcopyToSIMD(self, B, N, simd::index::Unroll<3>{r}, op);
      }
    }
  } else if constexpr (AbstractVector<P>) {
#endif
    ptrdiff_t L = IsOne<decltype(N)> ? M : N;
    constexpr bool isstatic =
      IsOne<decltype(N)> ? StaticInt<decltype(M)> : StaticInt<decltype(N)>;
    if constexpr (!std::is_copy_assignable_v<PT> && assign) {
      POLYMATHIVDEP
      for (ptrdiff_t j = 0; j < L; ++j)
        if constexpr (std::convertible_to<decltype(B), PT>) self[j] = auto{B};
        else self[j] = auto{B[j]};
    } else if constexpr (isstatic) {
      POLYMATHFULLUNROLL
      for (ptrdiff_t j = 0; j < L; ++j)
        utils::assign(self, B, utils::NoRowIndex{}, j, op);
    } else {
      POLYMATHIVDEP
      for (ptrdiff_t j = 0; j < L; ++j)
        utils::assign(self, B, utils::NoRowIndex{}, j, op);
    }
  } else {
    ptrdiff_t R = ptrdiff_t(M), C = ptrdiff_t(N);
    POLYMATHNOVECTORIZE
    for (ptrdiff_t i = 0; i < R; ++i) {
      if constexpr (!std::is_copy_assignable_v<PT> && assign) {
        POLYMATHIVDEP
        for (ptrdiff_t j = 0; j < C; ++j)
          if constexpr (std::convertible_to<decltype(B), PT>)
            self[i, j] = auto{B};
          else if constexpr (RowVector<decltype(B)>) self[i, j] = auto{B[j]};
          else if constexpr (ColVector<decltype(B)>) self[i, j] = auto{B[i]};
          else self[i, j] = auto{B[i, j]};
      } else {
        POLYMATHIVDEP
        for (ptrdiff_t j = 0; j < C; ++j) utils::assign(self, B, i, j, op);
      }
    }
  }
}

} // namespace poly::math

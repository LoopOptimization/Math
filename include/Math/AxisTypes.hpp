#pragma once
#include "Utilities/Invariant.hpp"
#include <cstddef>
#include <iostream>
#include <limits>
#include <type_traits>

/// LinAlg
///
/// This is the namespace for all mathematical functions.
/// Semantics:
/// We generally around structs holding pointer/size/stride/capacity information
/// by value.
/// For assignments that copy the underlying data, use `<<`
/// E.g., `A << B + C;`
/// Updating assignment operators like `+=`, `-=`, and, `*=` are supported.
///
/// The choice of `<<` over `=` for copying data is so that `operator=`
/// can be the usual copy assignment operator, which is useful to use when
/// we pass arrays/pointers by value to functions that truncate their size.
///
/// Operations like `+` and `-` are elementwise, while `*` performs matrix
/// multiplication. All operations are lazy, building up expression templates
/// that are evaluated upon assignments, e.g. `<<` or `+=`.
///
/// All the PtrVector/PtrMatrix types are trivially destructible, copyable, etc
/// Their lifetimes are governed by the Arena or RAII type used to back
/// them.
namespace poly::math {

using utils::invariant;

template <ptrdiff_t M = -1> struct Length {
  // static constexpr ptrdiff_t len = M;
  static_assert(M >= 0);
  [[gnu::artificial, gnu::always_inline]] explicit inline constexpr
  operator ptrdiff_t() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] explicit inline constexpr
  operator bool() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] static inline constexpr auto
  staticint() {
    return std::integral_constant<ptrdiff_t, M>{};
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr
  operator Length<-1>() const;
};
template <> struct Length<-1> {
  enum class len : ptrdiff_t {};
  len M;
  [[gnu::artificial, gnu::always_inline]] explicit inline constexpr
  operator ptrdiff_t() const {
    auto m = static_cast<ptrdiff_t>(M);
    invariant(m >= 0);
    return m;
  }
  [[gnu::artificial, gnu::always_inline]] explicit inline constexpr
  operator bool() const {
    return static_cast<ptrdiff_t>(M);
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++() -> Length & {
    M = static_cast<len>(static_cast<ptrdiff_t>(M) + 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--() -> Length & {
    M = static_cast<len>(static_cast<ptrdiff_t>(M) - 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++(int) -> Length {
    Length tmp{*this};
    M = static_cast<len>(static_cast<ptrdiff_t>(M) + 1z);
    return tmp;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--(int) -> Length {
    Length tmp{*this};
    M = static_cast<len>(static_cast<ptrdiff_t>(M) + 1z);
    return tmp;
  }
};

template <ptrdiff_t M>
[[gnu::artificial,
  gnu::always_inline]] inline constexpr Length<M>::operator Length<-1>() const {
  return {static_cast<Length<-1>::len>(M)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(ptrdiff_t x, Length<> y) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Length<> y, ptrdiff_t x) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Length<> x, Length<> y) -> bool {
  return ptrdiff_t(x) == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(ptrdiff_t x, Length<> y) -> std::strong_ordering {
  return x <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Length<> x, ptrdiff_t y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> y;
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Length<> x, Length<> y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> ptrdiff_t(y);
}

template <ptrdiff_t M = -1, std::integral I = ptrdiff_t> struct Capacity {
  static_assert(M >= 0);
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator ptrdiff_t() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr
  operator Capacity<-1, I>() const;
};
template <std::integral I> struct Capacity<-1, I> {
  enum class cap : I {};
  cap M;
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator I() const {
    auto m = static_cast<I>(M);
    invariant(m >= 0);
    return m;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return static_cast<I>(M);
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++() -> Capacity & {
    M = static_cast<cap>(static_cast<I>(M) + 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--() -> Capacity & {
    M = static_cast<cap>(static_cast<I>(M) - 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++(int) -> Capacity {
    Capacity tmp{*this};
    M = static_cast<cap>(static_cast<I>(M) + 1z);
    return tmp;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--(int) -> Capacity {
    Capacity tmp{*this};
    M = static_cast<cap>(static_cast<I>(M) + 1z);
    return tmp;
  }
};
template <ptrdiff_t M, std::integral I>
[[gnu::artificial, gnu::always_inline]] inline constexpr Capacity<
  M, I>::operator Capacity<-1, I>() const {
  return {static_cast<Capacity<-1, I>::cap>(M)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(ptrdiff_t x, Capacity<> y) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Capacity<> y, ptrdiff_t x) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Capacity<> x, Capacity<> y) -> bool {
  return ptrdiff_t(x) == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(ptrdiff_t x, Capacity<> y) -> std::strong_ordering {
  return x <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Capacity<> x, ptrdiff_t y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> y;
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Length<> x, Capacity<> y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Capacity<> x, Length<> y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Capacity<> x, Capacity<> y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> ptrdiff_t(y);
}

template <ptrdiff_t M = -1> struct Row {
  static_assert(M >= 0);
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator ptrdiff_t() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr
  operator Row<-1>() const;
};
template <> struct Row<-1> {
  enum class row : ptrdiff_t {};
  row M;
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator ptrdiff_t() const {
    auto m = static_cast<ptrdiff_t>(M);
    invariant(m >= 0);
    return m;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return static_cast<ptrdiff_t>(M);
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++() -> Row & {
    M = static_cast<row>(static_cast<ptrdiff_t>(M) + 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--() -> Row & {
    M = static_cast<row>(static_cast<ptrdiff_t>(M) - 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++(int) -> Row {
    Row tmp{*this};
    M = static_cast<row>(static_cast<ptrdiff_t>(M) + 1z);
    return tmp;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--(int) -> Row {
    Row tmp{*this};
    M = static_cast<row>(static_cast<ptrdiff_t>(M) + 1z);
    return tmp;
  }
};
static_assert(sizeof(Row<>) == sizeof(ptrdiff_t));
template <ptrdiff_t M>
[[gnu::artificial,
  gnu::always_inline]] inline constexpr Row<M>::operator Row<-1>() const {
  return {static_cast<Row<-1>::row>(M)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(ptrdiff_t x, Row<> y) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Row<> y, ptrdiff_t x) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Row<> x, Row<> y) -> bool {
  return ptrdiff_t(x) == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(ptrdiff_t x, Row<> y) -> std::strong_ordering {
  return x <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Row<> x, ptrdiff_t y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> y;
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Row<> x, Row<> y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> ptrdiff_t(y);
}
inline auto operator<<(std::ostream &os, Row<> x) -> std::ostream & {
  return os << "Row<>{" << ptrdiff_t(x) << "}";
}
template <ptrdiff_t M = -1> struct Col {
  static_assert(M >= 0);
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator ptrdiff_t() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr
  operator Col<-1>() const;
};
template <> struct Col<-1> {
  enum class col : ptrdiff_t {};
  col M;
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator ptrdiff_t() const {
    auto m = static_cast<ptrdiff_t>(M);
    invariant(m >= 0);
    return m;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return static_cast<ptrdiff_t>(M);
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++() -> Col & {
    M = static_cast<col>(static_cast<ptrdiff_t>(M) + 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--() -> Col & {
    M = static_cast<col>(static_cast<ptrdiff_t>(M) - 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++(int) -> Col {
    Col tmp{*this};
    M = static_cast<col>(static_cast<ptrdiff_t>(M) + 1z);
    return tmp;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--(int) -> Col {
    Col tmp{*this};
    M = static_cast<col>(static_cast<ptrdiff_t>(M) - 1z);
    return tmp;
  }
};
static_assert(sizeof(Col<>) == sizeof(ptrdiff_t));
template <ptrdiff_t M>
[[gnu::artificial,
  gnu::always_inline]] inline constexpr Col<M>::operator Col<-1>() const {
  return {static_cast<Col<-1>::col>(M)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(ptrdiff_t x, Col<> y) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Col<> y, ptrdiff_t x) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Col<> x, Col<> y) -> bool {
  return ptrdiff_t(x) == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(ptrdiff_t x, Col<> y) -> std::strong_ordering {
  return x <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Col<> x, ptrdiff_t y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> y;
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Col<> x, Col<> y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator*(Row<> r, Col<> c) -> ptrdiff_t {
  return ptrdiff_t(r) * ptrdiff_t(c);
}
inline auto operator<<(std::ostream &os, Col<> x) -> std::ostream & {
  return os << "Col<>{" << ptrdiff_t(x) << "}";
}
template <ptrdiff_t M = -1> struct RowStride {
  static_assert(M >= 0);
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator ptrdiff_t() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return M;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr
  operator RowStride<-1>() const;
};
template <> struct RowStride<-1> {
  enum class stride : ptrdiff_t {};
  stride M;
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator ptrdiff_t() const {
    auto m = static_cast<ptrdiff_t>(M);
    invariant(m >= 0);
    return m;
  }
  [[gnu::artificial, gnu::always_inline]] inline explicit constexpr
  operator bool() const {
    return static_cast<ptrdiff_t>(M);
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++() -> RowStride & {
    M = static_cast<stride>(static_cast<ptrdiff_t>(M) + 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--() -> RowStride & {
    M = static_cast<stride>(static_cast<ptrdiff_t>(M) - 1z);
    return *this;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator++(int) -> RowStride {
    RowStride tmp{*this};
    M = static_cast<stride>(static_cast<ptrdiff_t>(M) + 1z);
    return tmp;
  }
  [[gnu::artificial, gnu::always_inline]] inline constexpr auto
  operator--(int) -> RowStride {
    RowStride tmp{*this};
    M = static_cast<stride>(static_cast<ptrdiff_t>(M) - 1z);
    return tmp;
  }
};
static_assert(sizeof(RowStride<>) == sizeof(ptrdiff_t));
template <ptrdiff_t M>
[[gnu::artificial,
  gnu::always_inline]] inline constexpr RowStride<M>::operator RowStride<-1>()
  const {
  return {static_cast<RowStride<-1>::stride>(M)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(ptrdiff_t x, RowStride<> y) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(RowStride<> y, ptrdiff_t x) -> bool {
  return x == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(RowStride<> x, RowStride<> y) -> bool {
  return ptrdiff_t(x) == ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(ptrdiff_t x, RowStride<> y) -> std::strong_ordering {
  return x <=> ptrdiff_t(y);
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(RowStride<> x, ptrdiff_t y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> y;
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(RowStride<> x, RowStride<> y) -> std::strong_ordering {
  return ptrdiff_t(x) <=> ptrdiff_t(y);
}
inline auto operator<<(std::ostream &os, RowStride<> x) -> std::ostream & {
  return os << "RowStride<>{" << ptrdiff_t(x) << "}";
}

// constexpr auto max(Row M, Col N) -> ptrdiff_t {
//   return std::max(ptrdiff_t(M), ptrdiff_t(N));
// }
// constexpr auto max(Col N, RowStride X) -> RowStride {
//   return RowStride{std::max(ptrdiff_t(N), ptrdiff_t(X))};
// }
// constexpr auto min(Col N, Col X) -> Col {
//   return Col{std::max(Col::V(N), Col::V(X))};
// }
// constexpr auto min(Row N, Col X) -> ptrdiff_t {
//   return std::min(ptrdiff_t(N), ptrdiff_t(X));
// }

template <ptrdiff_t M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
standardizeRangeBound(Row<M> x) {
  if constexpr (M == -1) return ptrdiff_t(x);
  else return std::integral_constant<ptrdiff_t, M>{};
}
template <ptrdiff_t M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
standardizeRangeBound(Col<M> x) {
  if constexpr (M == -1) return ptrdiff_t(x);
  else return std::integral_constant<ptrdiff_t, M>{};
}

template <ptrdiff_t M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
unwrapRow(Row<M> x) {
  if constexpr (M == -1) return ptrdiff_t(x);
  else return std::integral_constant<ptrdiff_t, M>{};
}
template <ptrdiff_t M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
unwrapCol(Col<M> x) {
  if constexpr (M == -1) return ptrdiff_t(x);
  else return std::integral_constant<ptrdiff_t, M>{};
}
template <ptrdiff_t M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
unwrapStride(RowStride<M> x) {
  if constexpr (M == -1) return ptrdiff_t(x);
  else return std::integral_constant<ptrdiff_t, M>{};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
unwrapRow(auto x) {
  return x;
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
unwrapCol(auto x) {
  return x;
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
unwrapStride(auto x) {
  return x;
}

template <ptrdiff_t C, ptrdiff_t X>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator==(Col<C> c, RowStride<X> x) -> bool {
  return ptrdiff_t(c) == ptrdiff_t(x);
}
template <ptrdiff_t C, ptrdiff_t X>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator<=>(Col<C> c, RowStride<X> x) -> std::strong_ordering {
  return ptrdiff_t(c) <=> ptrdiff_t(x);
}

[[gnu::artificial, gnu::always_inline]] inline constexpr auto
row(ptrdiff_t x) -> Row<> {
  invariant(x >= 0);
  return Row<-1>{static_cast<Row<-1>::row>(x)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
col(ptrdiff_t x) -> Col<> {
  invariant(x >= 0);
  return Col<-1>{static_cast<Col<-1>::col>(x)};
}
template <ptrdiff_t M>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
col(Length<M> x) -> Col<M> {
  if constexpr (M != -1) return Col<M>{};
  else return Col<-1>{static_cast<Col<-1>::col>(ptrdiff_t(x))};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
rowStride(ptrdiff_t x) -> RowStride<> {
  invariant(x >= 0);
  return RowStride<-1>{static_cast<RowStride<-1>::stride>(x)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
length(ptrdiff_t x) -> Length<> {
  invariant(x >= 0);
  return Length<-1>{static_cast<Length<-1>::len>(x)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
capacity(ptrdiff_t x) -> Capacity<> {
  invariant(x >= 0);
  return Capacity<-1>{static_cast<Capacity<-1>::cap>(x)};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
capacity(size_t x) -> Capacity<> {
  invariant(x <= size_t(std::numeric_limits<ptrdiff_t>::max()));
  return capacity(ptrdiff_t(x));
}
template <std::integral I, I x>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
row(std::integral_constant<I, x>) -> Row<ptrdiff_t(x)> {
  static_assert(x >= 0);
  return {};
}
template <std::integral I, I x>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
col(std::integral_constant<I, x>) -> Col<ptrdiff_t(x)> {
  static_assert(x >= 0);
  return {};
}
template <std::integral I, I x>
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
rowStride(std::integral_constant<I, x>) -> RowStride<ptrdiff_t(x)> {
  static_assert(x >= 0);
  return {};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator+(Row<> a, Row<> b) -> Row<> {
  return {static_cast<Row<-1>::row>(ptrdiff_t(a) + ptrdiff_t(b))};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator+(Col<> a, Col<> b) -> Col<> {
  return {static_cast<Col<-1>::col>(ptrdiff_t(a) + ptrdiff_t(b))};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator-(Row<> a, Row<> b) -> Row<> {
  return {static_cast<Row<-1>::row>(ptrdiff_t(a) - ptrdiff_t(b))};
}
[[gnu::artificial, gnu::always_inline]] inline constexpr auto
operator-(Col<> a, Col<> b) -> Col<> {
  return {static_cast<Col<-1>::col>(ptrdiff_t(a) - ptrdiff_t(b))};
}

template <ptrdiff_t M> constexpr auto asrow(Length<M> len) -> Row<M> {
  if constexpr (M != -1) return {};
  else return {static_cast<Row<-1>::row>(ptrdiff_t(len))};
}
template <ptrdiff_t M> constexpr auto asrow(Col<M> len) -> Row<M> {
  if constexpr (M != -1) return {};
  else return {static_cast<Row<-1>::row>(ptrdiff_t(len))};
}
template <ptrdiff_t M> constexpr auto ascol(Length<M> len) -> Col<M> {
  if constexpr (M != -1) return {};
  else return {static_cast<Col<-1>::col>(ptrdiff_t(len))};
}
template <ptrdiff_t M> constexpr auto ascol(Row<M> len) -> Col<M> {
  if constexpr (M != -1) return {};
  else return {static_cast<Col<-1>::col>(ptrdiff_t(len))};
}
template <ptrdiff_t M>
constexpr auto asrowStride(Length<M> len) -> RowStride<M> {
  if constexpr (M != -1) return {};
  else return {static_cast<RowStride<-1>::stride>(ptrdiff_t(len))};
}

} // namespace poly::math

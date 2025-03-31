#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#include "Macros.hxx"
#ifndef USE_MODULE
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <ostream>
#include <type_traits>

#include "Math/AxisTypes.cxx"
#include "Utilities/Invariant.cxx"
#else
export module Range;

import AxisTypes;
import Invariant;
import STL;
#endif

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif
template <ptrdiff_t M>
TRIVIAL inline constexpr auto standardizeRangeBound(math::Row<M> x) {
  if constexpr (M == -1) return ptrdiff_t(x);
  else return std::integral_constant<ptrdiff_t, M>{};
}
template <ptrdiff_t M>
TRIVIAL inline constexpr auto standardizeRangeBound(math::Col<M> x) {
  if constexpr (M == -1) return ptrdiff_t(x);
  else return std::integral_constant<ptrdiff_t, M>{};
}

TRIVIAL constexpr auto standardizeRangeBound(auto x) { return x; }
TRIVIAL constexpr auto standardizeRangeBound(std::unsigned_integral auto x) {
  return ptrdiff_t(x);
}
TRIVIAL constexpr auto standardizeRangeBound(std::signed_integral auto x) {
  return ptrdiff_t(x);
}
using utils::invariant;
template <typename B, typename E> struct Range {
  [[no_unique_address]] B b_;
  [[no_unique_address]] E e_;
  [[nodiscard]] constexpr auto begin() const -> B { return b_; }
  [[nodiscard]] constexpr auto end() const -> E { return e_; }

private:
  TRIVIAL friend constexpr auto operator+(Range r, ptrdiff_t x) {
    return Range<ptrdiff_t,ptrdiff_t>(r.b_ + x, r.e_ + x);
  }
  TRIVIAL friend constexpr auto operator+(ptrdiff_t x, Range r) {
    return Range<ptrdiff_t,ptrdiff_t>(r.b_ + x, r.e_ + x);
  }
  TRIVIAL friend constexpr auto operator-(Range r, ptrdiff_t x) {
    return Range<ptrdiff_t,ptrdiff_t>(r.b_ - x, r.e_ - x);
  }
  TRIVIAL friend constexpr auto operator-(ptrdiff_t x, Range r) {
    return Range<ptrdiff_t,ptrdiff_t>(x - r.b_, x - r.e_);
  }
};
template <std::integral B, std::integral E> struct Range<B, E> {
  using value_type = std::common_type_t<B, E>;
  [[no_unique_address]] B b_;
  [[no_unique_address]] E e_;
  // wrapper that allows dereferencing
  struct Iterator {
    [[no_unique_address]] B i_;
    TRIVIAL constexpr auto operator==(E other) -> bool { return i_ == other; }
    TRIVIAL auto operator++() -> Iterator & {
      ++i_;
      return *this;
    }
    TRIVIAL auto operator++(int) -> Iterator { return Iterator{i_++}; }
    TRIVIAL auto operator--() -> Iterator & {
      --i_;
      return *this;
    }
    TRIVIAL auto operator--(int) -> Iterator { return Iterator{i_--}; }
    TRIVIAL auto operator*() -> B { return i_; }
  };
  TRIVIAL [[nodiscard]] constexpr auto begin() const -> Iterator {
    return Iterator{b_};
  }
  TRIVIAL [[nodiscard]] constexpr auto end() const -> E { return e_; }
  TRIVIAL [[nodiscard]] constexpr auto rbegin() const -> Iterator {
    return std::reverse_iterator{end()};
  }
  TRIVIAL [[nodiscard]] constexpr auto rend() const -> E {
    return std::reverse_iterator{begin()};
  }
  TRIVIAL [[nodiscard]] constexpr auto size() const { return e_ - b_; }
  friend auto operator<<(std::ostream &os, Range<B, E> r) -> std::ostream & {
    return os << "[" << r.b_ << ":" << r.e_ << ")";
  }
  template <std::integral BB, std::integral EE>
  TRIVIAL constexpr operator Range<BB, EE>() const {
    return Range<BB, EE>{static_cast<BB>(b_), static_cast<EE>(e_)};
  }
  TRIVIAL [[nodiscard]] constexpr auto view() const -> Range { return *this; }
  TRIVIAL [[nodiscard]] constexpr auto operator[](value_type i) const
    -> value_type {
    return b_ + i;
  }
  TRIVIAL [[nodiscard]] constexpr auto operator[](auto i) const {
    return i + b_;
  }

private:
  TRIVIAL friend inline constexpr auto operator+(Range r, ptrdiff_t x) {
    return Range{r.b_ + x, r.e_ + x};
  }
  TRIVIAL friend inline constexpr auto operator-(Range r, ptrdiff_t x) {
    return Range{r.b_ - x, r.e_ - x};
  }
  TRIVIAL friend inline constexpr auto operator+(ptrdiff_t x, Range r) {
    return Range{r.b_ + x, r.e_ + x};
  }
  TRIVIAL friend inline constexpr auto operator-(ptrdiff_t x, Range r) {
    return Range{x - r.b_, x - r.e_};
  }
};
template <typename B, typename E>
Range(B b, E e) -> Range<decltype(standardizeRangeBound(b)),
                         decltype(standardizeRangeBound(e))>;

TRIVIAL constexpr auto skipFirst(const auto &x) {
  auto b = x.begin();
  return Range{++b, x.end()};
}

template <typename T> struct StridedIterator {
  using value_type = std::remove_cvref_t<T>;
  T *ptr_;
  RowStride<> stride_;
  TRIVIAL constexpr auto operator==(const StridedIterator &other) const
    -> bool {
    return ptr_ == other.ptr_;
  }
  TRIVIAL constexpr auto operator!=(const StridedIterator &other) const
    -> bool {
    return ptr_ != other.ptr_;
  }
  TRIVIAL constexpr auto operator<(const StridedIterator &other) const -> bool {
    return ptr_ < other.ptr_;
  }
  TRIVIAL constexpr auto operator>(const StridedIterator &other) const -> bool {
    return ptr_ > other.ptr_;
  }
  TRIVIAL constexpr auto operator<=(const StridedIterator &other) const
    -> bool {
    return ptr_ <= other.ptr_;
  }
  TRIVIAL constexpr auto operator>=(const StridedIterator &other) const
    -> bool {
    return ptr_ >= other.ptr_;
  }
  TRIVIAL constexpr auto operator*() const -> T & { return *ptr_; }
  TRIVIAL constexpr auto operator->() const -> T * { return ptr_; }
  TRIVIAL constexpr auto operator++() -> StridedIterator & {
    ptr_ += ptrdiff_t(stride_);
    return *this;
  }
  TRIVIAL constexpr auto operator++(int) -> StridedIterator {
    auto tmp = *this;
    ptr_ += ptrdiff_t(stride_);
    return tmp;
  }
  TRIVIAL constexpr auto operator--() -> StridedIterator & {
    ptr_ -= ptrdiff_t(stride_);
    return *this;
  }
  TRIVIAL constexpr auto operator--(int) -> StridedIterator {
    auto tmp = *this;
    ptr_ -= ptrdiff_t(stride_);
    return tmp;
  }
  TRIVIAL constexpr auto operator+(ptrdiff_t x) const -> StridedIterator {
    return StridedIterator{ptr_ + x * ptrdiff_t(stride_), stride_};
  }
  TRIVIAL constexpr auto operator-(ptrdiff_t x) const -> StridedIterator {
    return StridedIterator{ptr_ - x * ptrdiff_t(stride_), stride_};
  }
  TRIVIAL constexpr auto operator+=(ptrdiff_t x) -> StridedIterator & {
    ptr_ += x * ptrdiff_t(stride_);
    return *this;
  }
  TRIVIAL constexpr auto operator-=(ptrdiff_t x) -> StridedIterator & {
    ptr_ -= x * ptrdiff_t(stride_);
    return *this;
  }
  TRIVIAL constexpr auto operator-(const StridedIterator &other) const
    -> ptrdiff_t {
    invariant(stride_ == other.stride_);
    return (ptr_ - other.ptr_) / ptrdiff_t(stride_);
  }
  TRIVIAL constexpr auto operator+(const StridedIterator &other) const
    -> ptrdiff_t {
    invariant(stride_ == other.stride_);
    return (ptr_ + other.ptr_) / ptrdiff_t(stride_);
  }
  TRIVIAL constexpr auto operator[](ptrdiff_t x) const -> T & {
    return ptr_[x * ptrdiff_t(stride_)];
  }
  TRIVIAL friend constexpr auto operator+(ptrdiff_t x,
                                          const StridedIterator &it)
    -> StridedIterator {
    return it + x;
  }
};
template <class T> StridedIterator(T *, RowStride<>) -> StridedIterator<T>;
static_assert(std::weakly_incrementable<StridedIterator<int64_t>>);
static_assert(std::input_or_output_iterator<StridedIterator<int64_t>>);
static_assert(std::indirectly_readable<StridedIterator<int64_t>>,
              "failed indirectly readable");
static_assert(std::indirectly_readable<StridedIterator<int64_t>>,
              "failed indirectly readable");
static_assert(std::output_iterator<StridedIterator<int64_t>, ptrdiff_t>,
              "failed output iterator");
static_assert(std::forward_iterator<StridedIterator<int64_t>>,
              "failed forward iterator");
static_assert(std::input_iterator<StridedIterator<int64_t>>,
              "failed input iterator");
static_assert(std::bidirectional_iterator<StridedIterator<int64_t>>,
              "failed bidirectional iterator");

static_assert(std::totally_ordered<StridedIterator<int64_t>>,
              "failed random access iterator");
static_assert(std::random_access_iterator<StridedIterator<int64_t>>,
              "failed random access iterator");
} // namespace math

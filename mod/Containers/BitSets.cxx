#ifdef USE_MODULE
module;
#else
#pragma once
#endif
#include "Macros.hxx"
#ifndef USE_MODULE
#include "Math/Array.cxx"
#include "Math/AxisTypes.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/MatrixDimensions.cxx"
#include "Utilities/ArrayPrint.cxx"
#include "Utilities/Invariant.cxx"
#include <algorithm>
#include <array>
#include <bit>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <limits>
#include <ranges>
#include <string>
#include <type_traits>
#else
export module BitSet;

import Array;
import AxisTypes;
import Invariant;
import ManagedArray;
import MatDim;
import Range;
import std;
#endif

#ifdef USE_MODULE
export namespace containers {
#else
namespace containers {
#endif
using utils::invariant;

struct EndSentinel {
  TRIVIAL [[nodiscard]] constexpr auto operator-(auto it) -> std::ptrdiff_t {
    std::ptrdiff_t i = 0;
    for (; it != EndSentinel{}; ++it, ++i) {}
    return i;
  }
  // overloaded operator== cannot be a static member function
  TRIVIAL constexpr auto operator==(EndSentinel) const -> bool { return true; }
};

template <typename T>
concept CanResize = requires(T t) { t.resize(0); };

template <std::unsigned_integral U> class BitSetIterator {
  [[no_unique_address]] const U *it_;
  [[no_unique_address]] const U *end_;
  [[no_unique_address]] U istate_;
  std::ptrdiff_t cstate0_{-1};
  std::ptrdiff_t cstate1_{0};

public:
  constexpr explicit BitSetIterator(const U *_it, const U *_end, U _istate)
    : it_{_it}, end_{_end}, istate_{_istate} {}
  using value_type = std::ptrdiff_t;
  using difference_type = std::ptrdiff_t;
  TRIVIAL constexpr auto operator*() const -> std::ptrdiff_t {
    return cstate0_ + cstate1_;
  }
  TRIVIAL constexpr auto operator++() -> BitSetIterator & {
    while (istate_ == 0) {
      if (++it_ == end_) return *this;
      istate_ = *it_;
      cstate0_ = -1;
      cstate1_ += 8 * sizeof(U);
    }
    std::ptrdiff_t tzp1 = std::countr_zero(istate_);
    cstate0_ += ++tzp1;
    istate_ = tzp1 >= 64 ? 0 : istate_ >> tzp1;
    return *this;
  }
  TRIVIAL constexpr auto operator++(int) -> BitSetIterator {
    BitSetIterator temp = *this;
    ++*this;
    return temp;
  }
  TRIVIAL constexpr auto operator==(EndSentinel) const -> bool {
    return it_ == end_ && (istate_ == 0);
  }
  TRIVIAL constexpr auto operator!=(EndSentinel) const -> bool {
    return it_ != end_ || (istate_ != 0);
  }
  TRIVIAL constexpr auto operator==(BitSetIterator j) const -> bool {
    return (it_ == j.it_) && (istate_ == j.istate_);
  }
  TRIVIAL friend constexpr auto operator==(EndSentinel,
                                           const BitSetIterator &bt) -> bool {
    return bt.it_ == bt.end_ && (bt.istate_ == 0);
  }
  TRIVIAL friend constexpr auto operator!=(EndSentinel,
                                           const BitSetIterator &bt) -> bool {
    return bt.it_ != bt.end_ || (bt.istate_ != 0);
  }
};

template <typename T>
concept BitCollection = requires(T t) {
  { std::size(t) } -> std::convertible_to<std::size_t>;
  { std::ssize(t) } -> std::convertible_to<std::ptrdiff_t>;
  { *t.begin() } -> std::convertible_to<std::uint64_t>;
};

/// A set of `std::ptrdiff_t` elements.
/// Initially constructed
template <BitCollection T = math::Vector<std::uint64_t, 1>> struct BitSet {
  using U = utils::eltype_t<T>;
  static constexpr U usize = 8 * sizeof(U);
  static constexpr U umask = usize - 1;
  static constexpr U ushift = std::countr_zero(usize);
  [[no_unique_address]] T data_{};
  // std::ptrdiff_t operator[](std::ptrdiff_t i) const {
  //     return data[i];
  // } // allow `getindex` but not `setindex`
  TRIVIAL constexpr explicit BitSet() = default;
  TRIVIAL constexpr explicit BitSet(T &&_data) : data_{std::move(_data)} {}
  TRIVIAL constexpr explicit BitSet(const T &_data) : data_{_data} {}
  TRIVIAL static constexpr auto numElementsNeeded(std::ptrdiff_t N)
    -> math::Length<> {
    return math::length(((N + usize - 1) >> ushift));
  }
  TRIVIAL constexpr explicit BitSet(std::ptrdiff_t N)
  requires(!std::is_trivially_destructible_v<T>)
    : data_{numElementsNeeded(N), 0} {}
  TRIVIAL static constexpr auto fromMask(U u) -> BitSet { return BitSet{T{u}}; }
  TRIVIAL constexpr void resizeData(std::ptrdiff_t N) {
    if constexpr (CanResize<T>) data_.resize(N);
    else invariant(N <= std::ssize(data_));
  }
  TRIVIAL constexpr void resize(std::ptrdiff_t N) {
    if constexpr (CanResize<T>) data_.resize(numElementsNeeded(N));
    else invariant(N <= std::ssize(data_) * usize);
  }
  TRIVIAL constexpr void maybeResize(std::ptrdiff_t N) {
    if constexpr (CanResize<T>) {
      math::Length<> M = numElementsNeeded(N);
      if (M > std::ssize(data_)) data_.resize(M);
    } else invariant(N <= std::ssize(data_) * std::ptrdiff_t(usize));
  }
  TRIVIAL static constexpr auto dense(std::ptrdiff_t N) -> BitSet {
    BitSet b{};
    math::Length M = numElementsNeeded(N);
    if (!M) return b;
    U maxval = std::numeric_limits<U>::max();
    if constexpr (CanResize<T>) b.data_.resizeForOverwrite(M);
    --M;
    for (std::ptrdiff_t i = 0z; i < M; ++i) b.data_[i] = maxval;
    if (std::ptrdiff_t rem = N & (usize - 1))
      b.data_[std::ptrdiff_t(M)] = (1z << rem) - 1z;
    return b;
  }
  TRIVIAL [[nodiscard]] constexpr auto maxValue() const -> std::ptrdiff_t {
    std::ptrdiff_t N = std::ssize(data_);
    return N ? ((usize * N) - std::countl_zero(data_[N - 1])) : 0;
  }
  TRIVIAL [[nodiscard]] constexpr auto count() const -> std::ptrdiff_t {
    std::ptrdiff_t c = 0;
    for (std::uint64_t x : data_) c += std::popcount(x);
    return c;
  }
  // BitSet::Iterator(std::vector<std::U> &seta)
  //     : set(seta), didx(0), offset(0), state(seta[0]), count(0) {};
  TRIVIAL [[nodiscard]] constexpr auto begin() const -> BitSetIterator<U> {
    const U *b(&*data_.begin());
    const U *e(&*data_.end());
    if (b == e) return BitSetIterator<U>{b, e, 0};
    BitSetIterator it{b, e, *b};
    return ++it;
  }
  TRIVIAL [[nodiscard]] static constexpr auto end() -> EndSentinel {
    return EndSentinel{};
  };
  TRIVIAL [[nodiscard]] constexpr auto front() const -> std::ptrdiff_t {
    for (std::ptrdiff_t i = 0; i < std::ssize(data_); ++i)
      if (data_[i]) return (usize * i) + std::countr_zero(data_[i]);
    return std::numeric_limits<std::ptrdiff_t>::max();
  }
  TRIVIAL static constexpr auto contains(math::PtrVector<U> data,
                                         std::ptrdiff_t x) -> U {
    if (data.empty()) return 0;
    std::ptrdiff_t d = x >> std::ptrdiff_t(ushift);
    U r = U(x) & umask;
    U mask = U(1) << r;
    return (data[d] & (mask));
  }
  /// Returns `true` if `i` is in the `BitSet`
  TRIVIAL [[nodiscard]] constexpr auto contains(std::ptrdiff_t i) const -> U {
    return contains(data_, i);
  }
  struct Contains {
    const T &d_;
    TRIVIAL constexpr auto operator()(std::ptrdiff_t i) const -> U {
      return contains(d_, i);
    }
  };
  TRIVIAL [[nodiscard]] constexpr auto contains() const -> Contains {
    return Contains{data_};
  }
  /// returns `true` if `x` already in the collection.
  /// Inserts it otherwise.
  TRIVIAL constexpr auto insert(std::ptrdiff_t x) -> bool {
    std::ptrdiff_t d = x >> std::ptrdiff_t(ushift);
    U r = U(x) & umask;
    U mask = U(1) << r;
    if (d >= std::ssize(data_)) resizeData(d + 1);
    bool contained = ((data_[d] & mask) != 0);
    if (!contained) data_[d] |= (mask);
    return contained;
  }
  TRIVIAL constexpr void uncheckedInsert(std::ptrdiff_t x) {
    std::ptrdiff_t d = x >> ushift;
    U r = U(x) & umask;
    U mask = U(1) << r;
    if (d >= std::ssize(data_)) resizeData(d + 1);
    data_[d] |= (mask);
  }
  // returns `true` the bitset contained `x`, i.e. if the
  // removal was succesful.
  TRIVIAL constexpr auto remove(std::ptrdiff_t x) -> bool {
    std::ptrdiff_t d = x >> ushift;
    U r = U(x) & umask;
    U mask = U(1) << r;
    bool contained = ((data_[d] & mask) != 0);
    if (contained) data_[d] &= (~mask);
    return contained;
  }
  TRIVIAL static constexpr void set(U &d, std::ptrdiff_t r, bool b) {
    U mask = U(1) << r;
    if (b == ((d & mask) != 0)) return;
    if (b) d |= mask;
    else d &= (~mask);
  }
  TRIVIAL static constexpr void set(math::MutPtrVector<U> data,
                                    std::ptrdiff_t x, bool b) {
    std::ptrdiff_t d = x >> ushift;
    U r = U(x) & umask;
    set(data[d], r, b);
  }

  class Reference {
    [[no_unique_address]] math::MutPtrVector<U> data_;
    [[no_unique_address]] std::ptrdiff_t i_;

  public:
    TRIVIAL constexpr explicit Reference(math::MutPtrVector<U> dd,
                                         std::ptrdiff_t ii)
      : data_(dd), i_(ii) {}
    TRIVIAL constexpr operator bool() const { return contains(data_, i_); }
    TRIVIAL constexpr auto operator=(bool b) -> Reference & {
      BitSet::set(data_, i_, b);
      return *this;
    }
  };

  TRIVIAL constexpr auto operator[](std::ptrdiff_t i) const -> bool {
    return contains(data_, i);
  }
  TRIVIAL constexpr auto operator[](std::ptrdiff_t i) -> Reference {
    maybeResize(i + 1);
    math::MutPtrVector<U> d{data_};
    return Reference{d, i};
  }
  TRIVIAL [[nodiscard]] constexpr auto size() const -> std::ptrdiff_t {
    std::ptrdiff_t s = 0;
    for (auto u : data_) s += std::popcount(u);
    return s;
  }
  TRIVIAL [[nodiscard]] constexpr auto empty() const -> bool {
    return std::ranges::all_of(data_, [](auto u) { return u == 0; });
  }
  TRIVIAL [[nodiscard]] constexpr auto any() const -> bool {
    return std::ranges::any_of(data_, [](auto u) { return u != 0; });
    // for (auto u : data)
    //   if (u) return true;
    // return false;
  }
  TRIVIAL constexpr void setUnion(const BitSet &bs) {
    std::ptrdiff_t O = std::ssize(bs.data_), N = std::ssize(data_);
    if (O > N) resizeData(O);
    for (std::ptrdiff_t i = 0; i < O; ++i) {
      U d = data_[i] | bs.data_[i];
      data_[i] = d;
    }
  }
  TRIVIAL constexpr auto operator&=(const BitSet &bs) -> BitSet & {
    if (std::ssize(bs.data_) < std::ssize(data_))
      resizeData(std::ssize(bs.data_));
    for (std::ptrdiff_t i = 0; i < std::ssize(data_); ++i)
      data_[i] &= bs.data_[i];
    return *this;
  }
  // &!
  TRIVIAL constexpr auto operator-=(const BitSet &bs) -> BitSet & {
    if (std::ssize(bs.data_) < std::ssize(data_))
      resizeData(std::ssize(bs.data_));
    for (std::ptrdiff_t i = 0; i < std::ssize(data_); ++i)
      data_[i] &= (~bs.data_[i]);
    return *this;
  }
  TRIVIAL constexpr auto operator|=(const BitSet &bs) -> BitSet & {
    if (std::ssize(bs.data_) > std::ssize(data_))
      resizeData(std::ssize(bs.data_));
    for (std::ptrdiff_t i = 0; i < std::ssize(bs.data_); ++i)
      data_[i] |= bs.data_[i];
    return *this;
  }
  TRIVIAL constexpr auto operator&(const BitSet &bs) const -> BitSet {
    BitSet r = *this;
    return r &= bs;
  }
  TRIVIAL constexpr auto operator|(const BitSet &bs) const -> BitSet {
    BitSet r = *this;
    return r |= bs;
  }
  TRIVIAL constexpr auto operator==(const BitSet &bs) const -> bool {
    return data_ == bs.data_;
  }
  TRIVIAL constexpr auto operator~() const -> BitSet {
    BitSet r = *this;
    for (U &u : r.data_) u = ~u;
    return r;
  }
  // Ranks higher elements as more important, thus iterating
  // backwards.
  TRIVIAL constexpr auto operator<=>(const BitSet &other) const
    -> std::strong_ordering {
    std::ptrdiff_t ntd = data_.size(), nod = other.data_.size();
    if (ntd != nod) {
      bool larger = ntd > nod;
      std::ptrdiff_t l = std::min(ntd, nod), L = std::max(ntd, nod);
      const T &d = larger ? data_ : other.data_;
      // The other is effectively all `0`, thus we may as well iterate forwards.
      for (std::ptrdiff_t i = l; i < L; ++i)
        if (d[i])
          return larger ? std::strong_ordering::greater
                        : std::strong_ordering::less;
      ntd = l;
    }
    for (std::ptrdiff_t i = ntd; i--;)
      if (auto cmp = data_[i] <=> other.data_[i]; cmp != 0) return cmp;
    return std::strong_ordering::equal;
  }

  void print() const {
    utils::print("BitSet[");
    auto it = begin();
    constexpr EndSentinel e = BitSet::end();
    if (it != e) {
      utils::print(*(it++));
      for (; it != e; ++it) utils::print(", ", *it);
    }
    utils::print("]");
  }
  TRIVIAL constexpr void clear() {
    std::fill_n(data_.begin(), std::ssize(data_), 0);
  }
  TRIVIAL [[nodiscard]] constexpr auto isEmpty() const -> bool {
    return std::ranges::all_of(data_, [](auto u) { return u == 0; });
    // for (auto u : data)
    //   if (u) return false;
    // return true;
  }
  TRIVIAL constexpr auto findFirstZero() -> std::ptrdiff_t {
    std::ptrdiff_t offset = 0;
    for (U x : data_) {
      U c = std::countr_one(x);
      offset += c;
      if (c != usize) return offset;
    }
    return usize;
  }
};

template <unsigned N>
using FixedSizeBitSet = BitSet<std::array<std::uint64_t, N>>;
// BitSet with length 64
using BitSet64 = FixedSizeBitSet<1>;
static_assert(std::is_trivially_destructible_v<BitSet64>);
static_assert(std::is_trivially_destructible_v<FixedSizeBitSet<2>>);
// static_assert(std::input_or_output_iterator<
//               decltype(std::declval<FixedSizeBitSet<2>>().begin())>);
static_assert(std::ranges::range<FixedSizeBitSet<2>>);

template <typename T, typename B = BitSet<>> struct BitSliceView {
  [[no_unique_address]] math::MutPtrVector<T> a_;
  [[no_unique_address]] const B &i_;
  struct Iterator {
    [[no_unique_address]] math::MutPtrVector<T> a_;
    [[no_unique_address]] BitSetIterator<std::uint64_t> it_;
    TRIVIAL constexpr auto operator==(EndSentinel) const -> bool {
      return it_ == EndSentinel{};
    }
    TRIVIAL constexpr auto operator++() -> Iterator & {
      ++it_;
      return *this;
    }
    TRIVIAL constexpr auto operator++(int) -> Iterator {
      Iterator temp = *this;
      ++it_;
      return temp;
    }
    TRIVIAL constexpr auto operator*() -> T & { return a_[*it_]; }
    TRIVIAL constexpr auto operator*() const -> const T & { return a_[*it_]; }
    TRIVIAL constexpr auto operator->() -> T * { return &a_[*it_]; }
    TRIVIAL constexpr auto operator->() const -> const T * { return &a_[*it_]; }
  };
  TRIVIAL constexpr auto begin() -> Iterator { return {a_, i_.begin()}; }
  struct ConstIterator {
    [[no_unique_address]] math::PtrVector<T> a_;
    [[no_unique_address]] BitSetIterator<std::uint64_t> it_;
    TRIVIAL constexpr auto operator==(EndSentinel) const -> bool {
      return it_ == EndSentinel{};
    }
    TRIVIAL constexpr auto operator==(ConstIterator c) const -> bool {
      return (it_ == c.it_) && (a_.data() == c.a_.data());
    }
    TRIVIAL constexpr auto operator++() -> ConstIterator & {
      ++it_;
      return *this;
    }
    TRIVIAL constexpr auto operator++(int) -> ConstIterator {
      ConstIterator temp = *this;
      ++it_;
      return temp;
    }
    TRIVIAL constexpr auto operator*() const -> const T & { return a_[*it_]; }
    TRIVIAL constexpr auto operator->() const -> const T * { return &a_[*it_]; }
  };
  TRIVIAL [[nodiscard]] constexpr auto begin() const -> ConstIterator {
    return {a_, i_.begin()};
  }
  TRIVIAL [[nodiscard]] constexpr auto end() const -> EndSentinel { return {}; }
  TRIVIAL [[nodiscard]] constexpr auto size() const -> std::ptrdiff_t {
    return i_.size();
  }
  TRIVIAL [[nodiscard]] friend constexpr auto operator-(EndSentinel, Iterator v)
    -> std::ptrdiff_t {
    return EndSentinel{} - v.it_;
  }
  TRIVIAL [[nodiscard]] friend constexpr auto operator-(EndSentinel,
                                                        ConstIterator v)
    -> std::ptrdiff_t {
    return EndSentinel{} - v.it_;
  }
};
template <typename T, typename B>
BitSliceView(math::MutPtrVector<T>, const B &) -> BitSliceView<T, B>;

static_assert(std::movable<BitSliceView<std::int64_t>::Iterator>);
static_assert(std::movable<BitSliceView<std::int64_t>::ConstIterator>);
} // namespace containers

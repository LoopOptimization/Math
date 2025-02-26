#ifdef USE_MODULE
module;
#else
#pragma once
#endif
#ifndef USE_MODULE
#include <cstddef>
#include <memory>

#include "Alloc/Arena.cxx"
#include "Utilities/Invariant.cxx"
#else
export module UnrolledList;

import Arena;
import Invariant;
import STL;
#endif

#ifdef USE_MODULE
export namespace containers {
#else
namespace containers {
#endif
using utils::invariant;
template <typename T> class UList {
  T data_[6]; // NOLINT(modernize-avoid-c-arrays)
  ptrdiff_t count_{0};
  UList<T> *next_{nullptr};

public:
  constexpr UList() = default;
  constexpr UList(T t) : count_(1) { data_[0] = t; }
  constexpr UList(T t, UList *n) : count_(1), next_(n) {
    std::construct_at(data_, t);
  }
  // constexpr UList(T t, UList *n) : count_(1), next_(n) { data_[0] = t; }
  constexpr UList(const UList &other) = default;
  [[nodiscard]] constexpr auto getHeadCount() const -> ptrdiff_t {
    return count_;
  }
  constexpr void forEach(const auto &f) {
    invariant(count_ <= std::ssize(data_));
    for (auto *L = this; L != nullptr; L = L->next_)
      for (ptrdiff_t i = 0, N = L->count_; i < N; i++) f(L->data_[i]);
  }
  constexpr void forEachRev(const auto &f) {
    invariant(count_ <= std::ssize(data_));
    for (auto *L = this; L != nullptr; L = L->next_)
      for (ptrdiff_t i = L->count_; i;) f(L->data_[--i]);
  }
  constexpr void forEachStack(const auto &f) {
    invariant(count_ <= std::ssize(data_));
    // the motivation of this implementation is that we use this to
    // deallocate the list, which may contain pointers that themselves
    // allocated this.
    UList<T> C{*this};
    while (true) {
      for (ptrdiff_t i = 0, N = C.count_; i < N; i++) f(C.data_[i]);
      if (C.next_ == nullptr) return;
      C = *C.next_;
    }
  }
  constexpr void forEachNoRecurse(const auto &f) {
    invariant(count_ <= std::ssize(data_));
    for (ptrdiff_t i = 0; i < count_; i++) f(data_[i]);
  }
  constexpr auto reduce(auto init, const auto &f) const {
    invariant(count_ <= std::ssize(data_));
    decltype(f(init, std::declval<T>())) acc = init;
    for (auto *L = this; L != nullptr; L = L->next_)
      for (ptrdiff_t i = 0, N = L->count_; i < N; i++)
        acc = f(acc, L->data_[i]);
    return acc;
  }
  constexpr auto transform_reduce(auto init, const auto &f) {
    invariant(count_ <= std::ssize(data_));
    decltype(f(init, std::declval<T &>())) acc = init;
    for (auto *L = this; L != nullptr; L = L->next_)
      for (ptrdiff_t i = 0, N = L->count_; i < N; i++)
        acc = f(acc, L->data_[i]);
    return acc;
  }
  constexpr void pushHasCapacity(T t) {
    invariant(count_ < std::ssize(data_));
    data_[count_++] = t;
  }
  /// unordered push
  template <class A>
  [[nodiscard]] constexpr auto push(A &alloc, T t) -> UList * {
    invariant(count_ <= std::ssize(data_));
    if (!isFull()) {
      data_[count_++] = t;
      return this;
    }
    UList<T> *other = alloc.allocate(1);
    std::construct_at(other, t, this);
    return other;
  }

  /// ordered push
  template <class A> constexpr void push_ordered(A &alloc, T t) {
    invariant(count_ <= std::ssize(data_));
    if (!isFull()) {
      data_[count_++] = t;
      return;
    }
    if (next_ == nullptr) {
      next_ = alloc.allocate(1);
      std::construct_at(next_, t);
    } else next_->push_ordered(alloc, t);
  }
  [[nodiscard]] constexpr auto contains(T t) const -> bool {
    invariant(count_ <= std::ssize(data_));
    for (const UList *L = this; L; L = L->getNext())
      for (size_t i = 0, N = L->getHeadCount(); i < N; ++i)
        if (data_[i] == t) return true;
    return false;
  }
  /// pushUnique(allocator, t)
  /// pushes `t` if it is unique
  template <class A>
  [[nodiscard]] constexpr auto pushUnique(A &alloc, T t) -> UList * {
    if (contains(t)) return this;
    return push(alloc, t);
  }

  // too dangerous
  // [[nodiscard]] constexpr auto push(T t) -> UList * {
  //   std::allocator<UList<T>> alloc;
  //   return push(alloc, t);
  // }
  // /// ordered push
  // constexpr void push_ordered(T t) {
  //   std::allocator<UList<T>> alloc;
  //   push_ordered(alloc, t);
  // }
  [[nodiscard]] constexpr auto push(alloc::Arena<> *alloc, T t) -> UList * {
    invariant(count_ <= std::ssize(data_));
    if (isFull()) return alloc->create<UList<T>>(t, this);
    data_[count_++] = t;
    return this;
  };
  constexpr void push_ordered(alloc::Arena<> *alloc, T t) {
    invariant(count_ <= std::ssize(data_));
    if (!isFull()) {
      data_[count_++] = t;
      return;
    }
    if (next_ == nullptr) next_ = alloc->create<UList<T>>(t);
    else next_->push_ordered(alloc, t);
  }
  constexpr auto copy(alloc::Arena<> *alloc) const -> UList * {
    UList<T> *L = alloc->create<UList<T>>();
    L->count_ = count_;
    std::copy(std::begin(data_), std::end(data_), std::begin(L->data_));
    if (next_) L->next_ = next_->copy(alloc);
    return L;
  }
  /// erase
  /// behavior is undefined if `x` doesn't point to this node
  constexpr void erase(T *x) {
    invariant(count_ <= std::ssize(data_));
    for (auto i = x, e = data_ + --count_; i != e; ++i) *i = *(i + 1);
  }
  /// eraseUnordered
  /// behavior is undefined if `x` doesn't point to this node
  constexpr void eraseUnordered(T *x) {
    invariant(count_ <= std::ssize(data_));
    *x = data_[--count_];
  }
  constexpr auto searchHead(T x) -> T * {
    for (auto *d = data_, *e = d + count_; d != e; ++d)
      if (*d == x) return d;
    return nullptr;
  }
  //
  constexpr void eraseUnordered(T x) {
    invariant(count_ || next_ != nullptr);
    if (!count_) next_->eraseUnordered(x);
    if (T *p = searchHead(x)) return eraseUnordered(p);
    // not in head -> search next until we find it;
    // move last here there.
    invariant(next_ != nullptr);
    next_->swapWith(x, std::move(data_[--count_]));
  }
  // search for `x`, swap with `y`.
  void swapWith(T x, T y) {
    for (UList *L = this; L; L->getNext()) {
      if (T *p = searchHead(x)) {
        *p = y;
        return;
      }
    }
  }
  [[nodiscard]] constexpr auto isFull() const -> bool {
    return count_ == std::ssize(data_);
  }
  [[nodiscard]] constexpr auto getNext() const -> UList * { return next_; }
  constexpr void clear() {
    invariant(count_ <= std::ssize(data_));
    count_ = 0;
    next_ = nullptr;
  }
  constexpr auto front() -> T & {
    invariant(count_ > 0);
    return data_[0];
  }
  constexpr auto only() -> T & {
    invariant(count_ == 1);
    return data_[0];
  }
  [[nodiscard]] constexpr auto front() const -> const T & {
    invariant(count_ > 0);
    return data_[0];
  }
  [[nodiscard]] constexpr auto only() const -> const T & {
    invariant(count_ == 1);
    return data_[0];
  }
  constexpr void append(UList *L) {
    UList *N = this;
    while (N->next_ != nullptr) N = N->next_;
    N->next_ = L;
  }
  [[nodiscard]] constexpr auto empty() const -> bool { return count_ == 0; }
  [[nodiscard]] constexpr auto operator==(const UList &other) const -> bool {
    if (count_ != other.count_) return false;
    for (ptrdiff_t i = 0; i < count_; i++)
      if (data_[i] != other.data_[i]) return false;
    if (next_ == nullptr && other.getNext() == nullptr) return true;
    if (next_ == nullptr || other.getNext() == nullptr) return false;
    return *next_ == *other.getNext();
  }
  constexpr auto operator[](ptrdiff_t i) -> T & {
    return (i < count_) ? data_[i] : next_->operator[](i - count_);
  }
  constexpr auto operator[](ptrdiff_t i) const -> const T & {
    return (i < count_) ? data_[i] : next_->operator[](i - count_);
  }
  constexpr auto operator=(const UList &other) -> UList & = default;
  struct End {};
  struct MutIterator {

    UList *list;
    ptrdiff_t index;
    constexpr auto operator==(End) const -> bool { return list == nullptr; }
    constexpr auto operator==(MutIterator other) const -> bool {
      return list == other.list && index == other.index;
    }
    constexpr auto operator++() -> MutIterator & {
      invariant(list != nullptr);
      if (++index == list->count_) {
        list = list->next_;
        index = 0;
      }
      return *this;
    }
    constexpr auto operator++(int) -> MutIterator {
      invariant(list != nullptr);
      auto ret = *this;
      ++*this;
      return ret;
    }
    constexpr auto operator*() -> T & {
      invariant(list != nullptr);
      return list->data_[index];
    }
    constexpr auto operator->() -> T * { return &**this; }
  };
  struct Iterator {
    const UList *list;
    ptrdiff_t index;
    constexpr auto operator==(End) const -> bool { return list == nullptr; }
    constexpr auto operator==(Iterator other) const -> bool {
      return list == other.list && index == other.index;
    }
    constexpr auto operator++() -> Iterator & {
      invariant(list != nullptr);
      if (++index == list->count_) {
        list = list->next_;
        index = 0;
      }
      return *this;
    }
    constexpr auto operator++(int) -> Iterator {
      invariant(list != nullptr);
      auto ret = *this;
      ++*this;
      return ret;
    }
    constexpr auto operator*() -> T {
      invariant(list != nullptr);
      return list->data_[index];
    }
    constexpr auto operator->() -> const T * { return &**this; }
  };
  [[nodiscard]] constexpr auto begin() const -> Iterator { return {this, 0}; }
  [[nodiscard]] constexpr auto begin() -> MutIterator { return {this, 0}; }
  static constexpr auto end() -> End { return {}; }
  constexpr auto dbegin() -> T * { return data_; }
  constexpr auto dend() -> T * { return data_ + count_; }
};

} // namespace containers

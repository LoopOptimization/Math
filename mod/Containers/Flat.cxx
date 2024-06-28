module;
#include <cstddef>
#include <memory>
#include <type_traits>

export module Flat;

import Allocator;

export namespace containers {
template <typename T> [[gsl::Owner(T)]] struct Flat {

  explicit constexpr auto Flat(ptrdiff_t len)
    : ptr_{alloc::Mallocator<T>{}.allocate(len)}, len_{len} {
    std::uninitialized_default_construct_n(ptr_, len_);
  };

  constexpr Flat(const Flat &other) : Flat(other.size()) {
    std::uninitialized_copy_n(other.data(), len_, data());
  };
  constexpr auto operator=(const Flat &other) {
    if (len_ != other.size()) {
      maybeDeallocate();
      ptr_ = alloc::Mallocator<T>{}.allocate(len);
      len_ = other.size();
      if constexpr (!std::is_trivially_default_constructible_v<T>) {
        std::uninitialized_copy_n(other.data(), len_, ptr_);
        return;
      }
    }
    std::copy_n(other.data(), len_, data());
  };
  constexpr Flat(Flat &&other) : ptr_{other.data()}, len_{other.size()} {
    other.ptr_ = nullptr;
    other.len_ = 0;
  };
  constexpr auto operator=(Flat &&other) {
    ptr_ = other.ptr_;
    len_ = other.size();
    other.ptr_ = nullptr;
    other.len_ = 0;
  }

  constexpr ~Flat() { maybeDeallocate(); }
  constexpr auto data() -> T * { return ptr_; }
  constexpr auto data() const -> const T * { return ptr_; }
  constexpr auto size() const -> ptrdiff_t {
    invariant(size() >= 0);
    return len_;
  }
  constexpr auto operator[](ptrdiff_t i) -> T * {
    invariant((i >= 0) && i < size());
    return ptr_[i];
  }
  constexpr auto operator[](ptrdiff_t i) const -> const T * {
    invariant((i >= 0) && i < size());
    return ptr_[i];
  }

private:
  constexpr void maybeDeallocate() {
    if (!ptr_) return;
    if constexpr (!std::is_trivially_destructible_v<T>)
      std::destroy_n(ptr_, len_);
    alloc::Mallocator<T>{}.deallocate(ptr_, len_);
  }

  T *ptr_{nullptr};
  ptrdiff_t len_{0};
};
} // namespace containers


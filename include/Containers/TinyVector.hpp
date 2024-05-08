#pragma once

#include "Containers/Storage.hpp"
#include "Math/AxisTypes.hpp"
#include "Utilities/Invariant.hpp"
#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <memory>
#include <type_traits>

namespace poly::containers {
using utils::invariant;

template <class T, size_t N, std::signed_integral L = ptrdiff_t>
class TinyVector {
  static_assert(N > 0);
  static_assert(std::numeric_limits<ptrdiff_t>::max() >= N);
  using Length = math::Length<-1, L>;
  Storage<T, N> data;
  Length len{};

public:
  using value_type = T;
  constexpr TinyVector(){}; // NOLINT (modernize-use-equals-default)
  constexpr TinyVector(const std::initializer_list<T> &list)
    : len{L(list.size())} {
    invariant(list.size() <= N);
    if constexpr (Storage<T, N>::trivial)
      std::copy_n(list.begin(), len, data.data());
    else std::uninitialized_move_n(list.begin(), len, data.data());
  }
  constexpr TinyVector(T t) : len{1} {
    std::construct_at(data.data(), std::move(t));
  }

  constexpr auto
  operator=(const std::initializer_list<T> &list) -> TinyVector & {
    invariant(list.size() <= ptrdiff_t(N));
    ptrdiff_t slen = std::ssize(list), oldLen = ptrdiff_t(len);
    len = math::length(slen);
    auto I = list.begin();
    if constexpr (Storage<T, N>::trivial) {
      // implicit lifetime type
      std::copy_n(I, slen, data.data());
    } else {
      ptrdiff_t J = std::min(len, ptrdiff_t(oldLen));
      // old values exist
      std::move(I, I + J, data.data());
      if (oldLen < len) {
        std::uninitialized_move_n(I + oldLen, slen - oldLen,
                                  data.data() + oldLen);
      } else if constexpr (!std::is_trivially_destructible_v<T>)
        if (oldLen > len) std::destroy_n(data.data() + len, oldLen - len);
    }
    return *this;
  }
  constexpr auto operator[](ptrdiff_t i) -> T & {
    invariant(i < len);
    return data.data()[i];
  }
  constexpr auto operator[](ptrdiff_t i) const -> const T & {
    invariant(i < len);
    return data.data()[i];
  }
  constexpr auto back() -> T & {
    invariant(len > 0);
    return data.data()[size() - 1z];
  }
  constexpr auto back() const -> const T & {
    invariant(len > 0);
    return data.data()[size() - 1z];
  }
  constexpr auto front() -> T & {
    invariant(len > 0);
    return data.data()[0z];
  }
  constexpr auto front() const -> const T & {
    invariant(len > 0);
    return data.data()[0z];
  }
  constexpr void push_back(const T &t) {
    invariant(len < math::length(L(N)));
    std::construct_at(data.data() + ptrdiff_t(len++), t);
  }
  constexpr void push_back(T &&t) {
    invariant(len < ptrdiff_t(N));
    std::construct_at(data.data() + ptrdiff_t(L(len++)), std::move(t));
  }
  template <class... Args> constexpr auto emplace_back(Args &&...args) -> T & {
    invariant(len < ptrdiff_t(N));
    return *std::construct_at(data.data() + ptrdiff_t(L(len++)),
                              std::forward<Args>(args)...);
  }
  constexpr void pop_back() {
    invariant(len > 0);
    --len;
    if constexpr (!std::is_trivially_destructible_v<T>)
      std::destroy_at(data.data() + size());
  }
  [[nodiscard]] constexpr auto pop_back_val() -> T {
    invariant(len > 0);
    return std::move(data.data()[ptrdiff_t(L(--len))]);
  }
  [[nodiscard]] constexpr auto size() const -> ptrdiff_t {
    auto l = ptrdiff_t(L(len));
    invariant(l <= ptrdiff_t(N));
    return l;
  }
  [[nodiscard]] constexpr auto empty() const -> bool { return len == 0z; }
  constexpr void clear() {
    len = 0;
    if constexpr (!std::is_trivially_destructible_v<T>)
      std::destroy_n(data.data(), ptrdiff_t(len));
  }

  constexpr auto begin() -> T * { return data.data(); }
  constexpr auto begin() const -> const T * { return data.data(); }
  constexpr auto end() -> T * { return data.data() + size(); }
  constexpr auto end() const -> const T * { return data.data() + size(); }
  constexpr void resize(L new_size) {
    // initialize new data
    for (ptrdiff_t i = size(); i < new_size; ++i)
      std::construct_at(data.data() + i);
    len = math::length(new_size);
  }
  constexpr void reserve(L space) {
    invariant(space >= 0);
    invariant(size_t(space) <= N);
  }
  constexpr ~TinyVector()
  requires(std::is_trivially_destructible_v<T>)
  = default;
  constexpr ~TinyVector()
  requires(!std::is_trivially_destructible_v<T>)
  {
    std::destroy_n(data.data(), len);
  }
  friend inline auto operator<<(std::ostream &os,
                                const TinyVector &x) -> std::ostream & {
    os << "[";
    if constexpr (std::same_as<T, int8_t> || std::same_as<T, uint8_t>) {
      if (!x.empty()) os << int(x[0]);
      for (L i = 1; i < x.size(); ++i) os << ", " << int(x[i]);
    } else {
      if (!x.empty()) os << x[0];
      for (L i = 1; i < x.size(); ++i) os << ", " << x[i];
    }
    return os << "]";
  }
};
} // namespace poly::containers

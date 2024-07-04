module;
#include <concepts>
#include <cstddef>

export module CompressReference;

import TypeCompression;

export namespace utils {
template <typename T> struct Reference {
  using C = utils::compressed_t<T>;
  static_assert(!std::same_as<C, T>);
  C *c;
  constexpr operator T() const { return T::decompress(c); }
  constexpr auto view() const -> T { return T::decompress(c); }
  // constexpr operator T &() const { return *t; }
  constexpr auto operator=(const T &t) -> Reference & {
    t.compress(c);
    return *this;
  }
  constexpr auto operator=(const C &x) -> Reference & {
    // Probably shouldn't be needed?
    // TODO: try removing this method
    *c = x;
    return *this;
  }
  constexpr auto operator==(const T &t) const -> bool {
    return T::decompress(c) == t;
  }
  [[gnu::always_inline]] constexpr auto operator+=(const auto &x) {
    T y{T::decompress(c)};
    y += x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator-=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y -= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator*=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y *= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator/=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y /= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator%=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y %= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator<<=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y <<= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator>>=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y >>= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator&=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y &= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator^=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y ^= x;
    y.compress(c);
    return y;
    // return *this;
  }
  [[gnu::always_inline]] constexpr auto operator|=(const auto &x) {
    // -> Reference & {
    T y{T::decompress(c)};
    y |= x;
    y.compress(c);
    return y;
    // return *this;
  }

  constexpr auto operator[](auto i) -> decltype(auto) {
    return c->operator[](i);
  }
  constexpr auto operator[](auto i) const -> decltype(auto) {
    return c->operator[](i);
  }
  constexpr auto operator[](auto i, auto j) -> decltype(auto) {
    return c->operator[](i, j);
  }
  constexpr auto operator[](auto i, auto j) const -> decltype(auto) {
    return c->operator[](i, j);
  }

private:
  // TODO: are these really needed / can we rely on implicit conversion?
  friend constexpr auto operator+(Reference x, const auto &y) {
    return T::decompress(x.c) + y;
  }
  friend constexpr auto operator-(Reference x, const auto &y) {
    return T::decompress(x.c) - y;
  }
  friend constexpr auto operator*(Reference x, const auto &y) {
    return T::decompress(x.c) * y;
  }
  friend constexpr auto operator/(Reference x, const auto &y) {
    return T::decompress(x.c) / y;
  }
  friend constexpr auto operator%(Reference x, const auto &y) {
    return T::decompress(x.c) % y;
  }
  friend constexpr auto operator>>(Reference x, const auto &y) {
    return T::decompress(x.c) >> y;
  }
  friend constexpr auto operator<<(Reference x, const auto &y) {
    return T::decompress(x.c) << y;
  }
  friend constexpr auto operator&(Reference x, const auto &y) {
    return T::decompress(x.c) & y;
  }
  friend constexpr auto operator^(Reference x, const auto &y) {
    return T::decompress(x.c) ^ y;
  }
  friend constexpr auto operator|(Reference x, const auto &y) {
    return T::decompress(x.c) | y;
  }
  friend constexpr auto operator>(Reference x, const auto &y) {
    return T::decompress(x.c) > y;
  }
  friend constexpr auto operator>=(Reference x, const auto &y) {
    return T::decompress(x.c) >= y;
  }
  friend constexpr auto operator<(Reference x, const auto &y) {
    return T::decompress(x.c) < y;
  }
  friend constexpr auto operator<=(Reference x, const auto &y) {
    return T::decompress(x.c) <= y;
  }
  friend constexpr auto operator==(Reference x, const auto &y) {
    return T::decompress(x.c) == y;
  }
  friend constexpr auto operator!=(Reference x, const auto &y) {
    return T::decompress(x.c) != y;
  }
  friend constexpr auto operator+(const auto &x, Reference y) {
    return x + T::decompress(y.c);
  }
  friend constexpr auto operator-(const auto &x, Reference y) {
    return x - T::decompress(y.c);
  }
  friend constexpr auto operator*(const auto &x, Reference y) {
    return x * T::decompress(y.c);
  }
  friend constexpr auto operator/(const auto &x, Reference y) {
    return x / T::decompress(y.c);
  }
  friend constexpr auto operator%(const auto &x, Reference y) {
    return x % T::decompress(y.c);
  }
  friend constexpr auto operator>>(const auto &x, Reference y) {
    return x >> T::decompress(y.c);
  }
  friend constexpr auto operator<<(const auto &x, Reference y) {
    return x << T::decompress(y.c);
  }
  friend constexpr auto operator&(const auto &x, Reference y) {
    return x & T::decompress(y.c);
  }
  friend constexpr auto operator^(const auto &x, Reference y) {
    return x ^ T::decompress(y.c);
  }
  friend constexpr auto operator|(const auto &x, Reference y) {
    return x | T::decompress(y.c);
  }
  friend constexpr auto operator<(const auto &x, Reference y) {
    return x < T::decompress(y.c);
  }
  friend constexpr auto operator<=(const auto &x, Reference y) {
    return x <= T::decompress(y.c);
  }
  friend constexpr auto operator>(const auto &x, Reference y) {
    return x > T::decompress(y.c);
  }
  friend constexpr auto operator>=(const auto &x, Reference y) {
    return x >= T::decompress(y.c);
  }
  friend constexpr auto operator==(const auto &x, Reference y) {
    return x == T::decompress(y.c);
  }
  friend constexpr auto operator!=(const auto &x, Reference y) {
    return x != T::decompress(y.c);
  }
  friend constexpr void swap(Reference x, Reference y) {
    std::swap(*x.c, *y.c);
  }
  friend constexpr auto value(Reference x) {
    if constexpr (requires(C *c) { c->value(); }) return value(x.c->value());
    else return value(*x.c);
  }
};

template <typename T>
[[gnu::always_inline]] constexpr auto ref(T *p, ptrdiff_t i) -> T & {
  return p[i];
}
template <typename T>
[[gnu::always_inline]] constexpr auto ref(const T *p,
                                          ptrdiff_t i) -> const T & {
  return p[i];
}
template <utils::Compressible T>
[[gnu::always_inline]] constexpr auto ref(utils::compressed_t<T> *p,
                                          ptrdiff_t i) -> Reference<T> {
  return Reference<T>{p + i};
}

template <utils::Compressible T>
[[gnu::always_inline]] constexpr auto ref(const utils::compressed_t<T> *p,
                                          ptrdiff_t i) -> T {
  return T::decompress(p + i);
}

} // namespace utils

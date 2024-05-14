#pragma once
#include <cstddef>
#include <type_traits>
#include <utility>

namespace poly::containers {

template <typename To, typename From>
concept ConvertibleFrom = std::is_convertible_v<From, To>;

template <class F, class S> struct Pair {
  [[no_unique_address]] F first;
  [[no_unique_address]] S second;
  template <ConvertibleFrom<F> A, ConvertibleFrom<S> B>
  constexpr operator Pair<A, B>() {
    return {A{first}, B{second}};
  }
  template <size_t I> constexpr auto get() -> auto & {
    if constexpr (I == 0) return first;
    else return second;
  }
  template <size_t I> constexpr auto get() const -> const auto & {
    if constexpr (I == 0) return first;
    else return second;
  }
  // template <typename T, typename U>
  // constexpr auto operator=(Pair<T, U> x)
  //   -> Pair &requires(std::assignable_from<F, T> &&std::assignable_from<S,
  //   U>) { first = x.first; second = x.second; return *this;
  // }
  friend void print_obj(std::ostream &os, const containers::Pair<F, S> &x) {
    os << "(" << x.first << ", " << x.second << ")";
  };
};
} // namespace poly::containers

template <typename F, typename S>
struct std::tuple_size<poly::containers::Pair<F, S>>
  : public std::integral_constant<size_t, 2> {};

template <typename F, typename S>
struct std::tuple_element<0, poly::containers::Pair<F, S>> {
  using type = F;
};
template <typename F, typename S>
struct std::tuple_element<1, poly::containers::Pair<F, S>> {
  using type = S;
};


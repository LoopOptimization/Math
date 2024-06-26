#pragma once

#include <type_traits>

namespace poly::utils {
template <typename T>
concept TriviallyCopyable =
  std::is_trivially_copyable_v<T> && std::is_trivially_destructible_v<T>;
template <typename T> struct InParameter {
  using type = const T &;
};
template <TriviallyCopyable T> struct InParameter<T> {
  using type = T;
};

/// This can be used like
/// auto foo_impl(inparam_t<T> x, ...);
/// template <typename T>
/// [[gnu::always_inline]] inline auto foo(const T& x){
///  return foo_impl<T>(x); // inparam_t blocks deduction
/// }
template <typename T> using inparam_t = typename InParameter<T>::type;

} // namespace poly::utils

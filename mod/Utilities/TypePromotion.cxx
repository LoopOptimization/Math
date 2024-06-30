module;

#include <concepts>
#include <type_traits>

export module TypePromotion;

export namespace utils {
template <typename T>
concept HasEltype = requires(T) {
  typename std::remove_reference_t<T>::value_type;
  // std::is_scalar_v<typename std::remove_reference_t<T>::value_type>;
};
} // namespace utils

template <typename A> struct GetEltype {
  using value_type = A;
};
template <utils::HasEltype A> struct GetEltype<A> {
  using value_type = typename A::value_type;
};

export namespace utils {

template <typename T>
using eltype_t = typename GetEltype<std::remove_reference_t<T>>::value_type;
template <class T, class C>
concept ElementOf = std::convertible_to<T, eltype_t<C>>;

} // namespace utils

template <typename A, typename B> struct PromoteEltype {

  using elta = utils::eltype_t<A>;
  using eltb = utils::eltype_t<B>;
  using value_type =
    std::conditional_t<std::convertible_to<A, eltb>, eltb,
                       std::conditional_t<std::convertible_to<B, elta>, elta,
                                          std::common_type_t<elta, eltb>>>;
};

export namespace utils {

template <typename A, typename B>
using promote_eltype_t =
  std::remove_cvref_t<typename PromoteEltype<A, B>::value_type>;

template <typename T>
concept StaticInt =
  std::is_same_v<T, std::integral_constant<typename T::value_type, T::value>>;

// static_assert(utils::StaticInt<std::integral_constant<ptrdiff_t, 3>>);
// static_assert(!utils::StaticInt<int64_t>);

} // namespace utils

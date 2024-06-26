module;

#ifdef NDEBUG
#ifdef __cpp_lib_unreachable
#include <utility>
#include <version>
#endif
#endif

export module invariant;

export namespace utils {

[[gnu::artificial, gnu::always_inline]] constexpr inline void
invariant(bool condition) {
#ifndef NDEBUG
  if (!condition) [[unlikely]]
    __builtin_trap();
#endif
}
template <typename T>
[[gnu::artificial, gnu::always_inline]] constexpr inline void invariant(T x,
                                                                        T y) {
#ifndef NDEBUG
  if (x != y) [[unlikely]]
    __builtin_trap();
#endif
}

[[gnu::artificial, gnu::always_inline]] constexpr inline void
assume(bool condition) {
#ifndef NDEBUG
  if (!condition) [[unlikely]]
    __builtin_trap();
#else
#if defined(__has_cpp_attribute) && __has_cpp_attribute(assume)
  [[assume(condition)]];
#else
  if (!condition) {
#ifdef __cpp_lib_unreachable
    std::unreachable();
#else
#ifdef __has_builtin
#if __has_builtin(__builtin_unreachable)
    __builtin_unreachable();
#endif
#endif
#endif // __cpp_lib_unreachable
#endif // assume
}
#endif // NDEBUG
}
template <typename T>
[[gnu::artificial, gnu::always_inline]] constexpr inline void assumeeq(T x,
                                                                       T y) {
#ifndef NDEBUG
  if (x != y) [[unlikely]]
    __builtin_trap();
#else
#if defined(__has_cpp_attribute) && __has_cpp_attribute(assume)
  [[assume(condition)]];
#else
  if (x != y) {
#ifdef __cpp_lib_unreachable
    std::unreachable();
#else
#ifdef __has_builtin
#if __has_builtin(__builtin_unreachable)
    __builtin_unreachable();
#endif
#endif
#endif // __cpp_lib_unreachable
#endif // assume
}
#endif // NDEBUG
}

} // namespace utils

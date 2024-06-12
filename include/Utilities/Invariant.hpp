#pragma once

#ifndef NDEBUG

namespace poly::utils {

[[gnu::artificial, gnu::always_inline]] constexpr inline void
invariant(bool condition) {
  if (!condition) [[unlikely]]
    __builtin_trap();
}
template <typename T>
[[gnu::artificial, gnu::always_inline]] constexpr inline void invariant(T x,
                                                                        T y) {
  if (x != y) [[unlikely]]
    __builtin_trap();
}
// we want gdb-friendly builtin trap
#define ASSERT(condition) ::poly::utils::invariant(condition)
} // namespace poly::utils
#else // ifdef NDEBUG
#ifdef __cpp_lib_unreachable
#include <utility>
#include <version>
#endif
namespace poly::utils {
// [[gnu::artificial, gnu::always_inline]] constexpr inline void invariant(bool)
// {}

// template <typename T>
// [[gnu::artificial, gnu::always_inline]] constexpr inline void invariant(T, T)
// {}
[[gnu::artificial, gnu::always_inline]] constexpr inline void
invariant(bool condition) {

#ifdef __has_cpp_attribute
#if __has_cpp_attribute(assume)
  [[assume(condition)]];
#endif
#endif
  if (!condition) {
#ifdef __cpp_lib_unreachable
    std::unreachable();
#else
#ifdef __has_builtin
#if __has_builtin(__builtin_unreachable)
    __builtin_unreachable();
#endif
#endif
#endif
  }
}
template <typename T>
[[gnu::artificial, gnu::always_inline]] constexpr inline void invariant(T x,
                                                                        T y) {
#ifdef __has_cpp_attribute
#if __has_cpp_attribute(assume)
  [[assume(x == y)]];
#endif
#endif
  if (x != y) {
#ifdef __cpp_lib_unreachable
    std::unreachable();
#else
#ifdef __has_builtin
#if __has_builtin(__builtin_unreachable)
    __builtin_unreachable();
#endif
#endif
#endif
  }
}
#define ASSERT(condition) ((void)0)
} // namespace poly::utils
#endif // ifdef NDEBUG

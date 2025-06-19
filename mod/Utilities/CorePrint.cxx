#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#ifndef USE_MODULE
#include <array>
#include <charconv>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <type_traits>
#else
export module CorePrint;
import std;
#endif

#ifdef USE_MODULE
export namespace utils {
#else
namespace utils {
#endif

constexpr auto getChars(std::integral auto x, bool ptr = false)
  -> std::array<char, 21> {
  std::array<char, 21> ret{};
  auto *b = ret.data();
  if (ptr) {
    ret[0] = '0';
    ret[1] = 'x';
    b += 2;
  }
  if (std::to_chars(b, ret.data() + ret.size() - 1, x, 16).ec ==
      std::errc::value_too_large)
    __builtin_trap();
  return ret;
}
constexpr auto getChars(double x) -> std::array<char, 21> {
  std::array<char, 21> ret{};
  if (std::to_chars(ret.data(), ret.data() + ret.size() - 1, x).ec !=
      std::errc::value_too_large)
    return ret;
  ret[0] = 'I';
  ret[1] = 'N';
  ret[2] = 'V';
  ret[3] = 'A';
  ret[4] = 'L';
  ret[5] = 'I';
  ret[6] = 'D';
  ret[7] = ' ';
  ret[8] = 'N';
  ret[9] = 'A';
  ret[10] = 'M';
  ret[11] = 'E';
  ret[12] = '\0';
  return ret;
}
constexpr auto getPtrChars(const void *p) -> std::array<char, 21> {
  return getChars(std::bit_cast<std::uint64_t>(p), true);
}

// Concept for objects with a .print() member function
template <typename T>
concept HasPrintMethod = requires(const T &x) {
  { x.print() } -> std::same_as<void>;
};

inline void print(char c) { std::putchar(c); }
inline void print(const char *c) { std::fputs(c, stdout); }
inline void print(const std::string &s) { std::fwrite(s.data(), 1, s.size(), stdout); }
inline void println(const char *c) { std::puts(c); }
inline void println(const std::string &s) { print(s); std::putchar('\n'); }

// Print function for integral and floating point types
template <typename T>
inline void print(const T &x)
requires(std::integral<T> || std::floating_point<T>)
{
  print(getChars(x).data());
}

// Print function for types with .print() member function
template <HasPrintMethod T> inline void print(const T &x) { x.print(); }

// Println variants
template <typename T>
inline void println(const T &x)
requires(std::integral<T> || std::floating_point<T>)
{
  println(getChars(x).data());
}

template <HasPrintMethod T> inline void println(const T &x) {
  x.print();
  print('\n');
}

// Variadic print functions using const references
template <typename... Args>
inline void print(const auto &x, const auto &y, const Args &...z) {
  print(x);
  print(y);
  (print(z), ...);
}

template <typename... Args>
inline void println(const auto &x, const auto &y, const Args &...z) {
  print(x);         // pop off first
  println(y, z...); // recurse
}

} // namespace utils

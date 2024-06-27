module;

#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <type_traits>

export module arrayprint;

template <std::integral T> consteval auto maxPow10() -> size_t {
  if constexpr (sizeof(T) == 1) return 3;
  else if constexpr (sizeof(T) == 2) return 5;
  else if constexpr (sizeof(T) == 4) return 10;
  else if constexpr (std::signed_integral<T>) return 19;
  else return 20;
}

template <std::unsigned_integral T> constexpr auto countDigits(T x) {
  std::array<T, maxPow10<T>() + 1> powers;
  powers[0] = 0;
  powers[1] = 10;
  for (ptrdiff_t i = 2; i < std::ssize(powers); i++)
    powers[i] = powers[i - 1] * 10;
  std::array<T, sizeof(T) * 8 + 1> bits;
  if constexpr (sizeof(T) == 8) {
    bits = {1,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,
            6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  9,  9,  9,  10, 10, 10, 10,
            11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16,
            16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 19, 20};
  } else if constexpr (sizeof(T) == 4) {
    bits = {1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4,  5,  5, 5,
            6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10};
  } else if constexpr (sizeof(T) == 2) {
    bits = {1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5};
  } else if constexpr (sizeof(T) == 1) {
    bits = {1, 1, 1, 1, 2, 2, 2, 3, 3};
  }
  T digits;
  if constexpr (std::same_as<T, char>)
    digits =
      bits[8 * sizeof(unsigned char) - std::countl_zero((unsigned char)x)];
  else digits = bits[8 * sizeof(T) - std::countl_zero(x)];
  return std::make_signed_t<T>(digits - (x < powers[digits - 1]));
}
template <std::signed_integral T> constexpr auto countDigits(T x) -> T {
  using U = std::make_unsigned_t<T>;
  if (x == std::numeric_limits<T>::min()) return T(sizeof(T) == 8 ? 20 : 11);
  return countDigits<U>(U(std::abs(x))) + T{x < 0};
}

template <typename T> inline auto countDigits(T *x) -> int {
  return countDigits(std::bit_cast<uintptr_t>(x));
}

export namespace utils {
template <typename T>
concept Printable = std::same_as<T, double> || requires(std::ostream &os, T x) {
  { os << x } -> std::same_as<std::ostream &>;
  { countDigits(x) };
};

static_assert(Printable<int64_t>);
void print_obj(std::ostream &os, Printable auto x) { os << x; };

inline auto printVector(std::ostream &os, auto B, auto E) -> std::ostream & {
  os << "[ ";
  if (B != E) {
    print_obj(os, *B);
    for (; ++B != E;) print_obj(os << ", ", *B);
  }
  os << " ]";
  return os;
}
} // namespace utils

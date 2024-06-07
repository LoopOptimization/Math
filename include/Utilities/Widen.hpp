#pragma once

#include <concepts>
#include <cstdint>
namespace poly::utils {
constexpr auto widen(std::signed_integral auto x) {
  static constexpr auto sz = sizeof(decltype(x));
  if constexpr (sz == 8) return __int128_t(x);
  else if constexpr (sz == 4) return int64_t(x);
  else if constexpr (sz == 2) return int32_t(x);
  else {
    static_assert(sz == 1);
    return int16_t(x);
  }
}
constexpr auto widen(std::unsigned_integral auto x) {
  static constexpr auto sz = sizeof(decltype(x));
  if constexpr (sz == 8) return __uint128_t(x);
  else if constexpr (sz == 4) return uint64_t(x);
  else if constexpr (sz == 2) return uint32_t(x);
  else {
    static_assert(sz == 1);
    return uint16_t(x);
  }
}
} // namespace poly::utils

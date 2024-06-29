module;

#include <cstddef>
#include <type_traits>

export module SIMD:Vec;

template <ptrdiff_t W, typename T>
using Vec_ [[gnu::vector_size(W * sizeof(T))]] = T;

export namespace simd {
template <ptrdiff_t W, typename T>
using Vec = std::conditional_t<W == 1, T, Vec_<W, T>>;
} // namespace simd
#ifdef __x86_64__
#ifdef __AVX512F__
inline constexpr ptrdiff_t REGISTERS = 32;
inline constexpr ptrdiff_t VECTORWIDTH = 64;
#else // not __AVX512F__
inline constexpr ptrdiff_t REGISTERS = 16;
#ifdef __AVX__
inline constexpr ptrdiff_t VECTORWIDTH = 32;
#else  // no AVX
inline constexpr ptrdiff_t VECTORWIDTH = 16;
#endif // no AVX
#endif
#else  // not __x86_64__
inline constexpr ptrdiff_t REGISTERS = 32;
inline constexpr ptrdiff_t VECTORWIDTH = 16;
#endif // __x86_64__


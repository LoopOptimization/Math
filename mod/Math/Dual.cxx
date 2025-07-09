module;
#include "Macros.hxx"

// Implementation unit for Dual module - explicit instantiations only
module Dual;

namespace math {

// Explicit instantiation definitions for common dual number types
// Using double with common SIMD widths minus 1 for partials

#ifdef __AVX__
// For 256-bit SIMD (4 doubles): 3 partials
template struct Dual<double, 4, false>;
template struct Dual<double, 4, true>;
template auto exp<double, 4>(const Dual<double, 4> &) -> Dual<double, 4>;
template auto log<double, 4>(const Dual<double, 4> &) -> Dual<double, 4>;
template auto log1p<double, 4>(const Dual<double, 4> &) -> Dual<double, 4>;
template auto log2<double, 4>(const Dual<double, 4> &) -> Dual<double, 4>;
template auto sigmoid<double, 4>(const Dual<double, 4> &) -> Dual<double, 4>;
template auto softplus<double, 4>(const Dual<double, 4> &) -> Dual<double, 4>;
#else

// Scalar dual numbers (single variable AD)
template struct Dual<double, 2, false>;
template struct Dual<double, 2, true>;
template auto exp<double, 2>(const Dual<double, 2> &) -> Dual<double, 2>;
template auto log<double, 2>(const Dual<double, 2> &) -> Dual<double, 2>;
template auto log2<double, 2>(const Dual<double, 2> &) -> Dual<double, 2>;
template auto log1p<double, 2>(const Dual<double, 2> &) -> Dual<double, 2>;
template auto sigmoid<double, 2>(const Dual<double, 2> &) -> Dual<double, 2>;
template auto softplus<double, 2>(const Dual<double, 2> &) -> Dual<double, 2>;
#endif

#ifdef __AVX512F__
// For 512-bit SIMD (8 doubles): 7 partials
template struct Dual<double, 7, false>;
template struct Dual<double, 7, true>;
template auto exp<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto log2<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto log1p<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto log<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto softplus<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto sigmoid<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;

// template struct Dual<double, 8, false>;
// template auto exp<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;
// template auto log2<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;
// template auto log1p<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;
// template auto log<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;
// template auto softplus<double, 8>(const Dual<double, 8> &) -> Dual<double,
// 8>; template auto sigmoid<double, 8>(const Dual<double, 8> &) -> Dual<double,
// 8>;
#endif

} // namespace math

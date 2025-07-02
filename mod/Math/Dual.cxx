module;
#include "Macros.hxx"

// Implementation unit for Dual module - explicit instantiations only
module Dual;

namespace math {

// Explicit instantiation definitions for common dual number types
// Using double with common SIMD widths minus 1 for partials

// For 256-bit SIMD (4 doubles): 3 partials
template struct Dual<double, 3, false>;
template struct Dual<double, 3, true>;

// For 512-bit SIMD (8 doubles): 7 partials
template struct Dual<double, 7, false>;
template struct Dual<double, 7, true>;

// Common general purpose size: 8 partials (as used in gradient function)
template struct Dual<double, 8, false>;

// Scalar dual numbers (single variable AD)
template struct Dual<double, 1, false>;
template struct Dual<double, 1, true>;

// Note: Operator instantiations removed because operators are commented out in
// interface

// Explicit instantiation of elementary functions for common types
template auto exp<double, 1>(const Dual<double, 1> &) -> Dual<double, 1>;
template auto exp<double, 3>(const Dual<double, 3> &) -> Dual<double, 3>;
template auto exp<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto exp<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;

template auto log<double, 1>(const Dual<double, 1> &) -> Dual<double, 1>;
template auto log<double, 3>(const Dual<double, 3> &) -> Dual<double, 3>;
template auto log<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto log<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;

template auto log2<double, 1>(const Dual<double, 1> &) -> Dual<double, 1>;
template auto log2<double, 3>(const Dual<double, 3> &) -> Dual<double, 3>;
template auto log2<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto log2<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;

template auto log1p<double, 1>(const Dual<double, 1> &) -> Dual<double, 1>;
template auto log1p<double, 3>(const Dual<double, 3> &) -> Dual<double, 3>;
template auto log1p<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto log1p<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;

template auto sigmoid<double, 1>(const Dual<double, 1> &) -> Dual<double, 1>;
template auto sigmoid<double, 3>(const Dual<double, 3> &) -> Dual<double, 3>;
template auto sigmoid<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto sigmoid<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;

template auto softplus<double, 1>(const Dual<double, 1> &) -> Dual<double, 1>;
template auto softplus<double, 3>(const Dual<double, 3> &) -> Dual<double, 3>;
template auto softplus<double, 7>(const Dual<double, 7> &) -> Dual<double, 7>;
template auto softplus<double, 8>(const Dual<double, 8> &) -> Dual<double, 8>;

// Note: Utility function instantiations removed due to templated auto
// parameters

} // namespace math

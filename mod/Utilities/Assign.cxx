#pragma once
#include <Math/Matrix.hpp>
#include <Utilities/TypePromotion.hpp>
#include <concepts>
#include <functional>
namespace poly::utils {

struct CopyAssign {};
struct NoRowIndex {};

template <typename D, typename S, typename Op>
[[gnu::artificial, gnu::always_inline]] inline constexpr void
assign(D &&d, const S &s, Op op) {
  if constexpr (std::same_as<Op, CopyAssign>) d = s;
  else if constexpr (std::same_as<Op, std::plus<>>) d += s;
  else if constexpr (std::same_as<Op, std::minus<>>) d -= s;
  else if constexpr (std::same_as<Op, std::multiplies<>>) d *= s;
  else if constexpr (std::same_as<Op, std::divides<>>) d /= s;
  else d = op(d, s);
}

template <typename D, typename S, typename R, typename C, typename Op>
[[gnu::artificial, gnu::always_inline]] inline constexpr void
assign(D &d, const S &s, R r, C c, Op op) {
  constexpr bool no_row_ind = std::same_as<R, NoRowIndex>;
  if constexpr (std::convertible_to<S, utils::eltype_t<D>>)
    if constexpr (no_row_ind) assign(d[c], s, op);
    else assign(d[r, c], s, op);
  else if constexpr (math::RowVector<S>)
    if constexpr (no_row_ind) assign(d[c], s[c], op);
    else assign(d[r, c], s[c], op);
  else if constexpr (math::ColVector<S>)
    if constexpr (no_row_ind) assign(d[c], s[c], op);
    else assign(d[r, c], s[r], op);
  else if constexpr (std::same_as<Op, CopyAssign>)
    if constexpr (no_row_ind) d[c] = s[c];
    else d[r, c] = s[r, c];
  else if constexpr (std::same_as<Op, std::plus<>>)
    if constexpr (no_row_ind) d[c] += s[c];
    else d[r, c] += s[r, c];
  else if constexpr (std::same_as<Op, std::minus<>>)
    if constexpr (no_row_ind) d[c] -= s[c];
    else d[r, c] -= s[r, c];
  else if constexpr (std::same_as<Op, std::multiplies<>>)
    if constexpr (no_row_ind) d[c] *= s[c];
    else d[r, c] *= s[r, c];
  else if constexpr (std::same_as<Op, std::divides<>>)
    if constexpr (no_row_ind) d[c] /= s[c];
    else d[r, c] /= s[r, c];
  else if constexpr (no_row_ind) d[c] = op(const_cast<const D &>(d)[c], s[c]);
  else d[r, c] = op(const_cast<const D &>(d)[r, c], s[r, c]);
}

} // namespace poly::utils


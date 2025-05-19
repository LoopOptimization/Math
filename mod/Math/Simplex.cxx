#ifdef USE_MODULE
module;
#else
#pragma once
#endif
#include <cassert>

#ifndef USE_MODULE
#include "Alloc/Arena.cxx"
#include "Alloc/Mallocator.cxx"
#include "Containers/BitSets.cxx"
#include "Math/Array.cxx"
#include "Math/ArrayConcepts.cxx"
#include "Math/AxisTypes.cxx"
#include "Math/Comparisons.cxx"
#include "Math/Constraints.cxx"
#include "Math/ExpressionTemplates.cxx"
#include "Math/GreatestCommonDivisor.cxx"
#include "Math/Indexing.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/MatrixDimensions.cxx"
#include "Math/NormalForm.cxx"
#include "Math/Ranges.cxx"
#include "Math/Rational.cxx"
#include "SIMD/Intrin.cxx"
#include "SIMD/Vec.cxx"
#include "Utilities/Invariant.cxx"
#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <ostream>
#else
export module Simplex;

import Allocator;
import Arena;
import Array;
import ArrayConcepts;
import AxisTypes;
import BitSet;
import Comparisons;
import Constraints;
import ExprTemplates;
import GCD;
import Invariant;
import ManagedArray;
import MatDim;
import NormalForm;
import Range;
import Rational;
import SIMD;
import std;
#endif

#ifndef NDEBUG
#define DEBUGUSED [[gnu::used]]
#else
#define DEBUGUSED
#endif

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif

// #define VERBOSESIMPLEX

/// Tableau for the Simplex algorithm.
/// We need a core Simplex type that is unmanaged
/// then for convenience, it would be nice to manage it.
/// Ideally, we could have a type hierarchy of
/// unmanaged -> managed
/// with some API to make the managed generic.
/// We also want the managed to be automatically demotable to unmanaged,
/// to avoid unnecessary specialization.
///
/// Slack variables are sorted first.
class Simplex {
  using index_t = int;
  using value_type = std::int64_t;

  template <std::integral T>
  static constexpr auto alignOffset(std::ptrdiff_t x) -> std::ptrdiff_t {
    --x;
    std::ptrdiff_t W = simd::VECTORWIDTH / sizeof(T); // simd::Width<T>;
    std::ptrdiff_t nW = -W;
    x += W;
    x &= nW;
    return x;
    // return (--x + simd::Width<T>)&(-simd::Width<T>);
  }
  [[gnu::returns_nonnull, nodiscard]] constexpr auto basicConsPointer() const
    -> index_t * {
    void *p = const_cast<char *>(memory_);
    return std::assume_aligned<simd::VECTORWIDTH>(static_cast<index_t *>(p));
  }
  [[gnu::returns_nonnull, nodiscard]] constexpr auto basicVarsPointer() const
    -> index_t * {
    std::ptrdiff_t offset = alignOffset<index_t>(reservedBasicConstraints());
    return std::assume_aligned<simd::VECTORWIDTH>(basicConsPointer() + offset);
  }
  // offset in bytes
  static constexpr auto tableauOffset(std::ptrdiff_t cons, std::ptrdiff_t vars)
    -> std::ptrdiff_t {
    std::ptrdiff_t coff = alignOffset<index_t>(cons);
    std::ptrdiff_t voff = alignOffset<index_t>(vars);
    std::ptrdiff_t offset = coff + voff;
    // std::ptrdiff_t offset =
    //   alignOffset<index_t>(cons) + alignOffset<index_t>(vars);
    return static_cast<std::ptrdiff_t>(sizeof(index_t)) * offset;
  }
  [[gnu::returns_nonnull, nodiscard]] constexpr auto tableauPointer() const
    -> value_type * {
    std::ptrdiff_t offset =
      tableauOffset(reservedBasicConstraints(), reservedBasicVariables());
    void *p = const_cast<char *>(memory_) + offset;
    return std::assume_aligned<simd::VECTORWIDTH>(static_cast<value_type *>(p));
  }
  Row<> num_constraints_{};
  Col<> num_vars_{};
  Capacity<> constraint_capacity_;
  RowStride<> var_capacity_p1_; // varCapacity + 1
#ifndef NDEBUG
  bool in_canonical_form_{false};
#endif
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#else
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) // FAM
  alignas(simd::VECTORWIDTH) char memory_[];
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#else
#pragma clang diagnostic pop
#endif

  // (varCapacity+1)%simd::Width<std::int64_t>==0
  static constexpr auto alignVarCapacity(RowStride<> rs) -> RowStride<> {
    static constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
    return stride((std::ptrdiff_t(rs) + W) & -W);
  }
  [[nodiscard]] static constexpr auto reservedTableau(std::ptrdiff_t cons,
                                                      std::ptrdiff_t vars)
    -> std::ptrdiff_t {
    return static_cast<std::ptrdiff_t>(sizeof(value_type)) *
           ((cons + 1) * vars);
  }
  static constexpr auto requiredMemory(std::ptrdiff_t cons, std::ptrdiff_t vars)
    -> std::size_t {
    std::ptrdiff_t base = static_cast<std::ptrdiff_t>(sizeof(Simplex)),
                   indices = tableauOffset(cons, vars),
                   tableau = reservedTableau(cons, vars);
    return static_cast<std::size_t>(base + indices + tableau);
  }

public:
  // tableau is constraint * var matrix w/ extra col for LHS
  // and extra row for objective function
  [[nodiscard]] constexpr auto reservedBasicConstraints() const
    -> std::ptrdiff_t {
    return std::ptrdiff_t(var_capacity_p1_) - 1;
  }
  [[nodiscard]] constexpr auto reservedBasicVariables() const
    -> std::ptrdiff_t {
    return std::ptrdiff_t(constraint_capacity_);
  }

  // [[nodiscard]] constexpr auto intsNeeded() const -> std::ptrdiff_t {
  //   return reservedTableau() + reservedBasicConstraints() +
  //          reservedBasicVariables();
  // }
  // We with to align every row of `getConstraints()` and `getTableau()`.
  // To do this, we
  // 0. (varCapacity+1) % simd::Width<std::int64_t> == 0
  // 1. offset...
  //
  // tableau has a stride of `varCapacityP1`, which is the maximum var capacity
  // +1. The `+1` is room for the LHS.
  //
  /// [ value | objective function ]
  /// [ LHS   | tableau            ]
  [[nodiscard]] constexpr auto getTableau() const -> PtrMatrix<value_type> {
    //
    return {tableauPointer(), StridedDims<>{
                                ++auto(num_constraints_),
                                ++auto(num_vars_),
                                var_capacity_p1_,
                              }};
  }
  // NOLINTNEXTLINE(readability-make-member-function-const)
  [[nodiscard]] constexpr auto getTableau() -> MutPtrMatrix<value_type> {
    return {tableauPointer(), StridedDims<>{
                                ++auto(num_constraints_),
                                ++auto(num_vars_),
                                var_capacity_p1_,
                              }};
  }
  [[nodiscard]] constexpr auto getConstraints() const -> PtrMatrix<value_type> {
    return {tableauPointer() + std::ptrdiff_t(var_capacity_p1_),
            StridedDims<>{
              num_constraints_,
              ++auto(num_vars_),
              var_capacity_p1_,
            }};
  }
  constexpr void zeroConstraints() {
    MutPtrVector<value_type>{
      tableauPointer() + std::ptrdiff_t(var_capacity_p1_),
      length(std::ptrdiff_t(num_constraints_) *
             std::ptrdiff_t(var_capacity_p1_)),
    }
      .zero();
  }
  // NOLINTNEXTLINE(readability-make-member-function-const)
  [[nodiscard]] constexpr auto getConstraints() -> MutPtrMatrix<value_type> {
    return {tableauPointer() + std::ptrdiff_t(var_capacity_p1_),
            StridedDims<>{
              num_constraints_,
              ++auto(num_vars_),
              var_capacity_p1_,
            }};
  }
  [[nodiscard]] constexpr auto getBasicConstraints() const
    -> PtrVector<index_t> {
    return {basicConsPointer(), aslength(num_vars_)};
  }
  // maps variables to constraints
  // if the constraint is < 0, then that means there is no associated basic
  // constraint, and the variable's value is `0`, i.e. it is non-basic
  // Otherwise, it is the index of the only non-zero constraint for that
  // variable.
  [[nodiscard]] constexpr auto getBasicConstraints() -> MutPtrVector<index_t> {
    return {basicConsPointer(), aslength(num_vars_)};
  }
  // maps constraints to variables
  [[nodiscard]] constexpr auto getBasicVariables() const -> PtrVector<index_t> {
    return {basicVarsPointer(), length(std::ptrdiff_t(num_constraints_))};
  }
  [[nodiscard]] constexpr auto getBasicVariables() -> MutPtrVector<index_t> {
    return {basicVarsPointer(), length(std::ptrdiff_t(num_constraints_))};
  }
  [[nodiscard]] constexpr auto getCost() const -> PtrVector<value_type> {
    return {tableauPointer(), length(std::ptrdiff_t(num_vars_) + 1z)};
  }
  // NOLINTNEXTLINE(readability-make-member-function-const)
  [[nodiscard]] constexpr auto getCost() -> MutPtrVector<value_type> {
    return {tableauPointer(), length(std::ptrdiff_t(num_vars_) + 1z)};
  }
  // maps variables to constraints; <0 if variable is not basic
  [[nodiscard]] constexpr auto getBasicConstraint(std::ptrdiff_t var) const
    -> index_t {
    return getBasicConstraints()[var];
  }
  // get the variable that is basic associated with this constraint;
  // must always be valid while the simplex is in canonical form.
  [[nodiscard]] constexpr auto getBasicVariable(std::ptrdiff_t constraint) const
    -> index_t {
    return getBasicVariables()[constraint];
  }
  [[nodiscard]] constexpr auto getObjectiveCoefficient(std::ptrdiff_t i) const
    -> value_type {
    return getCost()[++i];
  }
  [[nodiscard]] constexpr auto getObjectiveValue() -> value_type & {
    return getCost()[0];
  }
  [[nodiscard]] constexpr auto getObjectiveValue() const -> value_type {
    return getCost()[0];
  }
  constexpr void simplifySystem() {
#ifndef NDEBUG
    in_canonical_form_ = false;
#endif
    auto C{getConstraints()};
#ifndef NDEBUG
    for (std::ptrdiff_t r = 0, R = std::ptrdiff_t(numRows(C)); r < R; ++r)
      for (std::ptrdiff_t c = 0, N = std::ptrdiff_t(numCols(C)); c < N; ++c)
        invariant(C[r, c] != std::numeric_limits<std::int64_t>::min());
#endif
    NormalForm::solveSystemSkip(C);
    truncateConstraints(std::ptrdiff_t(NormalForm::numNonZeroRows(C)));
  }
#ifndef NDEBUG
  constexpr void assertCanonical() const {
    PtrMatrix<value_type> C{getTableau()};
    PtrVector<index_t> basic_vars{getBasicVariables()};
    PtrVector<index_t> basic_cons{getBasicConstraints()};
    for (std::ptrdiff_t v = 0; v < basic_cons.size();) {
      index_t c = basic_cons[v++];
      if (c < 0) continue;
      if (!allZero(C[_(1, 1 + c), v])) __builtin_trap();
      if (!allZero(C[_(2 + c, end), v])) __builtin_trap();
      if (std::ptrdiff_t(basic_vars[c]) != (v - 1)) __builtin_trap();
    }
    for (std::ptrdiff_t c = 1; c < C.numRow(); ++c) {
      index_t v = basic_vars[c - 1];
      if (std::ptrdiff_t(v) < basic_cons.size()) {
        invariant(c - 1, std::ptrdiff_t(basic_cons[v]));
        invariant(C[c, v + 1] >= 0);
      }
      invariant(C[c, 0] >= 0);
    }
  }
#endif
  [[nodiscard]] constexpr auto getConstants()
    -> MutStridedVector<std::int64_t> {
    return getTableau()[_(1, end), 0];
  }
  [[nodiscard]] constexpr auto getConstants() const
    -> StridedVector<std::int64_t> {
    return getTableau()[_(1, end), 0];
  }
  constexpr void truncateConstraints(std::ptrdiff_t i) {
    invariant(std::ptrdiff_t(num_constraints_) <=
              std::ptrdiff_t(constraint_capacity_));
    invariant(i >= 0z);
    invariant(i <= num_constraints_);
    num_constraints_ = row(i);
  }
  constexpr void setNumCons(std::ptrdiff_t i) {
    invariant(i <= constraint_capacity_);
    num_constraints_ = row(i);
  }
  constexpr void setNumVars(std::ptrdiff_t i) {
    invariant(i < var_capacity_p1_);
    num_vars_ = col(i);
  }
  constexpr void truncateVars(std::ptrdiff_t i) {
    invariant(i <= num_vars_);
    num_vars_ = col(i);
  }
  [[nodiscard]] constexpr auto getNumCons() const -> std::ptrdiff_t {
    invariant(num_constraints_ >= 0);
    return std::ptrdiff_t(num_constraints_);
  }
  [[nodiscard]] constexpr auto getNumVars() const -> std::ptrdiff_t {
    invariant(num_vars_ >= 0);
    return std::ptrdiff_t(num_vars_);
  }
  [[nodiscard]] constexpr auto getConCap() const -> Capacity<> {
    invariant(constraint_capacity_ >= 0);
    return constraint_capacity_;
  }
  [[nodiscard]] constexpr auto getVarCap() const -> RowStride<> {
    invariant(var_capacity_p1_ > 0);
    return --auto{var_capacity_p1_};
  }
  constexpr void deleteConstraint(std::ptrdiff_t c) {
    auto basic_cons = getBasicConstraints();
    auto basic_vars = getBasicVariables();
    auto constraints = getConstraints();
    --num_constraints_;
    if (auto basic_var = basic_vars[c]; basic_var >= 0)
      basic_cons[basic_var] = -1;
    if (c == num_constraints_) return;
    auto basic_var = basic_vars[std::ptrdiff_t(num_constraints_)];
    basic_vars[c] = basic_var;
    if (basic_var >= 0) basic_cons[basic_var] = index_t(c);
    constraints[c, _] << constraints[num_constraints_, _];
  }

  // AbstractVector
  struct Solution : Expr<Rational, Solution> {
    using value_type = Rational;
    // view of tableau dropping const column
    Valid<const Simplex> simplex_;
    Length<> skipped_vars_;
    Length<> num_vars_;
    class iterator { // NOLINT(readability-identifier-naming)
      const Solution *sol_;
      std::ptrdiff_t i_;

    public:
      using value_type = Rational;
      constexpr iterator(const Solution *s, std::ptrdiff_t j)
        : sol_(s), i_(j) {}
      constexpr iterator() = default;
      constexpr iterator(const iterator &) = default;
      constexpr auto operator=(const iterator &) -> iterator & = default;
      auto operator*() const -> Rational { return (*sol_)[i_]; }
      constexpr auto operator++() -> iterator & {
        ++i_;
        return *this;
      }
      constexpr auto operator++(int) -> iterator {
        auto tmp = *this;
        ++i_;
        return tmp;
      }
      constexpr auto operator--() -> iterator & {
        --i_;
        return *this;
      }
      constexpr auto operator--(int) -> iterator {
        auto tmp = *this;
        --i_;
        return tmp;
      }
      friend constexpr auto operator==(iterator a, iterator b) -> bool {
        return a.i_ == b.i_;
      }
      friend constexpr auto operator!=(iterator a, iterator b) -> bool {
        return a.i_ != b.i_;
      }
      constexpr auto operator-(iterator b) const -> std::ptrdiff_t {
        return std::ptrdiff_t(i_) - b.i_;
      }
      constexpr auto operator+(std::ptrdiff_t n) const -> iterator {
        return {sol_, i_ + n};
      }
    };
    [[nodiscard]] constexpr auto begin() const -> iterator { return {this, 0}; }
    [[nodiscard]] constexpr auto end() const -> iterator {
      return {this, std::ptrdiff_t(num_vars_ - skipped_vars_)};
    }

    [[nodiscard]] constexpr auto operator[](std::ptrdiff_t i) const
      -> Rational {
      invariant(i >= 0);
      return simplex_->getVarValue(i + std::ptrdiff_t(skipped_vars_));
    }
    [[nodiscard]] constexpr auto operator[](OffsetEnd k) const -> Rational {
      return simplex_->getVarValue(std::ptrdiff_t(num_vars_) - k.offset_);
    }
    [[nodiscard]] constexpr auto operator[](ScalarRelativeIndex auto i) const
      -> Rational {
      return (*this)[calcOffset(size(), i)];
    }
    template <typename B, typename E>
    constexpr auto operator[](Range<B, E> r) const -> Solution {
      return (*this)[canonicalizeRange(r, size())];
    }
    constexpr auto operator[](Range<std::ptrdiff_t, std::ptrdiff_t> r) const
      -> Solution {
      return {.simplex_ = simplex_,
              .skipped_vars_ = length(std::ptrdiff_t(skipped_vars_) + r.b_),
              .num_vars_ = length(std::ptrdiff_t(skipped_vars_) + r.e_)};
    }
    [[nodiscard]] constexpr auto size() const -> std::ptrdiff_t {
      return std::ptrdiff_t(num_vars_ - skipped_vars_);
    }
    [[nodiscard]] constexpr auto view() const -> Solution { return *this; };

    [[nodiscard]] constexpr auto denomLCM() const -> std::int64_t {
      std::int64_t l = 1;
      for (auto r : *this) l = lcm(l, r.denominator_);
      return l;
    }
#ifndef NDEBUG
    [[gnu::used]] void dump() const { std::cout << *this << '\n'; }
#endif
  private:
    friend auto operator<<(std::ostream &os, Solution sol) -> std::ostream & {
      os << "Simplex::Solution[";
      bool print_comma = false;
      for (Rational b : sol) {
        if (print_comma) os << ", ";
        print_comma = true;
        os << b;
      }
      os << "]";
      return os;
    }
  };
  DEBUGUSED [[nodiscard]] constexpr auto getSolution() const -> Solution {
    return {
      .simplex_ = this, .skipped_vars_ = {}, .num_vars_ = aslength(num_vars_)};
  }

  /// simplex.initiateFeasible() -> bool
  /// returns `true` if infeasible, `false ` if feasible
  /// The approach is to first put the equalities into HNF
  /// then, all diagonal elements are basic variables.
  /// For each non-diagonal element, we need to add an augment variable
  /// Then we try to set all augment variables to 0.
  /// If we fail, it is infeasible.
  /// If we succeed, then the problem is feasible, and we're in
  /// canonical form.
  [[nodiscard(
    "returns `true` if infeasible; should check when calling.")]] constexpr auto
  initiateFeasible() -> bool {
    // remove trivially redundant constraints
    simplifySystem();
    // [ I;  X ; b ]
    //
    // original number of variables
    const std::ptrdiff_t num_var = getNumVars();
    MutPtrMatrix<value_type> C{getConstraints()};
    MutPtrVector<index_t> basic_cons{getBasicConstraints()};
    // initialize to `-2`
    basic_cons << -2;
    // first pass, we make sure the equalities are >= 0
    // and we eagerly try and find columns with
    // only a single non-0 element.
    for (std::ptrdiff_t c = 0; c < C.numRow(); ++c) {
      std::int64_t &Ceq = C[c, 0];
      std::int64_t sign = 2 * (Ceq >= 0) - 1;
      Ceq *= sign;
      // was initialized to -2
      // if unset here (i.e. == -2) and >1, try to make basic
      // if set earlier and we're resetting (or < 0), set to -1
      for (std::ptrdiff_t v = 0; v < num_var; ++v)
        if (std::int64_t Ccv = C[c, v + 1] *= sign)
          basic_cons[v] =
            (((basic_cons[v] == -2) && (Ccv > 0))) ? index_t(c) : -1;
    }
    // basicCons should now contain either `-1` or an integer >= 0
    // indicating which row contains the only non-zero element; we'll
    // now fill basicVars.
    //
    auto basic_vars{getBasicVariables()};
    basic_vars << -1;
    for (std::ptrdiff_t v = 0; v < num_var; ++v) {
      if (std::int64_t r = basic_cons[v]; r >= 0) {
        if (basic_vars[r] == -1) basic_vars[r] = index_t(v);
        else basic_cons[v] = -1;
      }
    }
#ifndef NDEBUG
    in_canonical_form_ = true;
#endif
    containers::BitSet<math::Vector<std::uint64_t, 8>> aug_vars{};
    // upper bound number of augmentVars is constraintCapacity
    // we push augment vars
    for (std::ptrdiff_t i = 0; i < basic_vars.size(); ++i)
      if (basic_vars[i] == -1) aug_vars.uncheckedInsert(i);
    return (!aug_vars.empty() && removeAugmentVars(aug_vars));
  }
  constexpr auto removeAugmentVars(
    const containers::BitSet<math::Vector<std::uint64_t, 8>> &augmentVars)
    -> bool {
    // TODO: try to avoid reallocating, via reserving enough ahead of time
    std::ptrdiff_t num_augment = augmentVars.size(),
                   old_num_var = std::ptrdiff_t(num_vars_);
    invariant(num_augment + std::ptrdiff_t(num_vars_) < var_capacity_p1_);
    num_vars_ = col(std::ptrdiff_t(num_vars_) + num_augment);
    MutPtrMatrix<value_type> C{getConstraints()};
    MutPtrVector<index_t> basic_vars{getBasicVariables()};
    MutPtrVector<index_t> basic_cons{getBasicConstraints()};
    MutPtrVector<value_type> costs{getCost()};
    costs << 0;
    C[_, _(old_num_var + 1, end)] << 0;
    {
      std::ptrdiff_t i = 0;
      for (std::ptrdiff_t a : augmentVars) {
        basic_vars[a] = index_t(i) + index_t(old_num_var);
        basic_cons[i + old_num_var] = index_t(a);
        C[a, old_num_var + (++i)] = 1;
        // we now zero out the implicit cost of `1`
        costs[_(begin, old_num_var + 1)] -= C[a, _(begin, old_num_var + 1)];
      }
    }
    assert(
      std::ranges::all_of(basic_vars, [](std::int64_t i) { return i >= 0; }));
    // false/0 means feasible
    // true/non-zero infeasible
    if (runCore()) return true;
    // check for any basic vars set to augment vars, and set them to some
    // other variable (column) instead.
    for (std::ptrdiff_t c = 0; c < C.numRow(); ++c) {
      if (std::ptrdiff_t(basic_vars[c]) >= old_num_var) {
        invariant(C[c, 0] == 0);
        invariant(c == basic_cons[basic_vars[c]]);
        invariant(C[c, basic_vars[c] + 1] >= 0);
        // find var to make basic in its place
        for (std::ptrdiff_t v = old_num_var; v != 0;) {
          // search for a non-basic variable
          // (basicConstraints<0)
          std::int64_t Ccv = C[c, v--];
          if (Ccv == 0 || (basic_cons[v] >= 0)) continue;
          if (Ccv < 0) C[c, _] *= -1;
          for (std::ptrdiff_t i = 0; i < C.numRow(); ++i)
            if (i != c) NormalForm::zeroWithRowOp(C, row(i), row(c), ++col(v));
          basic_vars[c] = index_t(v);
          basic_cons[v] = index_t(c);
          break;
        }
      }
    }
    // all augment vars are now 0
    num_vars_ = col(old_num_var);
#ifndef NDEBUG
    assertCanonical();
#endif
    return false;
  }

  // 1 based to match getBasicConstraints
  [[nodiscard]] static constexpr auto
  getEnteringVariable(PtrVector<std::int64_t> costs) -> std::optional<index_t> {
    // Bland's algorithm; guaranteed to terminate
    auto f = costs.begin(), l = costs.end();
    const auto *neg =
      std::find_if(f, l, [](std::int64_t c) -> bool { return c < 0; });
    return neg != l ? std::distance(f, neg) : std::optional<index_t>{};
  }
  struct NoFilter {
    static constexpr auto operator()(std::int64_t) -> bool { return false; }
  };
  // Searches rows (constraints) for the one to make non-zero for
  // this particular variable. It will be used to zero all other rows.
  // This makes the `entering_variable` basic.
  // The `leaving_variable`, i.e. the variable that will no longer
  // be basic, is `getBasicVariable(leaving_variable)`.
  //
  //
  // Tries to find the row with the lowest numerator/denominator
  // ratio, and returns that row. If any numerator == 0, it'll
  // stop early and return that row.
  [[nodiscard]] static constexpr auto
  getLeavingVariable(PtrMatrix<std::int64_t> C,
                     std::ptrdiff_t entering_variable, const auto &filter)
    -> std::optional<index_t> {
    // inits guarantee first valid is selected
    std::int64_t dj = -1, nj = 0;
    index_t j = 0;
    StridedVector<std::int64_t> cv{C[_, entering_variable + 1]};
    for (index_t i = 1; i < C.numRow(); ++i) {
      if (filter(i)) continue;
      std::int64_t di = cv[i];
      if (di <= 0) continue;
      std::int64_t ni = C[i, 0];
      if (!ni) return --i;
      invariant(ni > 0);
      if ((dj * ni) >= (di * nj)) continue;
      // we could consider something like:
      // so that in case of ties, we prefer having the maximum index
      // be the leaving variable.
      // That didn't really help the existing benchmarks.
      // auto basicVars = getBasicVariables(); // passed in as arg to static fun
      // if ((n * Cio) > (Civ * d)) continue;
      // if ((n * Cio) == (Civ * d) && (basicVars[i - 1] < basicVars[j - 1]))
      //   continue;
      dj = di;
      nj = ni;
      j = i;
    }
    // NOTE: if we fail to find a leaving variable, then `j = 0`,
    // and it will unsigned wrap to `std::ptrdiff_t(-1)`, which indicates
    // an empty `Optional<unsigned int>`
    return j ? --j : std::optional<index_t>{};
  }
  [[nodiscard]] static constexpr auto
  getLeavingVariable(PtrMatrix<std::int64_t> C,
                     std::ptrdiff_t entering_variable)
    -> std::optional<index_t> {
    return getLeavingVariable(C, entering_variable, NoFilter{});
  }

  constexpr auto makeBasic(MutPtrMatrix<std::int64_t> C, index_t enteringVar,
                           std::int64_t f) -> std::int64_t {
    std::optional<index_t> leave_opt = getLeavingVariable(C, enteringVar);
    if (!leave_opt) return 0; // unbounded
    index_t leaving_var = *leave_opt;
    for (std::ptrdiff_t i = 0; i < C.numRow(); ++i) {
      if (i == leaving_var + 1) continue;
      std::int64_t m = NormalForm::zeroWithRowOp(C, row(i), ++row(leaving_var),
                                                 ++col(enteringVar), i ? 0 : f);
      if (!i) f = m;
    }
    // update basic vars and constraints
    MutPtrVector<index_t> basic_vars{getBasicVariables()},
      basic_constraints{getBasicConstraints()};
    index_t old_basic_var = basic_vars[leaving_var];
    basic_vars[leaving_var] = enteringVar;
    basic_constraints[old_basic_var] = -1;
    basic_constraints[enteringVar] = leaving_var;
    return f;
  }
  // returns `true` if it failed to make basic
  constexpr auto makeBasic(MutPtrMatrix<std::int64_t> C, index_t ev) -> bool {
    std::optional<index_t> leave_opt = getLeavingVariable(C, ev);
    if (leave_opt) makeBasic(C, ev, *leave_opt);
    return !leave_opt;
  }
  constexpr void makeBasic(MutPtrMatrix<std::int64_t> C, index_t ev,
                           index_t l_var) {
    index_t leaving_variable = l_var++;
    for (std::ptrdiff_t i = 0; i < C.numRow(); ++i)
      if (i != l_var)
        NormalForm::zeroWithRowOp(C, row(i), row(l_var), ++col(ev));
    // update basic vars and constraints
    MutPtrVector<index_t> basic_vars{getBasicVariables()},
      basic_constraints{getBasicConstraints()};
    index_t old_basic_var = basic_vars[leaving_variable];
    basic_vars[leaving_variable] = ev;
    if (std::ptrdiff_t(old_basic_var) < basic_constraints.size())
      basic_constraints[old_basic_var] = -1;
    basic_constraints[ev] = leaving_variable;
  }
  // run the simplex algorithm, assuming basicVar's costs have been set to
  // 0
  constexpr auto runCore(std::int64_t f = 1) -> Rational {
#ifndef NDEBUG
    assert(in_canonical_form_);
#endif
    //     return runCore(getCostsAndConstraints(), f);
    // }
    // Rational runCore(MutPtrMatrix<std::int64_t> C, std::int64_t f = 1) {
    MutPtrMatrix<std::int64_t> C{getTableau()};
    do {
      // entering variable is the column
      std::optional<index_t> entering_variable =
        getEnteringVariable(C[0, _(1, end)]);
      if (!entering_variable) return Rational::create(C[0, 0], f);
      f = makeBasic(C, *entering_variable, f);
    } while (f);
    return std::numeric_limits<std::int64_t>::max(); // unbounded
  }
  // set basicVar's costs to 0, and then runCore()
  constexpr auto run() -> Rational {
#ifndef NDEBUG
    assert(in_canonical_form_);
    assertCanonical();
#endif
    MutPtrVector<index_t> basic_vars{getBasicVariables()};
    MutPtrMatrix<value_type> C{getTableau()};
    std::int64_t f = 1;
    // zero cost of basic variables to put in canonical form
    for (std::ptrdiff_t c = 0; c < basic_vars.size();) {
      std::int64_t v = basic_vars[c++];
      if ((std::ptrdiff_t(++v) < C.numCol()) && C[0, v])
        f = NormalForm::zeroWithRowOp(C, row(0), row(c), col(v), f);
    }
    return runCore(f);
  }
  // Returns the basic constraint
  // Will be negative if it failed.
  constexpr auto makeBasic(index_t var) -> index_t {
    MutPtrMatrix<value_type> C{getTableau()};
    if (index_t j = getBasicConstraint(var); j >= 0) return j;
    if (makeBasic(C, var)) return getBasicConstraint(var);
    return -1;
  }

  // don't touch variables lex > v
  constexpr void rLexCore(PtrMatrix<std::int64_t> C,
                          PtrVector<std::int64_t> costs, std::ptrdiff_t v) {
    invariant(v > 0);
    while (true) {
      // get new entering variable
      std::optional<index_t> entering_variable = getEnteringVariable(costs);
      // we break when no costs < 0
      if (!entering_variable) break;
      index_t ev = *entering_variable;
      std::optional<index_t> leave_opt = getLeavingVariable(C, ev);
      // or when no constraints were found
      if (!leave_opt) break;
      makeBasic(C, ev, *leave_opt);
    }
  }
  // Assumes all >v have already been lex-minimized
  // v starts at numVars-1
  // returns `false` if `0`, `true` if not zero
  // minimize v, not touching any variable lex > v
  constexpr auto rLexMin(std::ptrdiff_t v) -> bool {
#ifndef NDEBUG
    assert(in_canonical_form_);
#endif
    index_t c = getBasicConstraint(v);
    if (c < 0) return false;
    if (v == 0) return true;
    MutPtrMatrix<value_type> C{getTableau()};
    // we try to zero `v` or at least minimize it.
    // set cost to 1, and then try to minimize
    //
    // This means that
    // C[0,_] << 0;
    // C[0,1+v] = 1;
    // and then the first step of `run` is to zero the cost of basic
    // variables.
    // We jump to doing this immediately, except that we ignore
    // later `entering var` candidates:
    C[0, _(0, v)] << -C[++c, _(0, v)];
    rLexCore(C, C[0, _(1, v)], v);
    return makeZeroNonBasic(v);
  }
  // get the value of `var`
  [[nodiscard]] constexpr auto getVarValue(std::ptrdiff_t var) const
    -> Rational {
    std::int64_t j = getBasicConstraint(var);
    if (j < 0) return 0;
    PtrMatrix<std::int64_t> constraints = getConstraints();
    return Rational::create(constraints[j, 0], constraints[j, var + 1]);
  }
  /// makeZeroBasic(std::ptrdiff_t v) -> bool
  /// Tries to make `v` non-basic if `v` is zero.
  /// Returns `false` if `v` is zero, `true` otherwise
  constexpr auto makeZeroNonBasic(std::ptrdiff_t v) -> bool {
    utils::assume(v > 0);
    MutPtrMatrix<value_type> C{getTableau()};
    MutPtrVector<index_t> basic_vars{getBasicVariables()};
    MutPtrVector<index_t> basic_constraints{getBasicConstraints()};
    std::int64_t c = basic_constraints[v];
    if (c < 0) return false; // already not basic, hence v is zero
    std::int64_t cc = c++;
    // was not basic
    // v is basic, but not zero
    if (C[c, 0] != 0) return true;
#ifndef NDEBUG
    assertCanonical();
#endif
    // so v is basic and zero.
    // We're going to try to make it non-basic
    // We scan for a variable, trying to find a non-basic one
    std::ptrdiff_t ev = 0;
    for (; (basic_constraints[ev] >= 0) || (C[c, ev + 1] == 0);)
      if (v <= ++ev) return false;
    // we found `ev`

    // we have the invariant that C[c,0] >= 0, but multiplying
    // by `-1` is fine because `C[c,0] == 0`, guaranteed above
    std::ptrdiff_t evp1 = ev + 1;
    if (C[c, evp1] < 0) C[c, _] *= -1;
    for (std::ptrdiff_t i = 1; i < C.numRow(); ++i)
      if (i != std::ptrdiff_t(c))
        NormalForm::zeroWithRowOp(C, row(i), row(c), col(evp1));
    std::int64_t old_basic_var = basic_vars[cc];
    invariant(old_basic_var == std::int64_t(v));
    basic_vars[cc] = index_t(ev);
    basic_constraints[old_basic_var] = -1;
    basic_constraints[ev] = index_t(cc);
#ifndef NDEBUG
    assertCanonical();
#endif
    return false;
  }
  constexpr auto rLexMinLast(std::ptrdiff_t n) -> Solution {
#ifndef NDEBUG
    assert(in_canonical_form_);
    assertCanonical();
#endif
    for (std::ptrdiff_t v = getNumVars(), e = v - n; v != e;) rLexMin(--v);
#ifndef NDEBUG
    assertCanonical();
#endif
    return {.simplex_ = this,
            .skipped_vars_ = length(getNumVars() - n),
            .num_vars_ = length(getNumVars())};
  }
  constexpr auto rLexMinStop(std::ptrdiff_t skip_first,
                             std::ptrdiff_t skip_last = 0) -> Solution {
#ifndef NDEBUG
    assert(in_canonical_form_);
    assertCanonical();
#endif
    std::ptrdiff_t num_v = getNumVars() - skip_last;
    for (std::ptrdiff_t v = num_v; v != skip_first;) rLexMin(--v);
#ifndef NDEBUG
    assertCanonical();
#endif
    return {.simplex_ = this,
            .skipped_vars_ = length(skip_first),
            .num_vars_ = length(num_v)};
  }
  auto maximizeLast(std::ptrdiff_t num_last) -> Rational {
    utils::invariant(num_last > 0);
    MutPtrMatrix<value_type> C{getTableau()};
    std::ptrdiff_t num_v = std::ptrdiff_t(C.numCol());
    C[0, _(0, num_v - num_last)] << 0;
    C[0, _(num_v - num_last, num_v)] << -1;
    return run();
  }
  auto maximizeLastDrop(std::ptrdiff_t num_last) -> Rational {
    utils::invariant(num_last > 0);
    MutPtrMatrix<value_type> C{getTableau()};
    std::ptrdiff_t num_v = std::ptrdiff_t(C.numCol()) - 1;
    Rational f = maximizeLast(num_last);
    MutPtrVector<index_t> basic_cons{getBasicConstraints()},
      basic_vars{getBasicVariables()};
    std::ptrdiff_t last_row = std::ptrdiff_t(C.numRow());
    // now we must drop `num_last`
    // if they are basic, we subtract them from the constraints,
    // and drop the associated constraint
    std::ptrdiff_t v = num_v, num_r = num_v - num_last;
    do {
      --v;
      index_t c = basic_cons[v];
      if (c < 0) continue; // was 0
                           // std::int64_t n = C[c + 1, 0], d = C[c + 1, v + 1];
      // #ifndef NDEBUG
      //       if (!((C[_(last_row), v + 1] * n / d) * d)
      //              .isEqual(C[_(last_row), v + 1] * n))
      //         __builtin_trap();
      // #endif
      //       if (n) C[_(last_row), 0] -= (C[_(last_row), v + 1] * n) / d;
      C[c + 1, 0] = 0;
      if (--last_row == c + 1) continue;
      C[c + 1, _(v + 1)] << C[last_row, _(v + 1)];
      index_t last_row_var = basic_vars[last_row - 1];
      utils::invariant(basic_vars[c] == v);
      basic_vars[c] = last_row_var;
      basic_cons[last_row_var] = c;
    } while (v != num_r);
    truncateVars(num_r);
    truncateConstraints(--last_row);
    return f;
  }

  // reverse lexicographic ally minimize vars
  constexpr void rLexMin(Vector<Rational> &sol) {
    sol << rLexMinLast(sol.size());
  }
  // A(:,1:end)*x <= A(:,0)
  // B(:,1:end)*x == B(:,0)
  // returns a Simplex if feasible, and an empty `Optional` otherwise
  static constexpr auto positiveVariables(alloc::Arena<> *alloc,
                                          PtrMatrix<std::int64_t> A,
                                          PtrMatrix<std::int64_t> B)
    -> Optional<Simplex *> {
    invariant(A.numCol() == B.numCol());
    std::ptrdiff_t num_var = std::ptrdiff_t(A.numCol()) - 1,
                   num_slack = std::ptrdiff_t(A.numRow()),
                   num_strict = std::ptrdiff_t(B.numRow()),
                   num_con = num_slack + num_strict,
                   var_cap = num_var + num_slack;
    // see how many slack vars are infeasible as solution
    // each of these will require an augment variable
    for (std::ptrdiff_t i = 0; i < num_slack; ++i) var_cap += A[i, 0] < 0;
    // try to avoid reallocating
    auto checkpoint{alloc->checkpoint()};
    Simplex *simplex{Simplex::create(alloc, row(num_con),
                                     col(num_var + num_slack),
                                     capacity(num_con), stride(var_cap))};
    // construct:
    // [ I A
    //   0 B ]
    // then drop the extra variables
    slackEqualityConstraints(simplex->getConstraints()[_, _(1, end)],
                             A[_, _(1, end)], B[_, _(1, end)]);
    auto consts{simplex->getConstants()};
    consts[_(0, num_slack)] << A[_, 0];
    if (num_strict) consts[_(num_slack, num_slack + num_strict)] << B[_, 0];
    // for (std::ptrdiff_t i = 0; i < numSlack; ++i) consts[i] = A(i, 0);
    // for (std::ptrdiff_t i = 0; i < numStrict; ++i) consts[i + numSlack] =
    // B(i, 0);
    if (!simplex->initiateFeasible()) return simplex;
    alloc->rollback(checkpoint);
    return nullptr;
  }
  static constexpr auto positiveVariables(alloc::Arena<> *alloc,
                                          PtrMatrix<std::int64_t> A)
    -> Optional<Simplex *> {
    std::ptrdiff_t num_var = std::ptrdiff_t(A.numCol()) - 1,
                   num_slack = std::ptrdiff_t(A.numRow()), num_con = num_slack,
                   var_cap = num_var + num_slack;
    // see how many slack vars are infeasible as solution
    // each of these will require an augment variable
    for (std::ptrdiff_t i = 0; i < num_slack; ++i) var_cap += A[i, 0] < 0;
    // try to avoid reallocating
    auto checkpoint{alloc->checkpoint()};
    Simplex *simplex{Simplex::create(alloc, row(num_con),
                                     col(num_var + num_slack),
                                     capacity(num_con), stride(var_cap))};
    // construct:
    // [ I A ]
    // then drop the extra variables
    slackEqualityConstraints(simplex->getConstraints()[_, _(1, end)],
                             A[_, _(1, end)]);
    // auto consts{simplex.getConstants()};
    // for (std::ptrdiff_t i = 0; i < numSlack; ++i) consts[i] = A(i, 0);
    simplex->getConstants() << A[_, 0];
    if (!simplex->initiateFeasible()) return simplex;
    alloc->rollback(checkpoint);
    return nullptr;
  }

  constexpr void pruneBounds(alloc::Arena<> *alloc,
                             std::ptrdiff_t numSlack = 0) {
    auto p = alloc->scope();
    Simplex *simplex{Simplex::create(alloc, num_constraints_, num_vars_,
                                     constraint_capacity_,
                                     --auto{var_capacity_p1_})};
    // Simplex simplex{getNumCons(), getNumVars(), getNumSlack(), 0};
    for (std::ptrdiff_t c = 0; c < getNumCons(); ++c) {
      *simplex << *this;
      MutPtrMatrix<std::int64_t> constraints = simplex->getConstraints();
      std::int64_t bumped_bound = ++constraints[c, 0];
      MutPtrVector<std::int64_t> cost = simplex->getCost();
      for (std::ptrdiff_t v = numSlack; v < cost.size(); ++v)
        cost[v] = -constraints[c, v + 1];
      if (simplex->run() != bumped_bound) deleteConstraint(c--);
    }
  }

  // static constexpr auto toMask(PtrVector<std::int64_t> x) -> std::uint64_t {
  //   assert(x.size() <= 64);
  //   std::uint64_t m = 0;
  //   for (auto y : x) m = ((m << 1) | (y != 0));
  //   return m;
  // }
  // [[nodiscard]] constexpr auto getBasicTrueVarMask() const -> std::uint64_t {
  //   const std::ptrdiff_t numVarTotal = getNumVars();
  //   assert(numVarTotal <= 64);
  //   std::uint64_t m = 0;
  //   PtrVector<index_t> basicCons{getBasicConstraints()};
  //   for (std::ptrdiff_t i = numSlack; i < numVarTotal; ++i)
  //     m = ((m << 1) | (basicCons[i] > 0));
  //   return m;
  // }
  // check if a solution exists such that `x` can be true.
  // returns `true` if unsatisfiable
  [[nodiscard]] constexpr auto unSatisfiable(alloc::Arena<> alloc,
                                             PtrVector<std::int64_t> x,
                                             std::ptrdiff_t off) const -> bool {
    // is it a valid solution to set the first `x.size()` variables to
    // `x`? first, check that >= 0 constraint is satisfied
    if (!allGEZero(x)) return true;
    // approach will be to move `x.size()` variables into the
    // equality constraints, and then check if the remaining sub-problem
    // is satisfiable.
    const std::ptrdiff_t num_con = getNumCons(), num_var = getNumVars(),
                         num_fix = x.size();
    Simplex *sub_simp{
      Simplex::create(&alloc, row(num_con), col(num_var - num_fix))};
    // subSimp.tableau(0, 0) = 0;
    // subSimp.tableau(0, 1) = 0;
    auto fC{getTableau()};
    auto sC{sub_simp->getTableau()};
    sC[_, 0] << fC[_, 0] - fC[_, _(1 + off, 1 + off + num_fix)] * x.t();
    // sC(_, 0) = fC(_, 0);
    // for (std::ptrdiff_t i = 0; i < numFix; ++i)
    //     sC(_, 0) -= x(i) * fC(_, i + 1 + off);
    sC[_, _(1, 1 + off)] << fC[_, _(1, 1 + off)];
    sC[_, _(1 + off, end)] << fC[_, _(1 + off + num_fix, end)];
    // returns `true` if unsatisfiable
    return sub_simp->initiateFeasible();
  }
  [[nodiscard]] constexpr auto satisfiable(alloc::Arena<> alloc,
                                           PtrVector<std::int64_t> x,
                                           std::ptrdiff_t off) const -> bool {
    return !unSatisfiable(alloc, x, off);
  }
  // check if a solution exists such that `x` can be true.
  // zeros remaining rows
  [[nodiscard]] constexpr auto
  unSatisfiableZeroRem(alloc::Arena<> alloc, PtrVector<std::int64_t> x,
                       std::ptrdiff_t off, std::ptrdiff_t numRow) const
    -> bool {
    // is it a valid solution to set the first `x.size()` variables to
    // `x`? first, check that >= 0 constraint is satisfied
    if (!allGEZero(x)) return true;
    // approach will be to move `x.size()` variables into the
    // equality constraints, and then check if the remaining sub-problem
    // is satisfiable.
    invariant(numRow <= getNumCons());
    const std::ptrdiff_t num_fix = x.size();
    Simplex *sub_simp{Simplex::create(&alloc, row(numRow), col(off++))};
    auto fC{getConstraints()};
    auto sC{sub_simp->getConstraints()};
    sC[_, 0] << fC[_(begin, numRow), 0] -
                  fC[_(begin, numRow), _(off, off + num_fix)] * x.t();
    sC[_, _(1, off)] << fC[_(begin, numRow), _(1, off)];
    return sub_simp->initiateFeasible();
  }
  /// indsFree gives how many variables are free to take  any >= 0 value
  /// indOne is var ind greater than indsFree that's pinned to 1
  /// (i.e., indsFree + indOne == index of var pinned to 1)
  /// numRow is number of rows used, extras are dropped
  // [[nodiscard]] constexpr auto
  [[nodiscard]] auto unSatisfiableZeroRem(alloc::Arena<> alloc,
                                          std::ptrdiff_t iFree,
                                          std::array<std::ptrdiff_t, 2> inds,
                                          std::ptrdiff_t numRow) const -> bool {
    invariant(numRow <= getNumCons());
    Simplex *sub_simp{Simplex::create(&alloc, row(numRow), col(iFree++))};
    auto fC{getConstraints()};
    auto sC{sub_simp->getConstraints()};
    auto r = _(0, numRow);
    sC[_, 0] << fC[r, 0] - (fC[r, inds[0] + iFree] + fC[r, inds[1] + iFree]);
    sC[_, _(1, iFree)] << fC[r, _(1, iFree)];
    return sub_simp->initiateFeasible();
  }
  [[nodiscard]] constexpr auto
  satisfiableZeroRem(alloc::Arena<> alloc, PtrVector<std::int64_t> x,
                     std::ptrdiff_t off, std::ptrdiff_t numRow) const -> bool {
    return !unSatisfiableZeroRem(alloc, x, off, numRow);
  }
  void printResult(std::ptrdiff_t numSlack = 0) {
    auto C{getConstraints()};
    auto basic_vars{getBasicVariables()};
    for (std::ptrdiff_t i = 0; i < basic_vars.size(); ++i) {
      std::ptrdiff_t v = basic_vars[i];
      if (v <= numSlack) continue;
      if (C[i, 0]) {
        if (++v < C.numCol()) {
          std::cout << "v_" << v - numSlack << " = " << C[i, 0] << " / "
                    << C[i, v] << "\n";
        } else {
          std::cout << "v_" << v << " = " << C[i, 0] << "\n";
          __builtin_trap();
        }
      }
    }
  }
  static constexpr auto create(alloc::Arena<> *alloc, Row<> numCon,
                               Col<> numVar) -> Valid<Simplex> {
    return create(alloc, numCon, numVar, capacity(std::ptrdiff_t(numCon)),
                  stride(std::ptrdiff_t(numVar) + std::ptrdiff_t(numCon)));
  }
  static constexpr auto create(alloc::Arena<> *alloc, Row<> numCon,
                               Col<> numVar, Capacity<> conCap,
                               RowStride<> varCap) -> Valid<Simplex> {
    varCap = alignVarCapacity(varCap);
    auto c_cap = std::ptrdiff_t(conCap), v_cap = std::ptrdiff_t(varCap);
    std::size_t mem_needed = requiredMemory(c_cap, v_cap);
    auto *mem = (Simplex *)alloc->allocate<alignof(Simplex)>(mem_needed);
    mem->num_constraints_ = numCon;
    mem->num_vars_ = numVar;
    mem->constraint_capacity_ = conCap;
    mem->var_capacity_p1_ = varCap;
    return mem;
  }

  static auto operator new(std::size_t count, Capacity<> conCap,
                           RowStride<> varCap) -> void * {
    auto cC = std::ptrdiff_t(conCap), vC = std::ptrdiff_t(varCap);
    std::size_t mem_needed = requiredMemory(cC, vC);
    // void *p = ::operator new(count * memNeeded);
    return alloc::malloc(count * mem_needed,
                         std::align_val_t(alignof(Simplex)));
    // return ::operator new(count * memNeeded,
    //                       std::align_val_t(alignof(Simplex)));
  }
  static void operator delete(void *ptr, std::size_t sz) {
    alloc::free(ptr, sz, std::align_val_t(alignof(Simplex)));
    // ::operator delete(ptr, std::align_val_t(alignof(Simplex)));
  }

  static auto create(Row<> numCon, Col<> numVar) -> std::unique_ptr<Simplex> {
    auto nc = std::ptrdiff_t(numCon);
    return create(numCon, numVar, capacity(nc),
                  stride(std::ptrdiff_t(numVar) + nc));
  }
  static auto create(Row<> numCon, Col<> numVar, Capacity<> conCap,
                     RowStride<> varCap) -> std::unique_ptr<Simplex> {
    varCap = alignVarCapacity(varCap);
    auto *ret = new (conCap, varCap) Simplex;
    ret->num_constraints_ = numCon;
    ret->num_vars_ = numVar;
    ret->constraint_capacity_ = conCap;
    ret->var_capacity_p1_ = varCap;
    return std::unique_ptr<Simplex>(ret);
  }

  static constexpr auto
  create(alloc::Arena<> *alloc, Row<> numCon,
         Col<> numVar, // NOLINT(bugprone-easily-swappable-parameters)
         std::ptrdiff_t numSlack) -> Valid<Simplex> {
    std::ptrdiff_t con_cap = std::ptrdiff_t(numCon),
                   var_cap = std::ptrdiff_t(numVar) + numSlack + con_cap;
    return create(alloc, numCon, numVar, capacity(con_cap), stride(var_cap));
  }
  constexpr auto copy(alloc::Arena<> *alloc) const -> Valid<Simplex> {
    Valid<Simplex> res = create(alloc, row(getNumCons()), col(getNumVars()),
                                getConCap(), getVarCap());
    *res << *this;
    return res;
  }
  constexpr auto operator<<(const Simplex &other) -> Simplex & {
    setNumCons(other.getNumCons());
    setNumVars(other.getNumVars());
    getTableau() << other.getTableau();
    getBasicVariables() << other.getBasicVariables();
    getBasicConstraints() << other.getBasicConstraints();
    return *this;
  }
  friend auto operator<<(std::ostream &os, const Simplex &s) -> std::ostream & {
    os << "Basic Variables: " << s.getBasicVariables();
    os << "Basic Constraints: " << s.getBasicConstraints();
    os << "Constraints:\n" << s.getConstraints();
    return os;
  }
#ifndef NDEBUG
  [[gnu::used]] void dump() const { std::cout << *this; }
#endif
};

static_assert(AbstractVector<Simplex::Solution>);

static_assert(AbstractVector<PtrVector<Rational>>);
// static_assert(AbstractVector<ElementwiseBinaryOp<
//                 PtrVector<Rational>, PtrVector<Rational>, std::minus<>>>);
static_assert(std::movable<Simplex::Solution::iterator>);
static_assert(std::indirectly_readable<Simplex::Solution::iterator>);
static_assert(std::forward_iterator<Simplex::Solution::iterator>);
static_assert(alignof(Simplex) == simd::VECTORWIDTH);
} // namespace math

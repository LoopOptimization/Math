#pragma once
#include "Alloc/Arena.hpp"
#include "Math/Array.hpp"
#include "Math/AxisTypes.hpp"
#include "Math/Constraints.hpp"
#include "Math/Indexing.hpp"
#include "Math/Math.hpp"
#include "Math/MatrixDimensions.hpp"
#include "Math/NormalForm.hpp"
#include "Math/Rational.hpp"
#include "SIMD/Intrin.hpp"
#include "Utilities/Invariant.hpp"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <ostream>

namespace poly::math {
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
  using index_type = int;
  using value_type = int64_t;

  template <std::integral T>
  static constexpr auto alignOffset(ptrdiff_t x) -> ptrdiff_t {
    --x;
    ptrdiff_t W = simd::Width<T>;
    ptrdiff_t nW = -W;
    x += W;
    x &= nW;
    return x;
    // return (--x + simd::Width<T>)&(-simd::Width<T>);
  }
  [[gnu::returns_nonnull, nodiscard]] constexpr auto
  basicConsPointer() const -> index_type * {
    void *p = const_cast<char *>(memory);
    return std::assume_aligned<simd::VECTORWIDTH>(static_cast<index_type *>(p));
  }
  [[gnu::returns_nonnull, nodiscard]] constexpr auto
  basicVarsPointer() const -> index_type * {
    ptrdiff_t offset = alignOffset<index_type>(reservedBasicConstraints());
    return std::assume_aligned<simd::VECTORWIDTH>(basicConsPointer() + offset);
  }
  // offset in bytes
  static constexpr auto tableauOffset(ptrdiff_t cons,
                                      ptrdiff_t vars) -> ptrdiff_t {
    ptrdiff_t coff = alignOffset<index_type>(cons);
    ptrdiff_t voff = alignOffset<index_type>(vars);
    ptrdiff_t offset = coff + voff;
    // ptrdiff_t offset =
    //   alignOffset<index_type>(cons) + alignOffset<index_type>(vars);
    return static_cast<ptrdiff_t>(sizeof(index_type)) * offset;
  }
  [[gnu::returns_nonnull, nodiscard]] constexpr auto
  tableauPointer() const -> value_type * {
    ptrdiff_t offset =
      tableauOffset(reservedBasicConstraints(), reservedBasicVariables());
    void *p = const_cast<char *>(memory) + offset;
    return std::assume_aligned<simd::VECTORWIDTH>(static_cast<value_type *>(p));
  }
  Row<> numConstraints{};
  Col<> numVars{};
  Capacity<> constraintCapacity;
  RowStride<> varCapacityP1; // varCapacity + 1
#ifndef NDEBUG
  bool inCanonicalForm{false};
#endif
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#else
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) // FAM
  alignas(simd::VECTORWIDTH) char memory[];
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#else
#pragma clang diagnostic pop
#endif

  // (varCapacity+1)%simd::Width<int64_t>==0
  static constexpr auto alignVarCapacity(RowStride<> rs) -> RowStride<> {
    static constexpr ptrdiff_t W = simd::Width<int64_t>;
    return rowStride((ptrdiff_t(rs) + W) & -W);
  }
  [[nodiscard]] static constexpr auto
  reservedTableau(ptrdiff_t cons, ptrdiff_t vars) -> ptrdiff_t {
    return static_cast<ptrdiff_t>(sizeof(value_type)) *
           ((ptrdiff_t(cons) + 1) * vars);
  }
  static constexpr auto requiredMemory(ptrdiff_t cons,
                                       ptrdiff_t vars) -> size_t {
    ptrdiff_t base = static_cast<ptrdiff_t>(sizeof(Simplex)),
              indices = tableauOffset(cons, vars),
              tableau = reservedTableau(cons, vars);
    return static_cast<size_t>(base + indices + tableau);
  }

public:
  // tableau is constraint * var matrix w/ extra col for LHS
  // and extra row for objective function
  [[nodiscard]] constexpr auto reservedBasicConstraints() const -> ptrdiff_t {
    return ptrdiff_t(varCapacityP1) - 1;
  }
  [[nodiscard]] constexpr auto reservedBasicVariables() const -> ptrdiff_t {
    return ptrdiff_t(constraintCapacity);
  }

  // [[nodiscard]] constexpr auto intsNeeded() const -> ptrdiff_t {
  //   return reservedTableau() + reservedBasicConstraints() +
  //          reservedBasicVariables();
  // }
  // We with to align every row of `getConstraints()` and `getTableau()`.
  // To do this, we
  // 0. (varCapacity+1) % simd::Width<int64_t> == 0
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
                                ++auto(numConstraints),
                                ++auto(numVars),
                                varCapacityP1,
                              }};
  }
  // NOLINTNEXTLINE(readability-make-member-function-const)
  [[nodiscard]] constexpr auto getTableau() -> MutPtrMatrix<value_type> {
    return {tableauPointer(), StridedDims<>{
                                ++auto(numConstraints),
                                ++auto(numVars),
                                varCapacityP1,
                              }};
  }
  [[nodiscard]] constexpr auto getConstraints() const -> PtrMatrix<value_type> {
    return {tableauPointer() + ptrdiff_t(varCapacityP1), StridedDims<>{
                                                           numConstraints,
                                                           ++auto(numVars),
                                                           varCapacityP1,
                                                         }};
  }
  // NOLINTNEXTLINE(readability-make-member-function-const)
  [[nodiscard]] constexpr auto getConstraints() -> MutPtrMatrix<value_type> {
    return {tableauPointer() + ptrdiff_t(varCapacityP1), StridedDims<>{
                                                           numConstraints,
                                                           ++auto(numVars),
                                                           varCapacityP1,
                                                         }};
  }
  [[nodiscard]] constexpr auto
  getBasicConstraints() const -> PtrVector<index_type> {
    return {basicConsPointer(), aslength(numVars)};
  }
  [[nodiscard]] constexpr auto
  getBasicConstraints() -> MutPtrVector<index_type> {
    return {basicConsPointer(), aslength(numVars)};
  }
  [[nodiscard]] constexpr auto
  getBasicVariables() const -> PtrVector<index_type> {
    return {basicVarsPointer(), length(ptrdiff_t(numConstraints))};
  }
  [[nodiscard]] constexpr auto getBasicVariables() -> MutPtrVector<index_type> {
    return {basicVarsPointer(), length(ptrdiff_t(numConstraints))};
  }
  [[nodiscard]] constexpr auto getCost() const -> PtrVector<value_type> {
    return {tableauPointer(), length(ptrdiff_t(numVars) + 1z)};
  }
  // NOLINTNEXTLINE(readability-make-member-function-const)
  [[nodiscard]] constexpr auto getCost() -> MutPtrVector<value_type> {
    return {tableauPointer(), length(ptrdiff_t(numVars) + 1z)};
  }
  [[nodiscard]] constexpr auto
  getBasicConstraint(ptrdiff_t i) const -> index_type {
    return getBasicConstraints()[i];
  }
  [[nodiscard]] constexpr auto
  getBasicVariable(ptrdiff_t i) const -> index_type {
    return getBasicVariables()[i];
  }
  [[nodiscard]] constexpr auto
  getObjectiveCoefficient(ptrdiff_t i) const -> value_type {
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
    inCanonicalForm = false;
#endif
    auto C{getConstraints()};
#ifndef NDEBUG
    for (ptrdiff_t r = 0, R = ptrdiff_t(numRows(C)); r < R; ++r)
      for (ptrdiff_t c = 0, N = ptrdiff_t(numCols(C)); c < N; ++c)
        invariant(C[r, c] != std::numeric_limits<int64_t>::min());
#endif
    NormalForm::solveSystemSkip(C);
    truncateConstraints(ptrdiff_t(NormalForm::numNonZeroRows(C)));
  }
#ifndef NDEBUG
  constexpr void assertCanonical() const {
    PtrMatrix<value_type> C{getTableau()};
    PtrVector<index_type> basicVars{getBasicVariables()};
    PtrVector<index_type> basicCons{getBasicConstraints()};
    for (ptrdiff_t v = 0; v < basicCons.size();) {
      index_type c = basicCons[v++];
      if (c < 0) continue;
      if (!allZero(C[_(1, 1 + c), v])) __builtin_trap();
      if (!allZero(C[_(2 + c, end), v])) __builtin_trap();
      if (ptrdiff_t(basicVars[c]) != (v - 1)) __builtin_trap();
    }
    for (ptrdiff_t c = 1; c < C.numRow(); ++c) {
      index_type v = basicVars[c - 1];
      if (ptrdiff_t(v) < basicCons.size()) {
        invariant(c - 1, ptrdiff_t(basicCons[v]));
        invariant(C[c, v + 1] >= 0);
      }
      invariant(C[c, 0] >= 0);
    }
  }
#endif
  [[nodiscard]] constexpr auto getConstants() -> MutStridedVector<int64_t> {
    return getTableau()[_(1, end), 0];
  }
  [[nodiscard]] constexpr auto getConstants() const -> StridedVector<int64_t> {
    return getTableau()[_(1, end), 0];
  }
  constexpr void truncateConstraints(ptrdiff_t i) {
    invariant(ptrdiff_t(numConstraints) <= ptrdiff_t(constraintCapacity));
    invariant(i >= 0z);
    invariant(i <= numConstraints);
    numConstraints = row(i);
  }
  constexpr void setNumCons(ptrdiff_t i) {
    invariant(i <= constraintCapacity);
    numConstraints = row(i);
  }
  constexpr void setNumVars(ptrdiff_t i) {
    invariant(i < varCapacityP1);
    numVars = col(i);
  }
  constexpr void truncateVars(ptrdiff_t i) {
    invariant(i <= numVars);
    numVars = col(i);
  }
  [[nodiscard]] constexpr auto getNumCons() const -> ptrdiff_t {
    invariant(numConstraints >= 0);
    return ptrdiff_t(numConstraints);
  }
  [[nodiscard]] constexpr auto getNumVars() const -> ptrdiff_t {
    invariant(numVars >= 0);
    return ptrdiff_t(numVars);
  }
  [[nodiscard]] constexpr auto getConCap() const -> Capacity<> {
    invariant(constraintCapacity >= 0);
    return constraintCapacity;
  }
  [[nodiscard]] constexpr auto getVarCap() const -> RowStride<> {
    invariant(varCapacityP1 > 0);
    return --auto{varCapacityP1};
  }
  constexpr void deleteConstraint(ptrdiff_t c) {
    auto basicCons = getBasicConstraints();
    auto basicVars = getBasicVariables();
    auto constraints = getConstraints();
    --numConstraints;
    if (auto basicVar = basicVars[c]; basicVar >= 0) basicCons[basicVar] = -1;
    if (c == numConstraints) return;
    auto basicVar = basicVars[ptrdiff_t(numConstraints)];
    basicVars[c] = basicVar;
    if (basicVar >= 0) basicCons[basicVar] = index_type(c);
    constraints[c, _] << constraints[numConstraints, _];
  }

  // AbstractVector
  struct Solution {
    using value_type = Rational;
    // view of tableau dropping const column
    Valid<const Simplex> simplex;
    Length<> skippedVars;
    Length<> numVars;
    class iterator { // NOLINT(readability-identifier-naming)
      const Solution *sol;
      ptrdiff_t i;

    public:
      using value_type = Rational;
      constexpr iterator(const Solution *s, ptrdiff_t j) : sol(s), i(j) {}
      constexpr iterator() = default;
      constexpr iterator(const iterator &) = default;
      constexpr auto operator=(const iterator &) -> iterator & = default;
      auto operator*() const -> Rational { return (*sol)[i]; }
      constexpr auto operator++() -> iterator & {
        ++i;
        return *this;
      }
      constexpr auto operator++(int) -> iterator {
        auto tmp = *this;
        ++i;
        return tmp;
      }
      constexpr auto operator--() -> iterator & {
        --i;
        return *this;
      }
      constexpr auto operator--(int) -> iterator {
        auto tmp = *this;
        --i;
        return tmp;
      }
      friend constexpr auto operator==(iterator a, iterator b) -> bool {
        return a.i == b.i;
      }
      friend constexpr auto operator!=(iterator a, iterator b) -> bool {
        return a.i != b.i;
      }
      constexpr auto operator-(iterator b) const -> ptrdiff_t {
        return ptrdiff_t(i) - ptrdiff_t(b.i);
      }
      constexpr auto operator+(ptrdiff_t n) const -> iterator {
        return {sol, i + ptrdiff_t(n)};
      }
    };
    [[nodiscard]] constexpr auto begin() const -> iterator { return {this, 0}; }
    [[nodiscard]] constexpr auto end() const -> iterator {
      return {this, ptrdiff_t(numVars - skippedVars)};
    }

    [[nodiscard]] constexpr auto operator[](ptrdiff_t i) const -> Rational {
      invariant(ptrdiff_t(i) >= 0);
      i += ptrdiff_t(skippedVars);
      int64_t j = simplex->getBasicConstraint(i);
      if (j < 0) return 0;
      PtrMatrix<int64_t> constraints = simplex->getConstraints();
      return Rational::create(constraints[j, 0], constraints[j, i + 1]);
    }
    [[nodiscard]] constexpr auto operator[](OffsetEnd k) const -> Rational {
      ptrdiff_t i = ptrdiff_t(simplex->numVars) - k.offset;
      int64_t j = simplex->getBasicConstraint(i);
      if (j < 0) return 0;
      PtrMatrix<int64_t> constraints = simplex->getConstraints();
      return Rational::create(constraints[j, 0], constraints[j, i + 1]);
    }
    [[nodiscard]] constexpr auto
    operator[](ScalarRelativeIndex auto i) const -> Rational {
      return (*this)[calcOffset(size(), i)];
    }
    template <typename B, typename E>
    constexpr auto operator[](Range<B, E> r) const -> Solution {
      return (*this)[canonicalizeRange(r, size())];
    }
    constexpr auto operator[](Range<ptrdiff_t, ptrdiff_t> r) const -> Solution {
      return {simplex, length(ptrdiff_t(skippedVars) + r.b),
              length(ptrdiff_t(skippedVars) + r.e)};
    }
    [[nodiscard]] constexpr auto size() const -> ptrdiff_t {
      return ptrdiff_t(numVars - skippedVars);
    }
    [[nodiscard]] constexpr auto view() const -> Solution { return *this; };

    [[nodiscard]] constexpr auto denomLCM() const -> int64_t {
      int64_t l = 1;
      for (auto r : *this) l = lcm(l, r.denominator);
      return l;
    }
    friend inline auto operator<<(std::ostream &os,
                                  Solution sol) -> std::ostream & {
      os << "Simplex::Solution[";
      bool printComma = false;
      for (auto b : sol) {
        if (printComma) os << ", ";
        printComma = true;
        os << b;
      }
      os << "]";
      return os;
    }
#ifndef NDEBUG
    [[gnu::used]] void dump() const { std::cout << *this; }
#endif
  };
  [[nodiscard]] constexpr auto getSolution() const -> Solution {
    return {this, {}, aslength(numVars)};
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
    const auto numVar = ptrdiff_t(getNumVars());
    MutPtrMatrix<value_type> C{getConstraints()};
    MutPtrVector<index_type> basicCons{getBasicConstraints()};
    basicCons << -2;
    // first pass, we make sure the equalities are >= 0
    // and we eagerly try and find columns with
    // only a single non-0 element.
    for (ptrdiff_t c = 0; c < C.numRow(); ++c) {
      int64_t &Ceq = C[c, 0];
      int64_t sign = 2 * (Ceq >= 0) - 1;
      Ceq *= sign;
      for (ptrdiff_t v = 0; v < numVar; ++v)
        if (int64_t Ccv = C[c, v + 1] *= sign)
          basicCons[v] =
            (((basicCons[v] == -2) && (Ccv > 0))) ? index_type(c) : -1;
    }
    // basicCons should now contain either `-1` or an integer >= 0
    // indicating which row contains the only non-zero element; we'll
    // now fill basicVars.
    //
    auto basicVars{getBasicVariables()};
    basicVars << -1;
    for (ptrdiff_t v = 0; v < numVar; ++v) {
      if (int64_t r = basicCons[v]; r >= 0) {
        if (basicVars[r] == -1) basicVars[r] = index_type(v);
        else basicCons[v] = -1;
      }
    }
#ifndef NDEBUG
    inCanonicalForm = true;
#endif
    Vector<unsigned> augVars{};
    // upper bound number of augmentVars is constraintCapacity
    for (ptrdiff_t i = 0; i < basicVars.size(); ++i)
      if (basicVars[i] == -1) augVars.push_back(i);
    return (!augVars.empty() && removeAugmentVars(augVars));
  }
  constexpr auto removeAugmentVars(PtrVector<unsigned> augmentVars) -> bool {
    // TODO: try to avoid reallocating, via reserving enough ahead of time
    ptrdiff_t numAugment = augmentVars.size(), oldNumVar = ptrdiff_t(numVars);
    invariant(numAugment + ptrdiff_t(numVars) < varCapacityP1);
    numVars = col(ptrdiff_t(numVars) + ptrdiff_t(numAugment));
    MutPtrMatrix<value_type> C{getConstraints()};
    MutPtrVector<index_type> basicVars{getBasicVariables()};
    MutPtrVector<index_type> basicCons{getBasicConstraints()};
    MutPtrVector<value_type> costs{getCost()};
    costs << 0;
    C[_, _(oldNumVar + 1, end)] << 0;
    for (ptrdiff_t i = 0; i < ptrdiff_t(augmentVars.size()); ++i) {
      ptrdiff_t a = augmentVars[i];
      basicVars[a] = index_type(i) + index_type(oldNumVar);
      basicCons[i + oldNumVar] = index_type(a);
      C[a, oldNumVar + 1 + i] = 1;
      // we now zero out the implicit cost of `1`
      costs[_(begin, oldNumVar + 1)] -= C[a, _(begin, oldNumVar + 1)];
    }
    ASSERT(std::all_of(basicVars.begin(), basicVars.end(),
                       [](int64_t i) { return i >= 0; }));
    // false/0 means feasible
    // true/non-zero infeasible
    if (runCore()) return true;
    // check for any basic vars set to augment vars, and set them to some
    // other variable (column) instead.
    for (ptrdiff_t c = 0; c < C.numRow(); ++c) {
      if (ptrdiff_t(basicVars[c]) >= ptrdiff_t(oldNumVar)) {
        invariant(C[c, 0] == 0);
        invariant(c == basicCons[basicVars[c]]);
        invariant(C[c, basicVars[c] + 1] >= 0);
        // find var to make basic in its place
        for (ptrdiff_t v = oldNumVar; v != 0;) {
          // search for a non-basic variable
          // (basicConstraints<0)
          int64_t Ccv = C[c, v--];
          if (Ccv == 0 || (basicCons[v] >= 0)) continue;
          if (Ccv < 0) C[c, _] *= -1;
          for (ptrdiff_t i = 0; i < C.numRow(); ++i)
            if (i != ptrdiff_t(c))
              NormalForm::zeroWithRowOp(C, row(i), row(c), ++col(v), 0);
          basicVars[c] = index_type(v);
          basicCons[v] = index_type(c);
          break;
        }
      }
    }
    // all augment vars are now 0
    numVars = col(oldNumVar);
#ifndef NDEBUG
    assertCanonical();
#endif
    return false;
  }

  // 1 based to match getBasicConstraints
  [[nodiscard]] static constexpr auto
  getEnteringVariable(PtrVector<int64_t> costs) -> Optional<int> {
    // Bland's algorithm; guaranteed to terminate
    auto f = costs.begin(), l = costs.end();
    const auto *neg = std::find_if(f, l, [](int64_t c) { return c < 0; });
    if (neg == l) return {};
    return int(std::distance(f, neg));
  }
  [[nodiscard]] static constexpr auto
  getLeavingVariable(PtrMatrix<int64_t> C,
                     ptrdiff_t enteringVariable) -> Optional<unsigned int> {
    // inits guarantee first valid is selected
    int64_t n = -1, d = 0;
    unsigned int j = 0;
    for (ptrdiff_t i = 1; i < C.numRow(); ++i) {
      int64_t Civ = C[i, enteringVariable + 1];
      if (Civ <= 0) continue;
      int64_t Cio = C[i, 0];
      if (Cio == 0) return --i;
      invariant(Cio > 0);
      if ((n * Cio) >= (Civ * d)) continue;
      // we could consider something like:
      // so that in case of ties, we prefer having the maximum index
      // be the leaving variable.
      // That didn't really help the existing benchmarks.
      // auto basicVars = getBasicVariables(); // passed in as arg to static fun
      // if ((n * Cio) > (Civ * d)) continue;
      // if ((n * Cio) == (Civ * d) && (basicVars[i - 1] < basicVars[j - 1]))
      //   continue;
      n = Civ;
      d = Cio;
      j = i;
    }
    // NOTE: if we fail to find a leaving variable, then `j = 0`,
    // and it will unsigned wrap to `ptrdiff_t(-1)`, which indicates
    // an empty `Optional<unsigned int>`
    return --j;
  }
  constexpr auto makeBasic(MutPtrMatrix<int64_t> C, int64_t f,
                           int enteringVar) -> int64_t {
    Optional<unsigned int> leaveOpt = getLeavingVariable(C, enteringVar);
    if (!leaveOpt) return 0; // unbounded
    auto leavingVar = ptrdiff_t(*leaveOpt);
    for (ptrdiff_t i = 0; i < C.numRow(); ++i) {
      if (i == leavingVar + 1) continue;
      int64_t m = NormalForm::zeroWithRowOp(C, row(i), ++row(leavingVar),
                                            ++col(enteringVar), i ? 0 : f);
      if (!i) f = m;
    }
    // update basic vars and constraints
    MutPtrVector<index_type> basicVars{getBasicVariables()};
    int64_t oldBasicVar = basicVars[leavingVar];
    basicVars[leavingVar] = enteringVar;
    MutPtrVector<index_type> basicConstraints{getBasicConstraints()};
    basicConstraints[oldBasicVar] = -1;
    basicConstraints[enteringVar] = index_type(leavingVar);
    return f;
  }
  // run the simplex algorithm, assuming basicVar's costs have been set to
  // 0
  constexpr auto runCore(int64_t f = 1) -> Rational {
#ifndef NDEBUG
    ASSERT(inCanonicalForm);
#endif
    //     return runCore(getCostsAndConstraints(), f);
    // }
    // Rational runCore(MutPtrMatrix<int64_t> C, int64_t f = 1) {
    MutPtrMatrix<int64_t> C{getTableau()};
    do {
      // entering variable is the column
      Optional<int> enteringVariable = getEnteringVariable(C[0, _(1, end)]);
      if (!enteringVariable) return Rational::create(C[0, 0], f);
      f = makeBasic(C, f, *enteringVariable);
    } while (f);
    return std::numeric_limits<int64_t>::max(); // unbounded
  }
  // set basicVar's costs to 0, and then runCore()
  constexpr auto run() -> Rational {
#ifndef NDEBUG
    ASSERT(inCanonicalForm);
    assertCanonical();
#endif
    MutPtrVector<index_type> basicVars{getBasicVariables()};
    MutPtrMatrix<value_type> C{getTableau()};
    int64_t f = 1;
    // zero cost of basic variables to put in canonical form
    for (ptrdiff_t c = 0; c < basicVars.size();) {
      int64_t v = basicVars[c++];
      if ((ptrdiff_t(++v) < C.numCol()) && C[0, v])
        f = NormalForm::zeroWithRowOp(C, row(0), row(c), col(v), f);
    }
    return runCore(f);
  }

  // don't touch variables lex > v
  constexpr void rLexCore(ptrdiff_t v) {
    MutPtrMatrix<value_type> C{getTableau()};
    MutPtrVector<index_type> basicVars{getBasicVariables()};
    MutPtrVector<index_type> basicConstraints{getBasicConstraints()};
    invariant(v > 0);
    while (true) {
      // get new entering variable
      Optional<int> enteringVariable = getEnteringVariable(C[0, _(1, v)]);
      if (!enteringVariable) break;
      auto ev = *enteringVariable;
      auto leaveOpt = getLeavingVariable(C, ev);
      if (!leaveOpt) break;
      auto lVar = ptrdiff_t(*leaveOpt);
      ptrdiff_t leavingVariable = lVar++;
      for (ptrdiff_t i = 0; i < C.numRow(); ++i)
        if (i != ptrdiff_t(lVar))
          NormalForm::zeroWithRowOp(C, row(i), row(lVar), ++col(ev), 0);
      // update basic vars and constraints
      int64_t oldBasicVar = basicVars[leavingVariable];
      basicVars[leavingVariable] = ev;
      if (ptrdiff_t(oldBasicVar) < basicConstraints.size())
        basicConstraints[oldBasicVar] = -1;
      basicConstraints[ev] = index_type(leavingVariable);
    }
  }
  // Assumes all >v have already been lex-minimized
  // v starts at numVars-1
  // returns `false` if `0`, `true` if not zero
  // minimize v, not touching any variable lex > v
  constexpr auto rLexMin(ptrdiff_t v) -> bool {
#ifndef NDEBUG
    ASSERT(inCanonicalForm);
#endif
    MutPtrMatrix<value_type> C{getTableau()};
    MutPtrVector<index_type> basicConstraints{getBasicConstraints()};
    int64_t c = basicConstraints[v];
    if (c < 0) return false;
    if (v == 0) return true;
    // we try to zero `v` or at least minimize it.
    // set cost to 1, and then try to alkalize
    // set v and all > v to 0
    C[0, _(0, 1 + v)] << -C[++c, _(0, 1 + v)];
    C[0, _(1 + v, end)] << 0;
    rLexCore(v);
    return makeZeroBasic(v);
  }
  /// makeZeroBasic(ptrdiff_t v) -> bool
  /// Tries to make `v` non-basic if `v` is zero.
  /// Returns `false` if `v` is zero, `true` otherwise
  constexpr auto makeZeroBasic(ptrdiff_t v) -> bool {
    MutPtrMatrix<value_type> C{getTableau()};
    MutPtrVector<index_type> basicVars{getBasicVariables()};
    MutPtrVector<index_type> basicConstraints{getBasicConstraints()};
    int64_t c = basicConstraints[v];
    int64_t cc = c++;
    // was not basic
    // not basic, v is  zero
    if (cc < 0) return false;
    // v is basic, but not zero
    if (C[c, 0] != 0) return true;
#ifndef NDEBUG
    assertCanonical();
#endif
    // so v is basic and zero.
    // We're going to try to make it non-basic
    for (ptrdiff_t ev = 0; ev < v;) {
      auto evm1 = ev++;
      if ((basicConstraints[evm1] >= 0) || (C[c, ev] == 0)) continue;
      if (C[c, ev] < 0) C[c, _] *= -1;
      for (ptrdiff_t i = 1; i < C.numRow(); ++i)
        if (i != ptrdiff_t(c))
          NormalForm::zeroWithRowOp(C, row(i), row(c), col(ev), 0);
      int64_t oldBasicVar = basicVars[cc];
      invariant(oldBasicVar == int64_t(v));
      basicVars[cc] = index_type(evm1);
      // if (ptrdiff_t(oldBasicVar) < basicConstraints.size())
      basicConstraints[oldBasicVar] = -1;
      basicConstraints[evm1] = index_type(cc);
      break;
    }
#ifndef NDEBUG
    assertCanonical();
#endif
    return false;
  }
  constexpr auto rLexMinLast(ptrdiff_t n) -> Solution {
#ifndef NDEBUG
    ASSERT(inCanonicalForm);
    assertCanonical();
#endif
    for (ptrdiff_t v = getNumVars(), e = v - n; v != e;) rLexMin(--v);
#ifndef NDEBUG
    assertCanonical();
#endif
    return {this, length(getNumVars() - n), length(getNumVars())};
  }
  constexpr auto rLexMinStop(ptrdiff_t skippedVars) -> Solution {
#ifndef NDEBUG
    ASSERT(inCanonicalForm);
    assertCanonical();
#endif
    for (ptrdiff_t v = getNumVars(); v != skippedVars;) rLexMin(--v);
#ifndef NDEBUG
    assertCanonical();
#endif
    return {this, length(skippedVars), length(getNumVars())};
  }

  // reverse lexicographic ally minimize vars
  constexpr void rLexMin(Vector<Rational> &sol) {
    sol << rLexMinLast(sol.size());
  }
  // A(:,1:end)*x <= A(:,0)
  // B(:,1:end)*x == B(:,0)
  // returns a Simplex if feasible, and an empty `Optional` otherwise
  static constexpr auto
  positiveVariables(Arena<> *alloc, PtrMatrix<int64_t> A,
                    PtrMatrix<int64_t> B) -> Optional<Simplex *> {
    invariant(A.numCol() == B.numCol());
    ptrdiff_t numVar = ptrdiff_t(A.numCol()) - 1,
              numSlack = ptrdiff_t(A.numRow()),
              numStrict = ptrdiff_t(B.numRow()), numCon = numSlack + numStrict,
              varCap = numVar + numSlack;
    // see how many slack vars are infeasible as solution
    // each of these will require an augment variable
    for (ptrdiff_t i = 0; i < numSlack; ++i) varCap += A[i, 0] < 0;
    // try to avoid reallocating
    auto checkpoint{alloc->checkpoint()};
    Simplex *simplex{Simplex::create(alloc, row(numCon), col(numVar + numSlack),
                                     capacity(numCon), rowStride(varCap))};
    // construct:
    // [ I A
    //   0 B ]
    // then drop the extra variables
    slackEqualityConstraints(simplex->getConstraints()[_, _(1, end)],
                             A[_, _(1, end)], B[_, _(1, end)]);
    auto consts{simplex->getConstants()};
    consts[_(0, numSlack)] << A[_, 0];
    if (numStrict) consts[_(numSlack, numSlack + numStrict)] << B[_, 0];
    // for (ptrdiff_t i = 0; i < numSlack; ++i) consts[i] = A(i, 0);
    // for (ptrdiff_t i = 0; i < numStrict; ++i) consts[i + numSlack] = B(i, 0);
    if (!simplex->initiateFeasible()) return simplex;
    alloc->rollback(checkpoint);
    return nullptr;
  }
  static constexpr auto positiveVariables(Arena<> *alloc, PtrMatrix<int64_t> A)
    -> Optional<Simplex *> {
    ptrdiff_t numVar = ptrdiff_t(A.numCol()) - 1,
              numSlack = ptrdiff_t(A.numRow()), numCon = numSlack,
              varCap = numVar + numSlack;
    // see how many slack vars are infeasible as solution
    // each of these will require an augment variable
    for (ptrdiff_t i = 0; i < numSlack; ++i) varCap += A[i, 0] < 0;
    // try to avoid reallocating
    auto checkpoint{alloc->checkpoint()};
    Simplex *simplex{Simplex::create(alloc, row(numCon), col(numVar + numSlack),
                                     capacity(numCon), rowStride(varCap))};
    // construct:
    // [ I A ]
    // then drop the extra variables
    slackEqualityConstraints(simplex->getConstraints()[_, _(1, end)],
                             A[_, _(1, end)]);
    // auto consts{simplex.getConstants()};
    // for (ptrdiff_t i = 0; i < numSlack; ++i) consts[i] = A(i, 0);
    simplex->getConstants() << A[_, 0];
    if (!simplex->initiateFeasible()) return simplex;
    alloc->rollback(checkpoint);
    return nullptr;
  }

  constexpr void pruneBounds(Arena<> *alloc, ptrdiff_t numSlack = 0) {
    auto p = alloc->scope();
    Simplex *simplex{Simplex::create(alloc, numConstraints, numVars,
                                     constraintCapacity,
                                     --auto{varCapacityP1})};
    // Simplex simplex{getNumCons(), getNumVars(), getNumSlack(), 0};
    for (ptrdiff_t c = 0; c < getNumCons(); ++c) {
      *simplex << *this;
      MutPtrMatrix<int64_t> constraints = simplex->getConstraints();
      int64_t bumpedBound = ++constraints[c, 0];
      MutPtrVector<int64_t> cost = simplex->getCost();
      for (ptrdiff_t v = numSlack; v < cost.size(); ++v)
        cost[v] = -constraints[c, v + 1];
      if (simplex->run() != bumpedBound) deleteConstraint(c--);
    }
  }

  constexpr void dropVariable(ptrdiff_t i) {
    // We remove a variable by isolating it, and then dropping the
    // constraint. This allows us to preserve canonical form
    MutPtrVector<index_type> basicConstraints{getBasicConstraints()};
    MutPtrMatrix<value_type> C{getConstraints()};
    // ensure sure `i` is basic
    if (basicConstraints[i] < 0) makeBasic(C, 0, index_type(i));
    ptrdiff_t ind = basicConstraints[i];
    ptrdiff_t lastRow = ptrdiff_t(C.numRow()) - 1;
    if (lastRow != ind) swap(C, row(ind), row(lastRow));
    truncateConstraints(lastRow);
  }
  constexpr void removeExtraVariables(ptrdiff_t i) {
    for (ptrdiff_t j = getNumVars(); j > i;) {
      dropVariable(--j);
      truncateVars(j);
    }
  }
  // static constexpr auto toMask(PtrVector<int64_t> x) -> uint64_t {
  //   assert(x.size() <= 64);
  //   uint64_t m = 0;
  //   for (auto y : x) m = ((m << 1) | (y != 0));
  //   return m;
  // }
  // [[nodiscard]] constexpr auto getBasicTrueVarMask() const -> uint64_t {
  //   const ptrdiff_t numVarTotal = getNumVars();
  //   assert(numVarTotal <= 64);
  //   uint64_t m = 0;
  //   PtrVector<index_type> basicCons{getBasicConstraints()};
  //   for (ptrdiff_t i = numSlack; i < numVarTotal; ++i)
  //     m = ((m << 1) | (basicCons[i] > 0));
  //   return m;
  // }
  // check if a solution exists such that `x` can be true.
  // returns `true` if unsatisfiable
  [[nodiscard]] constexpr auto unSatisfiable(Arena<> alloc,
                                             PtrVector<int64_t> x,
                                             ptrdiff_t off) const -> bool {
    // is it a valid solution to set the first `x.size()` variables to
    // `x`? first, check that >= 0 constraint is satisfied
    if (!allGEZero(x)) return true;
    // approach will be to move `x.size()` variables into the
    // equality constraints, and then check if the remaining sub-problem
    // is satisfiable.
    const ptrdiff_t numCon = getNumCons(), numVar = getNumVars(),
                    numFix = x.size();
    Simplex *subSimp{
      Simplex::create(&alloc, row(numCon), col(numVar - numFix))};
    // subSimp.tableau(0, 0) = 0;
    // subSimp.tableau(0, 1) = 0;
    auto fC{getTableau()};
    auto sC{subSimp->getTableau()};
    sC[_, 0] << fC[_, 0] - fC[_, _(1 + off, 1 + off + numFix)] * x.t();
    // sC(_, 0) = fC(_, 0);
    // for (ptrdiff_t i = 0; i < numFix; ++i)
    //     sC(_, 0) -= x(i) * fC(_, i + 1 + off);
    sC[_, _(1, 1 + off)] << fC[_, _(1, 1 + off)];
    sC[_, _(1 + off, end)] << fC[_, _(1 + off + numFix, end)];
    // returns `true` if unsatisfiable
    return subSimp->initiateFeasible();
  }
  [[nodiscard]] constexpr auto satisfiable(Arena<> alloc, PtrVector<int64_t> x,
                                           ptrdiff_t off) const -> bool {
    return !unSatisfiable(alloc, x, off);
  }
  // check if a solution exists such that `x` can be true.
  // zeros remaining rows
  [[nodiscard]] constexpr auto
  unSatisfiableZeroRem(Arena<> alloc, PtrVector<int64_t> x, ptrdiff_t off,
                       ptrdiff_t numRow) const -> bool {
    // is it a valid solution to set the first `x.size()` variables to
    // `x`? first, check that >= 0 constraint is satisfied
    if (!allGEZero(x)) return true;
    // approach will be to move `x.size()` variables into the
    // equality constraints, and then check if the remaining sub-problem
    // is satisfiable.
    invariant(numRow <= getNumCons());
    const ptrdiff_t numFix = x.size();
    Simplex *subSimp{Simplex::create(&alloc, row(numRow), col(off++))};
    auto fC{getConstraints()};
    auto sC{subSimp->getConstraints()};
    sC[_, 0] << fC[_(begin, numRow), 0] -
                  fC[_(begin, numRow), _(off, off + numFix)] * x.t();
    sC[_, _(1, off)] << fC[_(begin, numRow), _(1, off)];
    return subSimp->initiateFeasible();
  }
  /// indsFree gives how many variables are free to take  any >= 0 value
  /// indOne is var ind greater than indsFree that's pinned to 1
  /// (i.e., indsFree + indOne == index of var pinned to 1)
  /// numRow is number of rows used, extras are dropped
  // [[nodiscard]] constexpr auto
  [[nodiscard]] inline auto
  unSatisfiableZeroRem(Arena<> alloc, ptrdiff_t iFree,
                       std::array<ptrdiff_t, 2> inds,
                       ptrdiff_t numRow) const -> bool {
    invariant(numRow <= getNumCons());
    Simplex *subSimp{Simplex::create(&alloc, row(numRow), col(iFree++))};
    auto fC{getConstraints()};
    auto sC{subSimp->getConstraints()};
    auto r = _(0, numRow);
    sC[_, 0] << fC[r, 0] - (fC[r, inds[0] + iFree] + fC[r, inds[1] + iFree]);
    sC[_, _(1, iFree)] << fC[r, _(1, iFree)];
    return subSimp->initiateFeasible();
  }
  [[nodiscard]] constexpr auto
  satisfiableZeroRem(Arena<> alloc, PtrVector<int64_t> x, ptrdiff_t off,
                     ptrdiff_t numRow) const -> bool {
    return !unSatisfiableZeroRem(alloc, x, off, numRow);
  }
  void printResult(ptrdiff_t numSlack = 0) {
    auto C{getConstraints()};
    auto basicVars{getBasicVariables()};
    for (ptrdiff_t i = 0; i < basicVars.size(); ++i) {
      ptrdiff_t v = basicVars[i];
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
  static constexpr auto create(Arena<> *alloc, Row<> numCon,
                               Col<> numVar) -> Valid<Simplex> {
    return create(alloc, numCon, numVar, capacity(ptrdiff_t(numCon)),
                  rowStride(ptrdiff_t(numVar) + ptrdiff_t(numCon)));
  }
  static constexpr auto create(Arena<> *alloc, Row<> numCon, Col<> numVar,
                               Capacity<> conCap,
                               RowStride<> varCap) -> Valid<Simplex> {
    varCap = alignVarCapacity(varCap);
    auto cCap = ptrdiff_t(conCap), vCap = ptrdiff_t(varCap);
    size_t memNeeded = requiredMemory(cCap, vCap);
    auto *mem = (Simplex *)alloc->allocate<alignof(Simplex)>(memNeeded);
    mem->numConstraints = numCon;
    mem->numVars = numVar;
    mem->constraintCapacity = conCap;
    mem->varCapacityP1 = varCap;
    return mem;
  }

  static auto operator new(size_t count, Capacity<> conCap,
                           RowStride<> varCap) -> void * {
    auto cC = ptrdiff_t(conCap), vC = ptrdiff_t(varCap);
    size_t memNeeded = requiredMemory(cC, vC);
    // void *p = ::operator new(count * memNeeded);
    return poly::alloc::malloc(count * memNeeded,
                               std::align_val_t(alignof(Simplex)));
    // return ::operator new(count * memNeeded,
    //                       std::align_val_t(alignof(Simplex)));
  }
  static void operator delete(void *ptr, size_t sz) {
    poly::alloc::free(ptr, sz, std::align_val_t(alignof(Simplex)));
    // ::operator delete(ptr, std::align_val_t(alignof(Simplex)));
  }

  static auto create(Row<> numCon, Col<> numVar) -> std::unique_ptr<Simplex> {
    auto nc = ptrdiff_t(numCon);
    return create(numCon, numVar, capacity(nc),
                  rowStride(ptrdiff_t(numVar) + nc));
  }
  static auto create(Row<> numCon, Col<> numVar, Capacity<> conCap,
                     RowStride<> varCap) -> std::unique_ptr<Simplex> {
    varCap = alignVarCapacity(varCap);
    auto *ret = new (conCap, varCap) Simplex;
    ret->numConstraints = numCon;
    ret->numVars = numVar;
    ret->constraintCapacity = conCap;
    ret->varCapacityP1 = varCap;
    return std::unique_ptr<Simplex>(ret);
  }

  static constexpr auto
  create(Arena<> *alloc, Row<> numCon,
         Col<> numVar, // NOLINT(bugprone-easily-swappable-parameters)
         ptrdiff_t numSlack) -> Valid<Simplex> {
    ptrdiff_t conCap = ptrdiff_t(numCon),
              varCap = ptrdiff_t(numVar) + numSlack + conCap;
    return create(alloc, numCon, numVar, capacity(conCap), rowStride(varCap));
  }
  constexpr auto copy(Arena<> *alloc) const -> Valid<Simplex> {
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
  friend inline auto operator<<(std::ostream &os,
                                const Simplex &s) -> std::ostream & {
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
static_assert(AbstractVector<ElementwiseBinaryOp<
                PtrVector<Rational>, PtrVector<Rational>, std::minus<>>>);
static_assert(std::movable<Simplex::Solution::iterator>);
static_assert(std::indirectly_readable<Simplex::Solution::iterator>);
static_assert(std::forward_iterator<Simplex::Solution::iterator>);
static_assert(alignof(Simplex) == simd::VECTORWIDTH);
} // namespace poly::math

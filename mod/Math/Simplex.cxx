module;
#include "Macros.hxx"
module Simplex;
import Allocator;

namespace math {

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

// Private module fragment - static utility functions
namespace {
struct NoFilter {
  static constexpr auto operator()(std::int64_t) -> bool { return false; }
};

// 1 based to match getBasicConstraints
[[nodiscard]] TRIVIAL constexpr auto
getEnteringVariable(PtrVector<std::int64_t> costs)
  -> std::optional<Simplex::index_t> {
  // Bland's algorithm; guaranteed to terminate
  std::ptrdiff_t idx = find_first(costs, [](auto c) { return c < 0; });
  return idx >= 0 ? idx : std::optional<Simplex::index_t>{};
}

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
[[nodiscard]] TRIVIAL constexpr auto
getLeavingVariable(PtrMatrix<std::int64_t> C, std::ptrdiff_t entering_variable,
                   const auto &filter) -> std::optional<Simplex::index_t> {
  // inits guarantee first valid is selected
  std::int64_t dj = -1, nj = 0;
  Simplex::index_t j = 0;
  StridedVector<std::int64_t> cv{C[_, entering_variable + 1]};
  for (Simplex::index_t i = 1; i < C.numRow(); ++i) {
    if constexpr (!std::same_as<std::remove_cvref_t<decltype(filter)>,
                                NoFilter>)
      if (filter(i)) continue;
    std::int64_t di = cv[i];
    if (di <= 0) continue;
    std::int64_t ni = C[i, 0];
    if (!ni) return --i;
    invariant(ni > 0);
    if ((dj * ni) >= (di * nj)) continue;
    dj = di;
    nj = ni;
    j = i;
  }
  // NOTE: if we fail to find a leaving variable, then `j = 0`,
  // and it will unsigned wrap to `std::ptrdiff_t(-1)`, which indicates
  // an empty `Optional<unsigned int>`
  return j ? --j : std::optional<Simplex::index_t>{};
}

[[nodiscard]] TRIVIAL constexpr auto
getLeavingVariable(PtrMatrix<std::int64_t> C, std::ptrdiff_t entering_variable)
  -> std::optional<Simplex::index_t> {
  return getLeavingVariable(C, entering_variable, NoFilter{});
}
// (varCapacity+1)%simd::Width<std::int64_t>==0
TRIVIAL constexpr auto alignVarCapacity(RowStride<> rs) -> RowStride<> {
  std::ptrdiff_t r = std::ptrdiff_t(rs) + 1;
  static constexpr std::ptrdiff_t W = simd::Width<std::int64_t>;
  if constexpr (W <= 1) return stride(r);
  else return stride(__builtin_align_up(r, W));
}
[[nodiscard]] TRIVIAL constexpr auto reservedTableau(std::ptrdiff_t cons,
                                                     std::ptrdiff_t vars)
  -> std::ptrdiff_t {
  return static_cast<std::ptrdiff_t>(sizeof(Simplex::value_type)) *
         ((cons + 1) * vars);
}
template <std::ptrdiff_t N>
TRIVIAL constexpr auto align(std::ptrdiff_t x) -> std::ptrdiff_t {
  return (x + (N - 1z)) & (-N);
}
template <std::integral T>
TRIVIAL constexpr auto alignOffset(std::ptrdiff_t x) -> std::ptrdiff_t {
  static constexpr std::ptrdiff_t W =
    std::ptrdiff_t(std::max(simd::VECTORWIDTH / sizeof(T), alignof(T)));
  return align<W>(x);
}
// offset in bytes
TRIVIAL constexpr auto tableauOffset(std::ptrdiff_t cons, std::ptrdiff_t vars)
  -> std::ptrdiff_t {
  std::ptrdiff_t coff = alignOffset<Simplex::index_t>(cons);
  std::ptrdiff_t voff = alignOffset<Simplex::index_t>(vars);
  std::ptrdiff_t offset = coff + voff;
  // std::ptrdiff_t offset =
  //   alignOffset<index_t>(cons) + alignOffset<index_t>(vars);
  return static_cast<std::ptrdiff_t>(sizeof(Simplex::index_t)) * offset;
}
TRIVIAL constexpr auto requiredMemory(std::ptrdiff_t cons, std::ptrdiff_t vars)
  -> std::size_t {
  std::ptrdiff_t base = static_cast<std::ptrdiff_t>(sizeof(Simplex)),
                 indices = tableauOffset(cons, vars),
                 tableau = reservedTableau(cons, vars);
  static constexpr std::ptrdiff_t A = std::ptrdiff_t(
    std::max(simd::VECTORWIDTH / sizeof(std::int64_t), alignof(Simplex)));
  return static_cast<std::size_t>(align<A>(base + indices + tableau));
}
auto removeAugmentVars(
  Simplex *S,
  const containers::BitSet<math::Vector<std::uint64_t, 8>> &augmentVars)
  -> bool {
  std::ptrdiff_t num_augment = augmentVars.size(),
                 old_num_var = std::ptrdiff_t(S->getNumVars());
  S->setNumVars(old_num_var + num_augment);
  MutPtrMatrix<Simplex::value_type> C{S->getConstraints()};
  MutPtrVector<Simplex::index_t> basic_vars{S->getBasicVariables()};
  MutPtrVector<Simplex::index_t> basic_cons{S->getBasicConstraints()};
  MutPtrVector<Simplex::value_type> costs{S->getCost()};
  costs.zero();
  C[_, _(old_num_var + 1, end)] << 0;
  for (std::ptrdiff_t i = 0; std::ptrdiff_t a : augmentVars) {
    basic_vars[a] = Simplex::index_t(i) + Simplex::index_t(old_num_var);
    basic_cons[i + old_num_var] = Simplex::index_t(a);
    C[a, old_num_var + (++i)] = 1;
    // we now zero out the implicit cost of `1`
    costs[_(begin, old_num_var + 1)] -= C[a, _(begin, old_num_var + 1)];
  }
#ifndef NDEBUG
  if (anyLTZero(basic_vars)) __builtin_trap();
#endif
  // false/0 means feasible
  // true/non-zero infeasible
  if (S->runCore()) {
    S->truncateVars(old_num_var);
    return true;
  }
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
        NormalForm::zeroColumn(C, row(c), ++col(v));
        basic_vars[c] = Simplex::index_t(v);
        basic_cons[v] = Simplex::index_t(c);
        break;
      }
    }
  }
  // all augment vars are now 0
  S->truncateVars(old_num_var);
#ifndef NDEBUG
  S->assertCanonical();
#endif
  return false;
}
} // namespace

// Out-of-line Simplex member function definitions

[[gnu::returns_nonnull, nodiscard]] auto Simplex::basicConsPointer() const
  -> Simplex::index_t * {
  void *p = const_cast<char *>(memory_);
  return std::assume_aligned<simd::VECTORWIDTH>(
    static_cast<Simplex::index_t *>(p));
}
[[gnu::returns_nonnull, nodiscard]] auto Simplex::basicVarsPointer() const
  -> Simplex::index_t * {
  std::ptrdiff_t offset =
    alignOffset<Simplex::index_t>(reservedBasicConstraints());
  return std::assume_aligned<simd::VECTORWIDTH>(basicConsPointer() + offset);
}
[[gnu::returns_nonnull, nodiscard]] auto Simplex::tableauPointer() const
  -> value_type * {
  std::ptrdiff_t offset = tableauOffset(Simplex::reservedBasicConstraints(),
                                        Simplex::reservedBasicVariables());
  void *p = const_cast<char *>(memory_) + offset;
  return std::assume_aligned<simd::VECTORWIDTH>(static_cast<value_type *>(p));
}

// tableau is constraint * var matrix w/ extra col for LHS
// and extra row for objective function
[[nodiscard]] auto Simplex::reservedBasicConstraints() const -> std::ptrdiff_t {
  return std::ptrdiff_t(var_capacity_p1_) - 1;
}
[[nodiscard]] auto Simplex::reservedBasicVariables() const -> std::ptrdiff_t {
  return std::ptrdiff_t(constraint_capacity_);
}

//  [[nodiscard]]  auto Simplex::intsNeeded() const ->
// std::ptrdiff_t {
//   return reservedTableau() + reservedBasicConstraints() +
//          reservedBasicVariables();
// }
// We with to align every row of `getConstraints()` and
// `getTableau()`. To do this, we 0. (varCapacity+1) %
// simd::Width<std::int64_t> == 0
// 1. offset...
//
// tableau has a stride of `varCapacityP1`, which is the maximum var
// capacity +1. The `+1` is room for the LHS.
//
/// [ value | objective function ]
/// [ LHS   | tableau            ]
[[nodiscard]] auto Simplex::getTableau() const -> PtrMatrix<value_type> {
  //
  return {tableauPointer(), StridedDims<>{
                              ++auto(num_constraints_),
                              ++auto(num_vars_),
                              var_capacity_p1_,
                            }};
}
// NOLINTNEXTLINE(readability-make-member-function-const)
[[nodiscard]] auto Simplex::getTableau() -> MutPtrMatrix<value_type> {
  return {tableauPointer(), StridedDims<>{
                              ++auto(num_constraints_),
                              ++auto(num_vars_),
                              var_capacity_p1_,
                            }};
}
[[nodiscard]] auto Simplex::getConstraints() const -> PtrMatrix<value_type> {
  return {tableauPointer() + std::ptrdiff_t(var_capacity_p1_),
          StridedDims<>{
            num_constraints_,
            ++auto(num_vars_),
            var_capacity_p1_,
          }};
}
void Simplex::zeroConstraints() {
  MutPtrVector<value_type>{
    tableauPointer() + std::ptrdiff_t(var_capacity_p1_),
    length(std::ptrdiff_t(num_constraints_) * std::ptrdiff_t(var_capacity_p1_)),
  }
    .zero();
}
/// num constraints x hcat(LHS, vars)
// NOLINTNEXTLINE(readability-make-member-function-const)
[[nodiscard]] auto Simplex::getConstraints() -> MutPtrMatrix<value_type> {
  return {tableauPointer() + std::ptrdiff_t(var_capacity_p1_),
          StridedDims<>{
            num_constraints_,
            ++auto(num_vars_),
            var_capacity_p1_,
          }};
}
[[nodiscard]] auto Simplex::getBasicConstraints() const
  -> PtrVector<Simplex::index_t> {
  return {basicConsPointer(), aslength(num_vars_)};
}
// maps variables to constraints
// if the constraint is < 0, then that means there is no associated basic
// constraint, and the variable's value is `0`, i.e. it is non-basic
// Otherwise, it is the index of the only non-zero constraint for that
// variable.
[[nodiscard]] auto Simplex::getBasicConstraints()
  -> MutPtrVector<Simplex::index_t> {
  return {basicConsPointer(), aslength(num_vars_)};
}
// maps constraints to variables
[[nodiscard]] auto Simplex::getBasicVariables() const
  -> PtrVector<Simplex::index_t> {
  return {basicVarsPointer(), length(std::ptrdiff_t(num_constraints_))};
}
[[nodiscard]] auto Simplex::getBasicVariables()
  -> MutPtrVector<Simplex::index_t> {
  return {basicVarsPointer(), length(std::ptrdiff_t(num_constraints_))};
}
[[nodiscard]] auto Simplex::getCost() const -> PtrVector<value_type> {
  return {tableauPointer(), length(std::ptrdiff_t(num_vars_) + 1z)};
}
// NOLINTNEXTLINE(readability-make-member-function-const)
[[nodiscard]] auto Simplex::getCost() -> MutPtrVector<value_type> {
  return {tableauPointer(), length(std::ptrdiff_t(num_vars_) + 1z)};
}
// maps variables to constraints; <0 if variable is not basic
[[nodiscard]] auto Simplex::getBasicConstraint(std::ptrdiff_t var) const
  -> Simplex::index_t {
  return getBasicConstraints()[var];
}
// get the variable that is basic associated with this constraint;
// must always be valid while the simplex is in canonical form.
[[nodiscard]] auto Simplex::getBasicVariable(std::ptrdiff_t constraint) const
  -> Simplex::index_t {
  return getBasicVariables()[constraint];
}
[[nodiscard]] auto Simplex::getObjectiveCoefficient(std::ptrdiff_t i) const
  -> value_type {
  return getCost()[++i];
}
[[nodiscard]] auto Simplex::getObjectiveValue() -> value_type & {
  return getCost()[0];
}
[[nodiscard]] auto Simplex::getObjectiveValue() const -> value_type {
  return getCost()[0];
}
void Simplex::simplifySystem() {
#ifndef NDEBUG
  in_canonical_form_ = false;
#endif
  MutPtrMatrix<value_type> C{getConstraints()};
#ifndef NDEBUG
  for (std::ptrdiff_t r = 0, R = std::ptrdiff_t(numRows(C)); r < R; ++r)
    for (std::ptrdiff_t c = 0, N = std::ptrdiff_t(numCols(C)); c < N; ++c)
      invariant(C[r, c] != std::numeric_limits<std::int64_t>::min());
#endif
  NormalForm::solveSystemSkip(C);
  truncateConstraints(std::ptrdiff_t(NormalForm::numNonZeroRows(C)));
}
#ifndef NDEBUG
void Simplex::assertCanonical() const {
  PtrMatrix<value_type> C{getTableau()};
  PtrVector<Simplex::index_t> basic_vars{getBasicVariables()};
  PtrVector<Simplex::index_t> basic_cons{getBasicConstraints()};
  for (std::ptrdiff_t v = 0; v < basic_cons.size();) {
    Simplex::index_t c = basic_cons[v++];
    if (c < 0) continue;
    if (!allZero(C[_(1, 1 + c), v])) __builtin_trap();
    if (!allZero(C[_(2 + c, end), v])) __builtin_trap();
    if (std::ptrdiff_t(basic_vars[c]) != (v - 1)) __builtin_trap();
  }
  for (std::ptrdiff_t c = 1; c < C.numRow(); ++c) {
    Simplex::index_t v = basic_vars[c - 1];
    if (std::ptrdiff_t(v) < basic_cons.size()) {
      invariant(c - 1, std::ptrdiff_t(basic_cons[v]));
      invariant(C[c, v + 1] >= 0);
    }
    invariant(C[c, 0] >= 0);
  }
}
#endif
[[nodiscard]] auto Simplex::getConstants() -> MutStridedVector<std::int64_t> {
  return getTableau()[_(1, end), 0];
}
[[nodiscard]] auto Simplex::getConstants() const
  -> StridedVector<std::int64_t> {
  return getTableau()[_(1, end), 0];
}
void Simplex::truncateConstraints(std::ptrdiff_t i) {
  invariant(std::ptrdiff_t(num_constraints_) <=
            std::ptrdiff_t(constraint_capacity_));
  invariant(i >= 0z);
  invariant(i <= num_constraints_);
  num_constraints_ = row(i);
}
void Simplex::deleteConstraint(std::ptrdiff_t c) {
  auto basic_cons = getBasicConstraints();
  auto basic_vars = getBasicVariables();
  auto constraints = getConstraints();
  --num_constraints_;
  if (auto basic_var = basic_vars[c]; basic_var >= 0)
    basic_cons[basic_var] = -1;
  if (c == num_constraints_) return;
  auto basic_var = basic_vars[std::ptrdiff_t(num_constraints_)];
  basic_vars[c] = basic_var;
  if (basic_var >= 0) basic_cons[basic_var] = Simplex::index_t(c);
  constraints[c, _] << constraints[num_constraints_, _];
}

[[nodiscard]] auto Simplex::getSolution() const -> Simplex::Solution {
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
[[nodiscard("returns `true` if infeasible; should check when calling.")]] auto
Simplex::initiateFeasible() -> bool {
  // remove trivially redundant constraints
  simplifySystem();
  // [ I;  X ; b ]
  //
  // original number of variables
  const std::ptrdiff_t num_var = getNumVars();
  MutPtrMatrix<value_type> C{getConstraints()};
  MutPtrVector<Simplex::index_t> basic_cons{getBasicConstraints()};
  // initialize to `-2`
  basic_cons << -2;
  // first pass, we make sure the equalities are >= 0
  // and we eagerly try and find columns with
  // only a single non-0 element.
  for (std::ptrdiff_t c = 0; c < C.numRow(); ++c) {
    std::int64_t &Ceq = C[c, 0];
    std::int64_t sign = (2 * (Ceq >= 0)) - 1;
    Ceq *= sign;
    // was initialized to -2
    // if unset here (i.e. == -2) and >1, try to make basic
    // if set earlier and we're resetting (or < 0), set to -1
    for (std::ptrdiff_t v = 0; v < num_var; ++v)
      if (std::int64_t Ccv = C[c, v + 1] *= sign)
        basic_cons[v] =
          (((basic_cons[v] == -2) && (Ccv > 0))) ? Simplex::index_t(c) : -1;
  }
  // basicCons should now contain either `-1` or an integer >= 0
  // indicating which row contains the only non-zero element; we'll
  // now fill basicVars.
  //
  MutPtrVector<Simplex::index_t> basic_vars{getBasicVariables()};
  basic_vars << -1;
  for (std::ptrdiff_t v = 0; v < num_var; ++v) {
    if (std::int64_t r = basic_cons[v]; r >= 0) {
      if (basic_vars[r] == -1) basic_vars[r] = Simplex::index_t(v);
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
  return (!aug_vars.empty() && removeAugmentVars(this, aug_vars));
}

auto Simplex::makeBasic(MutPtrMatrix<std::int64_t> C,
                        Simplex::index_t enteringVar, std::int64_t f)
  -> std::int64_t {
  std::optional<Simplex::index_t> leave_opt =
    getLeavingVariable(C, enteringVar);
  if (!leave_opt) return 0; // unbounded
  Simplex::index_t leaving_var = *leave_opt;
  Row<> pivot_row = ++row(leaving_var);
  Col<> pivot_col = ++col(enteringVar);
  f = NormalForm::zeroColumn(C, pivot_row, pivot_col, f);
  // update basic vars and constraints
  MutPtrVector<Simplex::index_t> basic_vars{getBasicVariables()},
    basic_constraints{getBasicConstraints()};
  Simplex::index_t old_basic_var = basic_vars[leaving_var];
  basic_vars[leaving_var] = enteringVar;
  basic_constraints[old_basic_var] = -1;
  basic_constraints[enteringVar] = leaving_var;
  return f;
}
// returns `true` if it failed to make basic
auto Simplex::makeBasic(MutPtrMatrix<std::int64_t> C, Simplex::index_t ev)
  -> bool {
  std::optional<Simplex::index_t> leave_opt = getLeavingVariable(C, ev);
  if (leave_opt) makeBasic(C, ev, *leave_opt);
  return !leave_opt;
}
void Simplex::makeBasic(MutPtrMatrix<std::int64_t> C, Simplex::index_t ev,
                        Simplex::index_t l_var) {
  Simplex::index_t leaving_variable = l_var++;
  Row<> pivot_row = row(l_var);
  Col<> pivot_col = ++col(ev);
  NormalForm::zeroColumn(C, pivot_row, pivot_col);
  // update basic vars and constraints
  MutPtrVector<Simplex::index_t> basic_vars{getBasicVariables()},
    basic_constraints{getBasicConstraints()};
  Simplex::index_t old_basic_var = basic_vars[leaving_variable];
  basic_vars[leaving_variable] = ev;
  if (std::ptrdiff_t(old_basic_var) < basic_constraints.size())
    basic_constraints[old_basic_var] = -1;
  basic_constraints[ev] = leaving_variable;
}
// run the simplex algorithm, assuming basicVar's costs have been set to
// 0
auto Simplex::runCore(std::int64_t f) -> Rational {
#ifndef NDEBUG
  if (!in_canonical_form_) __builtin_trap();
#endif
  MutPtrMatrix<std::int64_t> C{getTableau()};
  do {
    // entering variable is the column
    std::optional<Simplex::index_t> entering_variable =
      getEnteringVariable(C[0, _(1, end)]);
    if (!entering_variable) return Rational::create(C[0, 0], f);
    f = makeBasic(C, *entering_variable, f);
  } while (f);
  return std::numeric_limits<std::int64_t>::max(); // unbounded
}
/// Set basicVar's costs to 0, and then runCore()
/// Essentially, minimize the `getCost()` expression
auto Simplex::run() -> Rational {
#ifndef NDEBUG
  if (!in_canonical_form_) __builtin_trap();
  assertCanonical();
#endif
  MutPtrVector<Simplex::index_t> basic_vars{getBasicVariables()};
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
auto Simplex::makeBasic(Simplex::index_t var) -> Simplex::index_t {
  MutPtrMatrix<value_type> C{getTableau()};
  if (Simplex::index_t j = getBasicConstraint(var); j >= 0) return j;
  if (makeBasic(C, var)) return getBasicConstraint(var);
  return -1;
}

// don't touch variables lex > v
void Simplex::rLexCore(PtrMatrix<std::int64_t> C, PtrVector<std::int64_t> costs,
                       std::ptrdiff_t v) {
  invariant(v > 0);
  while (true) {
    // get new entering variable
    std::optional<Simplex::index_t> entering_variable =
      getEnteringVariable(costs);
    // we break when no costs < 0
    if (!entering_variable) break;
    Simplex::index_t ev = *entering_variable;
    std::optional<Simplex::index_t> leave_opt = getLeavingVariable(C, ev);
    // or when no constraints were found
    if (!leave_opt) break;
    makeBasic(C, ev, *leave_opt);
  }
}
// Assumes all >v have already been lex-minimized
// v starts at numVars-1
// returns `false` if `0`, `true` if not zero
// minimize v, not touching any variable lex > v
auto Simplex::rLexMin(std::ptrdiff_t v) -> bool {
#ifndef NDEBUG
  if (!in_canonical_form_) __builtin_trap();
#endif
  Simplex::index_t c = getBasicConstraint(v);
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
[[nodiscard]] auto Simplex::getVarValue(std::ptrdiff_t var) const -> Rational {
  std::int64_t j = getBasicConstraint(var);
  if (j < 0) return 0;
  PtrMatrix<std::int64_t> constraints = getConstraints();
  return Rational::create(constraints[j, 0], constraints[j, var + 1]);
}
/// makeZeroBasic(std::ptrdiff_t v) -> bool
/// Tries to make `v` non-basic if `v` is zero.
/// Returns `true` if zero but we couldn't eliminate.
auto Simplex::makeZeroNonBasic(std::ptrdiff_t v) -> bool {
  utils::assume(v > 0);
  MutPtrMatrix<value_type> C{getTableau()};
  MutPtrVector<Simplex::index_t> basic_vars{getBasicVariables()};
  MutPtrVector<Simplex::index_t> basic_constraints{getBasicConstraints()};
  std::int64_t c = basic_constraints[v];
  if (c < 0) return false; // already not basic, hence v is zero
  std::int64_t cc = c++;
  // was not basic
  // v is basic, but not zero
  if (C[c, 0] != 0) return false;
#ifndef NDEBUG
  assertCanonical();
#endif
  // so v is basic and zero.
  // We're going to try to make it non-basic
  // We scan for a variable, trying to find a non-basic one
  std::ptrdiff_t ev = 0;
  for (; (basic_constraints[ev] >= 0) || (C[c, ev + 1] == 0);)
    if (v <= ++ev) return true;
  // we found `ev`

  // we have the invariant that C[c,0] >= 0, but multiplying
  // by `-1` is fine because `C[c,0] == 0`, guaranteed above
  std::ptrdiff_t evp1 = ev + 1;
  if (C[c, evp1] < 0) C[c, _] *= -1;
  // Subset matrix to exclude row 0, adjust pivot_row accordingly
  NormalForm::zeroColumn(C[_(1, end), _], row(c - 1), col(evp1));
  std::int64_t old_basic_var = basic_vars[cc];
  invariant(old_basic_var == std::int64_t(v));
  basic_vars[cc] = Simplex::index_t(ev);
  basic_constraints[old_basic_var] = -1;
  basic_constraints[ev] = Simplex::index_t(cc);
#ifndef NDEBUG
  assertCanonical();
#endif
  return false;
}
auto Simplex::rLexMinLast(std::ptrdiff_t n) -> Simplex::Solution {
#ifndef NDEBUG
  if (!in_canonical_form_) __builtin_trap();
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
auto Simplex::rLexMinStop(std::ptrdiff_t skip_first, std::ptrdiff_t skip_last)
  -> Simplex::Solution {
#ifndef NDEBUG
  if (!in_canonical_form_) __builtin_trap();
  assertCanonical();
#endif
  std::ptrdiff_t num_v = getNumVars() - skip_last;
  invariant(num_v > skip_first);
  for (std::ptrdiff_t v = num_v; v != skip_first;) rLexMin(--v);
#ifndef NDEBUG
  assertCanonical();
#endif
  return {.simplex_ = this,
          .skipped_vars_ = length(skip_first),
          .num_vars_ = length(num_v)};
}
auto Simplex::maximizeLast(std::ptrdiff_t num_last) -> Rational {
  utils::invariant(num_last > 0);
  MutPtrMatrix<value_type> C{getTableau()};
  std::ptrdiff_t num_v = std::ptrdiff_t(C.numCol());
  C[0, _(0, num_v - num_last)] << 0;
  C[0, _(num_v - num_last, num_v)] << -1;
  return run();
}
auto Simplex::maximizeLastDrop(std::ptrdiff_t num_last) -> Rational {
  utils::invariant(num_last > 0);
  MutPtrMatrix<value_type> C{getTableau()};
  std::ptrdiff_t num_v = std::ptrdiff_t(C.numCol()) - 1,
                 last_row = std::ptrdiff_t(C.numRow()) - 1;
  Rational f = maximizeLast(num_last);
  MutPtrVector<Simplex::index_t> basic_cons{getBasicConstraints()},
    basic_vars{getBasicVariables()};
  // now we must drop `num_last`
  // if they are basic, we subtract them from the constraints,
  // and drop the associated constraint
  std::ptrdiff_t v = num_v, num_r = num_v - num_last;
  do {
    Simplex::index_t c = basic_cons[--v];
    if (c < 0) continue; // was 0
    C[c + 1, 0] = 0;
    if (!makeZeroNonBasic(v)) continue;
    C[c + 1, _(v)] << C[last_row--, _(v)];
    Simplex::index_t last_row_var = basic_vars[last_row];
    utils::invariant(basic_vars[c] == v);
    basic_vars[c] = last_row_var;
    basic_cons[last_row_var] = c;
  } while (v != num_r);
  truncateVars(num_r);
  truncateConstraints(last_row);
  return f;
}

// reverse lexicographic ally minimize vars
void Simplex::rLexMin(Vector<Rational> &sol) { sol << rLexMinLast(sol.size()); }
// A(:,1:end)*x <= A(:,0)
// B(:,1:end)*x == B(:,0)
// returns a Simplex if feasible, and an empty `Optional` otherwise
auto Simplex::positiveVariables(alloc::Arena<> *alloc,
                                PtrMatrix<std::int64_t> A,
                                PtrMatrix<std::int64_t> B)
  -> Optional<Simplex &> {
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
  Simplex &simplex{Simplex::create(alloc, row(num_con),
                                   col(num_var + num_slack), capacity(num_con),
                                   stride(var_cap))};
  // construct:
  // [ I A
  //   0 B ]
  // then drop the extra variables
  slackEqualityConstraints(simplex.getConstraints()[_, _(1, end)],
                           A[_, _(1, end)], B[_, _(1, end)]);
  auto consts{simplex.getConstants()};
  consts[_(0, num_slack)] << A[_, 0];
  if (num_strict) consts[_(num_slack, num_slack + num_strict)] << B[_, 0];
  // for (std::ptrdiff_t i = 0; i < numSlack; ++i) consts[i] = A(i, 0);
  // for (std::ptrdiff_t i = 0; i < numStrict; ++i) consts[i + numSlack] =
  // B(i, 0);
  if (!simplex.initiateFeasible()) return simplex;
  alloc->rollback(checkpoint);
  return {};
}
auto Simplex::positiveVariables(alloc::Arena<> *alloc,
                                PtrMatrix<std::int64_t> A)
  -> Optional<Simplex &> {
  std::ptrdiff_t num_var = std::ptrdiff_t(A.numCol()) - 1,
                 num_slack = std::ptrdiff_t(A.numRow()), num_con = num_slack,
                 var_cap = num_var + num_slack;
  // see how many slack vars are infeasible as solution
  // each of these will require an augment variable
  for (std::ptrdiff_t i = 0; i < num_slack; ++i) var_cap += A[i, 0] < 0;
  // try to avoid reallocating
  auto checkpoint{alloc->checkpoint()};
  Simplex &simplex{Simplex::create(alloc, row(num_con),
                                   col(num_var + num_slack), capacity(num_con),
                                   stride(var_cap))};
  // construct:
  // [ I A ]
  // then drop the extra variables
  slackEqualityConstraints(simplex.getConstraints()[_, _(1, end)],
                           A[_, _(1, end)]);
  // auto consts{simplex.getConstants()};
  // for (std::ptrdiff_t i = 0; i < numSlack; ++i) consts[i] = A(i, 0);
  simplex.getConstants() << A[_, 0];
  if (!simplex.initiateFeasible()) return simplex;
  alloc->rollback(checkpoint);
  return {};
}

void Simplex::pruneBounds(alloc::Arena<> *alloc, std::ptrdiff_t numSlack) {
  auto p = alloc->scope();
  Simplex &simplex{Simplex::create(alloc, num_constraints_, num_vars_,
                                   constraint_capacity_,
                                   --auto{var_capacity_p1_})};
  // Simplex simplex{getNumCons(), getNumVars(), getNumSlack(), 0};
  for (std::ptrdiff_t c = 0; c < getNumCons(); ++c) {
    simplex << *this;
    MutPtrMatrix<std::int64_t> constraints = simplex.getConstraints();
    std::int64_t bumped_bound = ++constraints[c, 0];
    MutPtrVector<std::int64_t> cost = simplex.getCost();
    for (std::ptrdiff_t v = numSlack; v < cost.size(); ++v)
      cost[v] = -constraints[c, v + 1];
    if (simplex.run() != bumped_bound) deleteConstraint(c--);
  }
}

// static  auto toMask(PtrVector<std::int64_t> x) -> std::uint64_t {
//   if(x.size() > 64)__builtin_trap();
//   std::uint64_t m = 0;
//   for (auto y : x) m = ((m << 1) | (y != 0));
//   return m;
// }
// [[nodiscard]]  auto getBasicTrueVarMask() const -> std::uint64_t {
//   const std::ptrdiff_t numVarTotal = getNumVars();
//   if(numVarTotal > 64)__builtin_trap();
//   std::uint64_t m = 0;
//   PtrVector<Simplex::index_t> basicCons{getBasicConstraints()};
//   for (std::ptrdiff_t i = numSlack; i < numVarTotal; ++i)
//     m = ((m << 1) | (basicCons[i] > 0));
//   return m;
// }
// check if a solution exists such that `x` can be true.
// returns `true` if unsatisfiable
[[nodiscard]] auto Simplex::unSatisfiable(alloc::Arena<> alloc,
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
  Simplex &sub_simp{
    Simplex::create(&alloc, row(num_con), col(num_var - num_fix))};
  // subSimp.tableau(0, 0) = 0;
  // subSimp.tableau(0, 1) = 0;
  auto fC{getTableau()};
  auto sC{sub_simp.getTableau()};
  sC[_, 0] << fC[_, 0] - fC[_, _(1 + off, 1 + off + num_fix)] * x.t();
  // sC(_, 0) = fC(_, 0);
  // for (std::ptrdiff_t i = 0; i < numFix; ++i)
  //     sC(_, 0) -= x(i) * fC(_, i + 1 + off);
  sC[_, _(1, 1 + off)] << fC[_, _(1, 1 + off)];
  sC[_, _(1 + off, end)] << fC[_, _(1 + off + num_fix, end)];
  // returns `true` if unsatisfiable
  return sub_simp.initiateFeasible();
}
[[nodiscard]] auto Simplex::satisfiable(alloc::Arena<> alloc,
                                        PtrVector<std::int64_t> x,
                                        std::ptrdiff_t off) const -> bool {
  return !unSatisfiable(alloc, x, off);
}
// check if a solution exists such that `x` can be true.
// zeros remaining rows
[[nodiscard]] auto
Simplex::unSatisfiableZeroRem(alloc::Arena<> alloc, PtrVector<std::int64_t> x,
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
  Simplex &sub_simp{Simplex::create(&alloc, row(numRow), col(off++))};
  auto fC{getConstraints()};
  auto sC{sub_simp.getConstraints()};
  sC[_, 0] << fC[_(begin, numRow), 0] -
                fC[_(begin, numRow), _(off, off + num_fix)] * x.t();
  sC[_, _(1, off)] << fC[_(begin, numRow), _(1, off)];
  return sub_simp.initiateFeasible();
}
/// indsFree gives how many variables are free to take  any >= 0 value
/// indOne is var ind greater than indsFree that's pinned to 1
/// (i.e., indsFree + indOne == index of var pinned to 1)
/// numRow is number of rows used, extras are dropped
// [[nodiscard]]  auto
[[nodiscard]] auto
Simplex::unSatisfiableZeroRem(alloc::Arena<> alloc, std::ptrdiff_t iFree,
                              std::array<std::ptrdiff_t, 2> inds,
                              std::ptrdiff_t numRow) const -> bool {
  invariant(numRow <= getNumCons());
  Simplex &sub_simp{Simplex::create(&alloc, row(numRow), col(iFree++))};
  auto fC{getConstraints()};
  auto sC{sub_simp.getConstraints()};
  auto r = _(0, numRow);
  sC[_, 0] << fC[r, 0] - (fC[r, inds[0] + iFree] + fC[r, inds[1] + iFree]);
  sC[_, _(1, iFree)] << fC[r, _(1, iFree)];
  return sub_simp.initiateFeasible();
}
[[nodiscard]] auto
Simplex::satisfiableZeroRem(alloc::Arena<> alloc, PtrVector<std::int64_t> x,
                            std::ptrdiff_t off, std::ptrdiff_t numRow) const
  -> bool {
  return !unSatisfiableZeroRem(alloc, x, off, numRow);
}
void Simplex::printResult(std::ptrdiff_t numSlack) {
  auto C{getConstraints()};
  auto basic_vars{getBasicVariables()};
  for (std::ptrdiff_t i = 0; i < basic_vars.size(); ++i) {
    std::ptrdiff_t v = basic_vars[i];
    if (v <= numSlack) continue;
    if (C[i, 0]) {
      if (++v < C.numCol()) {
        utils::print("v_", v - numSlack, " = ", C[i, 0], " / ", C[i, v], "\n");
      } else {
        utils::print("v_", v, " = ", C[i, 0], "\n");
        __builtin_trap();
      }
    }
  }
}
auto Simplex::create(alloc::Arena<> *alloc, Row<> num_con, Col<> num_var)
  -> Simplex & {
  return create(alloc, num_con, num_var, capacity(std::ptrdiff_t(num_con)),
                stride(std::ptrdiff_t(num_var) + std::ptrdiff_t(num_con)));
}
auto Simplex::create(alloc::Arena<> *alloc, Row<> num_con, Col<> num_var,
                     Capacity<> con_cap, RowStride<> var_cap) -> Simplex & {
  var_cap = alignVarCapacity(var_cap);
  auto c_cap = std::ptrdiff_t(con_cap), v_cap = std::ptrdiff_t(var_cap);
  std::size_t mem_needed = requiredMemory(c_cap, v_cap);
  auto mem =
    static_cast<Simplex *>(alloc->allocate<alignof(Simplex)>(mem_needed));
  mem->num_constraints_ = num_con;
  mem->num_vars_ = num_var;
  mem->constraint_capacity_ = con_cap;
  mem->var_capacity_p1_ = var_cap;
  return *mem;
}

auto Simplex::malloc(Row<> numCon, Col<> numVar) -> Simplex * {
  auto nc = std::ptrdiff_t(numCon);
  return Simplex::malloc(numCon, numVar, capacity(nc),
                         stride(std::ptrdiff_t(numVar) + nc));
}
auto Simplex::malloc(Row<> numCon, Col<> numVar, Capacity<> conCap,
                     RowStride<> varCap) -> Simplex * {
  varCap = alignVarCapacity(varCap);
  std::ptrdiff_t cC = std::ptrdiff_t(conCap), vC = std::ptrdiff_t(varCap);
  std::size_t mem_needed = requiredMemory(cC, vC);
  Simplex *ret = static_cast<Simplex *>(
    alloc::malloc(mem_needed, std::align_val_t(alignof(Simplex))));
  ret->num_constraints_ = numCon;
  ret->num_vars_ = numVar;
  ret->constraint_capacity_ = conCap;
  ret->var_capacity_p1_ = varCap;
  return ret;
}

auto Simplex::Solution::operator[](std::ptrdiff_t i) const -> Rational {
  invariant(i >= 0);
  return simplex_->getVarValue(i + std::ptrdiff_t(skipped_vars_));
}
[[nodiscard]] auto Simplex::Solution::operator[](OffsetEnd k) const
  -> Rational {
  return simplex_->getVarValue(std::ptrdiff_t(num_vars_) - k.offset_);
}
[[nodiscard]] auto Simplex::Solution::iterator::operator*() const -> Rational {
  return (*sol_)[i_];
}

[[nodiscard]] auto Simplex::Solution::denomLCM() const -> std::int64_t {
  std::int64_t l = 1;
  for (auto r : *this) l = lcm(l, r.denominator_);
  return l;
}
DEBUGUSED void Simplex::Solution::print() const {
  utils::print("Simplex::Solution[");
  bool print_comma = false;
  for (Rational b : *this) {
    if (print_comma) utils::print(", ");
    print_comma = true;
    utils::print(b);
  }
  utils::print("]\n");
  utils::flush();
}

auto Simplex::create(alloc::Arena<> *alloc, Row<> numCon, Col<> numVar,
                     std::ptrdiff_t numSlack) -> Simplex & {
  std::ptrdiff_t con_cap = std::ptrdiff_t(numCon),
                 var_cap = std::ptrdiff_t(numVar) + numSlack + con_cap;
  return create(alloc, numCon, numVar, capacity(con_cap), stride(var_cap));
}
auto Simplex::copy(alloc::Arena<> *alloc) const -> Simplex & {
  Simplex &res = create(alloc, row(getNumCons()), col(getNumVars()),
                        getConCap(), getVarCap());
  res << *this;
#ifndef NDEBUG
  res.in_canonical_form_ = in_canonical_form_;
#endif
  return res;
}
auto Simplex::operator<<(const Simplex &other) -> Simplex & {
  setNumCons(other.getNumCons());
  setNumVars(other.getNumVars());
  getTableau() << other.getTableau();
  getBasicVariables() << other.getBasicVariables();
  getBasicConstraints() << other.getBasicConstraints();
  return *this;
}
void Simplex::print() const {
  utils::print("Basic Variables: ");
  utils::printVector(getBasicVariables().begin(), getBasicVariables().end());
  utils::print("Basic Constraints: ");
  utils::printVector(getBasicConstraints().begin(),
                     getBasicConstraints().end());
  utils::print("Constraints:\n");
  utils::printMatrix(getConstraints().data(),
                     std::ptrdiff_t(getConstraints().numRow()),
                     std::ptrdiff_t(getConstraints().numCol()),
                     std::ptrdiff_t(getConstraints().rowStride()));
}
#ifndef NDEBUG
[[gnu::used]] void Simplex::dump() const { print(); }
#endif

static_assert(AbstractVector<Simplex::Solution>);

static_assert(AbstractVector<PtrVector<Rational>>);
// static_assert(AbstractVector<ElementwiseBinaryOp<
//                 PtrVector<Rational>, PtrVector<Rational>, std::minus<>>>);
static_assert(std::movable<Simplex::Solution::iterator>);
static_assert(std::indirectly_readable<Simplex::Solution::iterator>);
static_assert(std::forward_iterator<Simplex::Solution::iterator>);
static_assert(alignof(Simplex) == simd::VECTORWIDTH);
} // namespace math

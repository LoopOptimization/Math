import Testing;
import Arena;
import Array;
import ArrayConcepts;
import ArrayParse;
import ArrayPrint;
import Comparisons;
import CorePrint;
import ManagedArray;
import MatDim;
import Rational;
import Simplex;
import std;

using namespace ::math;
using utils::operator""_mat;
using namespace testing;

auto simplexFromTableau(alloc::Arena<> *alloc, PtrMatrix<std::int64_t> tableau)
  -> Simplex & {
  std::ptrdiff_t numCon = std::ptrdiff_t(tableau.numRow()) - 1;
  std::ptrdiff_t numVar = std::ptrdiff_t(tableau.numCol()) - 1;
  Simplex &simp{Simplex::create(alloc, row(numCon), col(numVar))};
  static constexpr auto check_invalid = [](auto x) {
    return x == std::numeric_limits<std::int64_t>::min();
  };
  invariant(find_first(tableau, check_invalid) < 0);
  simp.getTableau() << tableau;
  PtrMatrix<Simplex::value_type> C{simp.getConstraints()};
  invariant(find_first(C, check_invalid) < 0);
  return simp;
}

auto main(int argc, const char **argv) -> int {
  cfg<override> = {.filter = argc > 1 ? argv[1] : ""};
  alloc::OwningArena<> managed_alloc;

  "DynsolveDynamicDistance"_test = [&managed_alloc] -> void {
    // clang-format off
    // Problem: mismatched dims require bounding constraints.
    auto t{"[0 0  0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 1 -1 0 -1 0 -1 0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0 -1 0  0  0  0;"
            "0 0  1 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  1 0  0 0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  1 0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  1  0 -1  0 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  1  0 -1 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0 -1 1  0 0  0 -1  0  1 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "1 0  0 0  0 0 -1 1  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0 -1 1  0 0  0 0 -1  0  1  0 0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0 1 -1 0 -1 0 -1 0  0  0  0  0 0 -1  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0 0  1 0  0 0  0 0  0  0  0  0 0  0 -1  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0 0  0 0  1 0  0 0  0  0  0  0 0  0  0 -1  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  1 0  0  0  0  0 0  0  0  0 -1;"
            "0 0  0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0 0  1  0 -1  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0 0  0  1  0 -1 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0 0  0 0 -1 1  0 0  0 -1  0  1 0  0  0  0  0;"
           "-1 0  0 0  0 0  0 0  0  0  0  0 0  0 0  0 0 -1 1  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0 0 -1 1  0 0  0 0 -1  0  1  0 0  0  0  0  0]"_mat};
    // clang-format on
    {
      alloc::Arena<> alloc = managed_alloc;
      Simplex &simp{simplexFromTableau(&alloc, t)};
      expect(fatal(!simp.initiateFeasible()));
      std::ptrdiff_t num_nuisance = 22;
      Simplex::Solution sol = simp.rLexMinStop(num_nuisance);
      // sol.print();
      expect(sol[::math::last] == 1);
      expect(eq(sol.size(), 5));
    }
    {
      alloc::Arena<> alloc = managed_alloc;
      Simplex &S{Simplex::create(&alloc, row(8), col(11))};
      MutPtrMatrix<Simplex::value_type> C{S.getConstraints()};
      C[_, 0] << 0;
      C[_, _(1, 10)] << t[_(2, 10), _(2, 11)];
      C[_, _(10, 12)] << 0;
      C[6, 10] = -1;
      C[6, 11] = 1;
      expect(fatal(!S.initiateFeasible()));
      MutPtrVector<Simplex::value_type> c{S.getCost()};
      c.zero();
      c[1] = -1;
      c[3] = -1;
      c[5] = -1;
      expect(S.run() <= 0);
    }
  };
  "DynsolveAppendEqualMatching"_test = [&managed_alloc] -> void {
    // clang-format off
    // Approach 0: append dims so they match
    // We also muset set the new dims equal to the old, which seems wrong.
    auto t{"[0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 1 -1 0 -1 0 -1 0  0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0  0 -1 0  0  0  0;"
            "0 0  1 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  1 0  0 0  0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  1 0  0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  1  0 -1  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  1  0 -1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
           "-1 0  0 0  0 0  0 0 -1  1  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0 -1 1  0 0  0  0  0 -1  0  1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "1 0  0 0  0 0 -1 1  1 -1  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 0 -1 1  0 0  0 0  0  0 -1  0  1  0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 1 -1 0 -1 0 -1 0  0  0  0  0  0  0 0 -1  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  1 0  0 0  0 0  0  0  0  0  0  0 0  0 -1  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  1 0  0 0  0  0  0  0  0  0 0  0  0 -1  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  1 0  0  0  0  0  0  0 0  0  0  0 -1;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  1  0 -1  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  1  0 -1 0  0  0  0  0;"
            "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  0 0 -1  1  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0 -1 1  0 0  0  0  0 -1  0  1 0  0  0  0  0;"
           "-1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 -1 1  1 -1  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0 -1 1  0 0  0 0  0  0 -1  0  1  0 0  0  0  0  0]"_mat};

           // "0 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  1 0  0  0  0  0  0  0 0  0  0  0 -1;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0  0 0 -1  1  0  0  0  0 0  0  0  0  0;"
           //"-1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 -1 1  1 -1  0  0  0  0 0  0  0  0  0;"

           // "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 0  0  0  0  0  0  0  0 0  0  0  0 1;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 0  0 -1  1  0  0  0  0 0  0  0  0 0;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 1 -1 -1  1  0  0  0  0 0  0  0  0 0;"

           // "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 0  0  0  0  0  0  0  0 0  0  0  0 1;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 0  0 -1  1  0  0  0  0 0  0  0  0 0;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0  0 0  0 0  0 0 1 -1 -1  1  0  0  0  0 0  0  0  0 0;"

           //  "0 0  0 0  0 0  0 0 0  0  0  0  0 0  0 0  0 0  1 0  0  0  0  0  0 0  0  0  0 -1;"
           //  "1 0  0 0  0 0  0 0 0  0  0  0  0 0  0 0  0 0  0 0 -1  0  0  0  0 0  0  0  0  0;"
           // "-1 0  0 0  0 0  0 0 0  0  0  0  0 0  0 0  0 0 -1 1  1  0  0  0  0 0  0  0  0  0;"

           // "-1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0 1  1  0  0  0  0 0  0  0  0 -1;"
           //  "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0 0 -1  0  0  0  0 0  0  0  0  0;"
           // "-1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 -1 1  1  0  0  0  0 0  0  0  0  0;"

           // Run through alg:                                e                               
           // "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 -1  0  0  0  0  0  0 0  0  0  0  1;" costs
           // "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  1  0  0  0  0  0  0 0  0  0  0 -1;" n = 0
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0  0 -1  0  0  0  0 0  0  0  0  0;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  1 -1 -1  0  0  0  0 0  0  0  0  0;"


           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 0 -1 -1  0  0  0  0 0  0  0  0  1;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 0  0 -1  0  0  0  0 0  0  0  0  0;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 1 -1 -1  0  0  0  0 0  0  0  0  0;"

           // Run through alg:                                e                               
           // "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 -1  0  0 0  0  0  0  0 0  0  0  0  1;" costs
           // "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  1  0  0 0  0  0  0  0 0  0  0  0 -1;" n = 0
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0  0 -1 1  0  0  0  0 0  0  0  0  0;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  1 -1 -1 1  0  0  0  0 0  0  0  0  0;"


           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 0 -1 -1 0  0  0  0  0 0  0  0  0  1;" costs
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 0 -1 -1 1  0  0  0  0 0  0  0  0  1;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 0  0 -1 1  0  0  0  0 0  0  0  0  0;"
           // "1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 1 -1 -1 1  0  0  0  0 0  0  0  0  0;"
           // The k_src >= k_dst equation can be made basic          ^  in place of the bound  ^
           // What about k_dst == 0? What if we add this constraint?
    // clang-format on
    alloc::Arena<> alloc = managed_alloc;
    Simplex &simp{simplexFromTableau(&alloc, t)};
    expect(fatal(!simp.initiateFeasible()));
    std::ptrdiff_t num_nuisance = 26;
    Simplex::Solution sol = simp.rLexMinStop(num_nuisance);
    // sol.print();
    expect(allZero(sol));
    expect(eq(sol.size(), 5));
  };

  "DynsolveZeroExtra"_test = [&managed_alloc] -> void {
    // clang-format off
    // Approach 0: append dims so they match
    // We also muset set the new dims equal to the old, which seems wrong.
    auto t{"[0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 1 -1 0 -1 0 -1 0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0 -1 0  0  0  0;"
            "0 0  1 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  1 0  0 0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  1 0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  1  0 -1  0 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  1  0 -1 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0 -1 1  0 0  0  0 -1  0  1 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "1 0  0 0  0 0 -1 1 -1  0  0  0  0 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 0 -1 1  0 0  0 0  0 -1  0  1  0 0  0 0  0 0  0 0  0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 1 -1 0 -1 0 -1 0  0  0  0  0  0 0 -1  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 0  1 0  0 0  0 0  0  0  0  0  0 0  0 -1  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  1 0  0 0  0  0  0  0  0 0  0  0 -1  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  1 0  0  0  0  0  0 0  0  0  0 -1;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0 0  0  1  0 -1  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0  0 0  0  0  1  0 -1 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0 -1 1  0 0  0  0 -1  0  1 0  0  0  0  0;"
           "-1 0  0 0  0 0  0 0  0  0  0  0  0 0  0 0  0 0 -1 1 -1  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0 0  0  0  0  0  0 0 -1 1  0 0  0 0  0 -1  0  1  0 0  0  0  0  0]"_mat};

       // "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1  0  0  0 0 0 0 0 0 0 0 0 1;"
       // "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  1 -1 -1  1 0 0 0 0 0 0 0 0 0;"
       // 
       // "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1  1  0  0  0  0 0  0  0  0 1;"
       // "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 -1  1  0  0  0  0 0  0  0  0 0;"
    // clang-format on
    alloc::Arena<> alloc = managed_alloc;
    Simplex &simp{simplexFromTableau(&alloc, t)};
    expect(fatal(!simp.initiateFeasible()));
    std::ptrdiff_t num_nuisance = 24;
    Simplex::Solution sol = simp.rLexMinStop(num_nuisance);
    // sol.print();
    expect(allZero(sol));
    expect(eq(sol.size(), 5));
  };
  "Dynsolve0"_test = [&managed_alloc] -> void {
    // clang-format off
    // Approach 0: append dims so they match
    // We also muset set the new dims equal to the old, which seems wrong.
    //                                          sat
    auto t{"[  0 0  0  0  0  0  0  0  0  0  0  0  0  0;"
            "  0 1 -1 -1 -1  0 -1 -1 -1  0  0  0  0  0;" // 1
            "  0 0  1  0  1  0  1  0  1  0  0  0  0  0;" // N
            "  0 0 -1 -1  0  1  0  0  0  0  0  0  1  0;" // n
            "  0 0  0  1 -1  0  0  0  0  0  1 -1  0  1;" // k
            "  0 0  0  0  0  0 -1 -1  0  1 -1  1 -1  0;" // n'
            "  0 0  0  0  0  0  0  1 -1  0  0  0  0 -1;" // k'
            "  1 0  0  0  0  0  0  0  0  0  0  0  1  1]"_mat};

    // clang-format on
    alloc::Arena<> alloc = managed_alloc;
    Simplex &simp{simplexFromTableau(&alloc, t)};
    expect(fatal(!simp.initiateFeasible()));
    Simplex::Solution sol = simp.rLexMinStop(11);
    expect(eq(sol.size(), 2));
    expect(sol[last - 1] == 1);
    expect(sol[last] == 0);
  };
  "DynsolveCheckEmpty"_test = [&managed_alloc] -> void {
    // clang-format off
    // Approach 0: append dims so they match
    // We also muset set the new dims equal to the old, which seems wrong.
    //                                          sat
    auto d{"[  0 -1 -1 -1  0 -1 -1 -1  0  0  0;"       // 1
            "  0  1  0  1  0  1  0  1  0  0  0;"       // N
            "  0 -1 -1  0  1  0  0  0  0  0  0;"       // n
            "  0  0  1 -1  0  0  0  0  0  1 -1;"       // k
            "  0  0  0  0  0 -1 -1  0  1 -1  1;"       // n'
            "  0  0  0  0  0  0  1 -1  0  0  0]"_mat}; // k'
    // here, we let n = n'
    auto n{"[  0 -1 -1 -1  0 -1 -1 -1  0  0  0  0  0;"       // 1
            "  0  1  0  1  0  1  0  1  0  0  0  0  0;"       // N
            "  0 -1 -1  0  1  0  0  0  0  0  0  1 -1;"       // n
            "  0  0  1 -1  0  0  0  0  0  1 -1  0  0;"       // k
            "  0  0  0  0  0 -1 -1  0  1 -1  1 -1  1;"       // n'
            "  0  0  0  0  0  0  1 -1  0  0  0  0  0]"_mat}; // k'
    // here, we let n = n'
    auto k{"[  0 -1 -1 -1  0 -1 -1 -1  0  0  0  0  0;"       // 1
            "  0  1  0  1  0  1  0  1  0  0  0  0  0;"       // N
            "  0 -1 -1  0  1  0  0  0  0  0  0  0  0;"       // n
            "  0  0  1 -1  0  0  0  0  0  1 -1  1 -1;"       // k
            "  0  0  0  0  0 -1 -1  0  1 -1  1  0  0;"       // n'
            "  0  0  0  0  0  0  1 -1  0  0  0 -1  1]"_mat}; // k'

    // clang-format on
    for (int i = 0; i < 3; ++i) {
      alloc::Arena<> alloc = managed_alloc;
      PtrMatrix<std::int64_t> tableau;
      switch (i) {
      case 0: tableau = d; break;
      case 1: tableau = n; break;
      default: tableau = k; break;
      }
      Simplex &S{simplexFromTableau(&alloc, tableau)};
      expect(fatal(!S.initiateFeasible()));
      MutPtrVector<Simplex::value_type> C{S.getCost()};
      C.zero();
      C[_(1, 4)] << -1;
      C[_(5, 8)] << -1;
      if (i) expect(S.run() > 0); // means infeasible
      else expect(S.run() <= 0);  // means feasible
    }
  };
  "Dynsolve1"_test = [&managed_alloc] -> void {
    // clang-format off
    // Approach 0: append dims so they match
    // We also muset set the new dims equal to the old, which seems wrong.
    //                                          sat
    auto t{"[0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;"
            "0  1  0 -1 -1  0 -1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;" // 1
            "0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;" // N
            "0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0;" // n  src
            "0  0  0  1 -1  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  1  0  0;" // k  src
            "0  0  0  0  0  1 -1  0 -1  1  0  0  0  0  0  0  0  0  0 -1  0  0  0;" // n' dst
            "0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0;" // k' dst
            "0  0  0  0  0  0  0  0  0  0  1  0 -1 -1  0 -1 -1  0  0  0  0 -1  0;" // 1
            "0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  0  0  0 -1;" // N
            "0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0 -1  0  0  0;" // n
            "0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  1 -1  0 -1  0  0;" // k
            "0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0 -1  1  1  0  0  0;" // n'
            "0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  1  0  0;" // k'
            "1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0]"_mat};

    // clang-format on
    alloc::Arena<> alloc = managed_alloc;
    Simplex &simp{simplexFromTableau(&alloc, t)};
    expect(fatal(!simp.initiateFeasible()));
    Simplex::Solution sol = simp.rLexMinStop(18);
    expect(eq(sol.size(), 4));
    expect(sol[last - 3] == 1);
    expect(sol[last - 2] == 0);
    expect(sol[last - 1] == 0);
    expect(sol[last] == 1);
  };
  "Dynsolve2"_test = [&managed_alloc] -> void {
    // clang-format off
    // Approach 0: append dims so they match
    // We also muset set the new dims equal to the old, which seems wrong.
    //                                         phi   sat
    auto t{"[0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;"
            "0  1 -1 -2  0  0  0  0  0  0  0  0  0  0  0;" // 1
            "0  0  0  1  0  0  0  0  0  0  0  0  0  0  0;" // N
            "0  0  1 -1  1 -1  0  0  0  0  0  0  1  0  0;" // k  src
            "0  0  0  0 -1  1  0  0  0  0  0 -1  0  0  0;" // n' dst
            "0  0  0  0  0  0  1 -1 -2  0  0  0  0 -1  0;" // 1
            "0  0  0  0  0  0  0  0  1  0  0  0  0  0 -1;" // N
            "0  0  0  0  0  0  0  1 -1  1 -1  0 -1  0  0;" // k  src
            "0  0  0  0  0  0  0  0  0 -1  1  1  0  0  0;" // n' dst
            "1  0  0  0  0  0  0  0  0  0  0  1  1  0  0]"_mat};

    // clang-format on
    alloc::Arena<> alloc = managed_alloc;
    Simplex &simp{simplexFromTableau(&alloc, t)};
    expect(fatal(!simp.initiateFeasible()));
    Simplex::Solution sol = simp.rLexMinStop(10);
    expect(eq(sol.size(), 4));
    expect(2 * sol[last - 3] == 1);
    expect(2 * sol[last - 2] == 1);
    expect(sol[last - 1] == 0);
    expect(sol[last] == 0);
  };

  "DynsolveDropExtra"_test = [&managed_alloc] -> void {
    // clang-format off
    // Approach 1: drop dims so they match
    auto t{"[0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 1 -1 0 -1 0  0  0  0  0 0  0 0  0 0  0  0  0  0 -1 0  0  0  0;"
            "0 0  1 0  0 0  0  0  0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  1 0  0  0  0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  1  0 -1  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0  1  0 -1 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0 -1 1  0 -1  0  1 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"

            "0 0 -1 1  0 0 -1  0  1  0 0  0 0  0 0  0  0  0  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0  0  0  0 1 -1 0 -1 0  0  0  0  0 0 -1  0  0  0;"
            "0 0  0 0  0 0  0  0  0  0 0  1 0  0 0  0  0  0  0 0  0 -1  0  0;"
            "0 0  0 0  0 0  0  0  0  0 0  0 0  1 0  0  0  0  0 0  0  0 -1  0;"
            "0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0  0  0  0 0  0  0  0 -1;"
            "0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  1  0 -1  0 0  0  0  0  0;"
            "0 0  0 0  0 0  0  0  0  0 0  0 0  0 0  0  1  0 -1 0  0  0  0  0;"
            "0 0  0 0  0 0  0  0  0  0 0  0 0 -1 1  0 -1  0  1 0  0  0  0  0;"
            "0 0  0 0  0 0  0  0  0  0 0 -1 1  0 0 -1  0  1  0 0  0  0  0  0]"_mat};

    // clang-format on
    alloc::Arena<> alloc = managed_alloc;
    Simplex &simp{simplexFromTableau(&alloc, t)};
    expect(fatal(!simp.initiateFeasible()));
    std::ptrdiff_t num_nuisance = 18;
    Simplex::Solution sol = simp.rLexMinStop(num_nuisance);
    // sol.print();
    expect(allZero(sol));
    expect(eq(sol.size(), 5));
  };
  "OrderedComparison"_test = [&managed_alloc] -> void {
    // clang-format off
    //  0 = m - a - s0
    // s0 = m - a
    // s0 >= 0 -> m >= a
    // 
    //           m  a  b  x  y s0 s1 s2 s3
    auto t0{"[0  0  0  0  0  0  0  0  0  0;"
            " 0  1 -1  0  0  0 -1  0  0  0;"
            " 0  1  0 -1  0  0  0 -1  0  0;"
            " 0 -1  0  0  1  0  0  0 -1  0;"
            " 0 -1  0  0  0  1  0  0  0 -1]"_mat};
    // 0 =  m -a - s_0
    // 1 = -m + x - s_2
    auto t1{"[0  0  0  0  0  0  0  0  0  0;"
            " 0  1 -1  0  0  0 -1  0  0  0;"
            " 0  1  0 -1  0  0  0 -1  0  0;"
            " 1 -1  0  0  1  0  0  0 -1  0;"
            " 0 -1  0  0  0  1  0  0  0 -1]"_mat};
    // clang-format on
    static constexpr std::ptrdiff_t a = 2;
    static constexpr std::ptrdiff_t x = 4;
    static constexpr std::ptrdiff_t y = 5;

    for (int i = 0; i < 2; ++i) {
      alloc::Arena<> alloc = managed_alloc;
      Simplex &S{simplexFromTableau(&alloc, i ? t1 : t0)};
      expect(fatal(!S.initiateFeasible()));
      MutPtrVector<Simplex::value_type> C{S.getCost()};
      // is x >= a ?
      //   m  a b x y
      // [ 0 -1 0 1 0  ]
      C.zero(); // zeros the vector
      // minimize x - a
      C[a] = -1;
      C[x] = 1;
      // means x-a >= 0, i.e. x >= a
      // run returns negative, so value <= 0 means proven
      expect(S.run() == -i);
      C.zero();
      // minimize x - y
      C[x] = 1;
      C[y] = -1;
      expect(S.run() == std::numeric_limits<std::int64_t>::max());
    }
  };
  "CheckEmptyFarkas"_test = [&managed_alloc] -> void {
    // for (int i = 0; i < I; ++i)
    //   for (int j = 0; j < i; ++j)
    //     A[i,j] = A[j,i]
    // clang-format off
    //         
    auto t{"[0 -1 -1 -1 -1  0 0 0 0  0  0  0  0;"
            "0  1  0  1  0  0 0 0 0  0  0  0  0;"
            "0 -1  1  0  0  1 0 0 0  1  0 -1  0;"
            "0  0 -1  0  0  0 1 0 0  0  1  0 -1;"
            "0  0  0 -1  1  0 0 1 0  0 -1  0  1;"
            "0  0  0  0 -1  0 0 0 1 -1  0  1  0]"_mat};
    // y
    // check psi = x - y >= 0 everywhere in polyhedra
    // auto t1{"[0  0  0  0  0  0;"       // cost
    //         "0  1  1 -1 -1  0;"       // m
    //         "0 -1  0  0  0  0;"       // a
    //         "0  0 -1  0  0  0;"       // b
    //         "0  0  0  1  0 -1;"       // x
    //         "0  0  0  0  1  1]"_mat}; // y
    // clang-format on

    for (int i = 0; i < 2; ++i) {
      alloc::Arena<> alloc = managed_alloc;
      Simplex &S{simplexFromTableau(&alloc, t)};
      expect(!S.initiateFeasible());
      // S->run() returns the negative answer (fix that?)
      // if b'y < 0, A'y=0, y>=0, then there is no solution Ax<=b
      // Thus, we minimize b'y, if a value < 0 (s a return >0),
      // then there is no feasible solution
      MutPtrVector<Simplex::value_type> C{S.getCost()};
      C.zero();
      if (i) {
        // we make it infeasible, by setting j <= i-1
        C[_(1, 5)] << -1;
        expect(S.run() > 0);
      } else {
        // we make it feasible, by setting j <= i
        C[1] = -1;
        C[3] = -1;
        expect(S.run() <= 0);
      }
    }
  };

  return 0;
}

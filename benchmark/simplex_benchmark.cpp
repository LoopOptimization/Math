import Nanobench;
import Arena;
import ManagedArray;
import MatrixBinaryIO;
import Simplex;
import std;

void BM_Simplex0(Bench &bench) {
  auto tableau =
    utils::readMatrixBinary(MATH_DATA_DIR "/simplex_tableau_0.binmat");
  alloc::OwningArena<> alloc;
  std::ptrdiff_t numCon = std::ptrdiff_t(tableau.numRow()) - 1;
  std::ptrdiff_t numVar = std::ptrdiff_t(tableau.numCol()) - 1;
  math::Simplex &simpBackup{
    math::Simplex::create(&alloc, math::row(numCon), math::col(numVar))};
  simpBackup.getTableau() << tableau;
  // Simplex simpBackup{tableau};
  math::Simplex &simp{
    math::Simplex::create(&alloc, math::row(simpBackup.getNumCons()),
                          math::col(simpBackup.getNumVars()))};
  // Vector<Rational> sol(37);
  bench.run("BM_Simplex0", [&] {
    simp << simpBackup;
    bool fail = simp.initiateFeasible();
#ifndef NDEBUG
    if (fail) __builtin_trap();
#endif
    if (!fail) simp.rLexMinLast(37);
  });
  alloc.reset();
}

void BM_Simplex1(Bench &bench) {
  auto tableau =
    utils::readMatrixBinary(MATH_DATA_DIR "/simplex_tableau_1.binmat");

  alloc::OwningArena<> alloc;
  std::ptrdiff_t numCon = std::ptrdiff_t(tableau.numRow()) - 1;
  std::ptrdiff_t numVar = std::ptrdiff_t(tableau.numCol()) - 1;
  math::Simplex &simpBackup{
    math::Simplex::create(&alloc, math::row(numCon), math::col(numVar), 0)};
  simpBackup.getTableau() << tableau;
  math::Simplex &simp{
    math::Simplex::create(&alloc, math::row(simpBackup.getNumCons()),
                          math::col(simpBackup.getNumVars()), 0)};
  bench.run("BM_Simplex1", [&] {
    simp << simpBackup;
    bool fail = simp.initiateFeasible();
#ifndef NDEBUG
    if (fail) __builtin_trap();
#endif
    if (!fail) simp.rLexMinLast(15);
  });
}

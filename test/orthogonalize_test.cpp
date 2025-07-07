import boost.ut;
import Arena;
import Array;
import ManagedArray;
import MatDim;
import NormalForm;
import std;

using namespace boost::ut;
using ::math::DenseMatrix, ::math::DenseDims, ::math::row, ::math::col;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
void testBasicAssertions() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(-2, 2);

  constexpr std::ptrdiff_t M = 7;
  constexpr std::ptrdiff_t N = 7;
  std::int64_t mem[M * N * 2];
  alloc::OwningArena<> alloc;
  ::math::MutDensePtrMatrix<std::int64_t> A{mem, DenseDims<>{row(M), col(N)}};
  constexpr std::size_t iters = 1000;
  for (std::size_t i = 0; i < iters; ++i) {
    for (auto &&a : A) a = distrib(gen);
    std::ptrdiff_t r = ::math::NormalForm::rank(alloc, A);
    auto orth = ::math::orthogonalize(&alloc, A);
    expect(fatal(orth.numCol() == A.numCol()));
    expect(fatal(orth.numRow() == r));
    ::math::MutDensePtrMatrix<std::int64_t> B{mem + (M * N),
                                            DenseDims<>{row(r), col(r)}};
    // note, A'A is not diagonal
    // but AA' is
    B << orth * orth.t();
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
    for (std::ptrdiff_t m = 0; m < r; ++m)
      for (std::ptrdiff_t n = 0; n < r; ++n)
        if (m != n) expect(fatal((B[m, n]) == 0));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  }
}

int main() {
  "OrthogonalizeMatricesTest BasicAssertions"_test = [] {
    testBasicAssertions();
  };
  return 0;
}

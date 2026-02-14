import Testing;
import Arena;
import Array;
import ManagedArray;
import MatDim;
import NormalForm;
import std;

using namespace testing;
using ::math::DenseMatrix, ::math::DenseDims, ::math::row, ::math::col;

namespace {
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
void testBasicAssertions(unsigned int rng_init) {
  std::mt19937 gen(rng_init);
  std::uniform_int_distribution<> distrib(-2, 2);

  static constexpr std::ptrdiff_t M = 7;
  static constexpr std::ptrdiff_t N = 7;
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
} // namespace

auto main(int argc, char *argv[]) -> int {
  unsigned int rng_init;
  if (argc == 1) {
    std::random_device rd;
    rng_init = rd();
  } else rng_init = std::atoi(argv[1]);
  "OrthogonalizeMatricesTest BasicAssertions"_test = [=] -> void {
    testBasicAssertions(rng_init);
  };
  return 0;
}

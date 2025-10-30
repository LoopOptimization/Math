import boost.ut;

import Arena;
import ArrayParse;
import AxisTypes;
import Comparisons;
import CorePrint;
import ManagedArray;
import MatDim;
import SmallSparseMatrix;
import StaticArray;
import std;
import TinyVector;
import Tuple;
import BaseUtils;
import UniformScaling;

using namespace ::math;
using utils::operator""_mat;
using namespace boost::ut;
int main() {
  // NOLINTNEXTLINE(modernize-use-trailing-return-type)
  "SparseIndexingTest BasicAssertions"_test = [] {
    SmallSparseMatrix<std::int64_t> sparseA(row(3), col(4));
    utils::print("&Asparse = ", utils::TinyString(&sparseA), '\n');
    sparseA[0, 1] = 5;
    sparseA[1, 3] = 3;
    sparseA[2, 0] = -1;
    sparseA[2, 1] = 4;
    sparseA[2, 2] = -2;
    IntMatrix<> A = sparseA;
    {
      IntMatrix<> A2(DenseDims<>{row(3), col(4)});
      MutPtrMatrix<std::int64_t> MA2 = A2;
      MA2 << sparseA;
      expect(A == A2);
    }
    for (std::ptrdiff_t i = 0; i < 3; ++i)
      for (std::ptrdiff_t j = 0; j < 4; ++j)
        expect((A[i, j]) == (sparseA[i, j]));
    DenseMatrix<std::int64_t> B(DenseDims<>{row(4), col(5)});
    expect(!B.isSquare());
    B[0, 0] = 3;
    B[0, 1] = -1;
    B[0, 2] = 0;
    B[0, 3] = -5;
    B[0, 4] = 1;
    B[1, 0] = -4;
    B[1, 1] = 5;
    B[1, 2] = -1;
    B[1, 3] = -1;
    B[1, 4] = -1;
    B[2, 0] = 1;
    B[2, 1] = 2;
    B[2, 2] = -5;
    B[2, 3] = 2;
    B[2, 4] = 3;
    B[3, 0] = -2;
    B[3, 1] = 1;
    B[3, 2] = 2;
    B[3, 3] = -3;
    B[3, 4] = 5;
    ManagedArray<std::int64_t, DenseDims<3>> C{DenseDims<3>{Row<3>{}, col(5)}};
    C[0, 0] = -20;
    C[0, 1] = 25;
    C[0, 2] = -5;
    C[0, 3] = -5;
    C[0, 4] = -5;
    C[1, 0] = -6;
    C[1, 1] = 3;
    C[1, 2] = 6;
    C[1, 3] = -9;
    C[1, 4] = 15;
    C[last, _] << "[-21 17 6 -3 -11]"_mat;
    expect(A.numRow() == (A * B).numRow());
    expect(B.numCol() == (A * B).numCol());
    expect(C == A * B);
    {
      IntMatrix<> C2{A * B};
      utils::print("C=");
      C.print();
      utils::print("\nC2=");
      C2.print();
      utils::print('\n');
      expect(C == C2);
      IntMatrix<> Bt{B.t()};
      {
        IntMatrix<> At{
          DenseDims<>{::math::asrow(A.numCol()), ::math::ascol(A.numRow())}};
        At[_(0, end), _(0, end)] << A.t();
        std::ptrdiff_t k =
          std::min(std::ptrdiff_t(A.numCol()), std::ptrdiff_t(B.numRow()));
        C2 += At.t() * Bt.t()[_(k), _];
        expect(C * 2 == C2);
        expect(C == At.t() * B);
        expect(C == A * Bt.t());
        expect(C == At.t() * Bt.t());
      }
      C2 -= A * Bt.t();
      expect(C == C2);
    }
    std::int64_t i = 0;
    IntMatrix<> D{C};
    utils::print("C=");
    C.print();
    utils::print('\n');
    static_assert(std::same_as<decltype(D[0, _]), MutPtrVector<std::int64_t>>);
    for (std::ptrdiff_t r : _(0, D.numRow())) D[r, _] += std::ptrdiff_t(r) + 1;
    for (auto r : C.eachRow()) {
      expect(r.size() == std::ptrdiff_t(C.numCol()));
      r += (++i);
    }
    expect(C == D);
    auto oldD{D};
    for (std::ptrdiff_t c : _(0, D.numCol()))
      D[_, c] += std::ptrdiff_t(c) + std::ptrdiff_t(D.numRow()) + 1;
    // for (auto c : C.eachCol()) c += (++i);
    // test structured binding
    for (auto &&[a, b, c] : C.eachCol()) {
      a += (++i);
      b += i;
      c += i;
    }
    expect(C == D);
    for (auto c : C.eachCol() | std::views::reverse) c -= (i--);
    expect(C == oldD);
    MutPtrVector<std::int64_t> v0 = C[0, _], v1 = C[1, _];
    containers::Pair(v0, v1) << containers::Pair(v1, v0);
    expect(v0 == oldD[1, _]);
    expect(v1 == oldD[0, _]);
  };

  // NOLINTNEXTLINE(modernize-use-trailing-return-type)
  "ExpressionTemplateTest BasicAssertions"_test = [] {
    // 6x8
    auto A{
      "[3 -5 1 10 -4 6 4 4; 4 6 3 -1 6 1 -4 0; -7 -2 0 0 -10 -2 3 7; 2 -7 -5 "
      "-5 -7 -5 1 -7; 2 -8 2 7 4 9 6 -3; -2 -8 -5 0 10 -4 5 -3]"_mat};

    auto A4{
      "[12 -20 4 40 -16 24 16 16; 16 24 12 -4 24 4 -16 0; -28 -8 0 0 -40 -8 "
      "12 28; 8 -28 -20 -20 -28 -20 4 -28; 8 -32 8 28 16 36 24 -12; -8 -32 "
      "-20 0 40 -16 20 -12]"_mat};
    static_assert(A4[0, _].size() == 8);
    static_assert(A4[_(::math::Row<3>()), _].size() == 24);
    // IntMatrix B{A*4};
    auto templateA4{A * 4};
    IntMatrix<> C{templateA4};
    IntMatrix<> B{A * 4};
    expect(A4 == B);
    expect(A4 == C);
    IntMatrix<> Z = A * 4 - A4;
    for (std::ptrdiff_t i = 0; i < Z.numRow(); ++i)
      for (std::ptrdiff_t j = 0; j < Z.numCol(); ++j) expect(!(Z[i, j]));
    auto D{
      "[-5 6 -1 -4 7 -9 6; -3 -5 -1 -2 -9 -4 -1; -4 7 -6 10 -2 2 9; -4 -7 -1 "
      "-7 5 9 -10; 5 -7 -5 -1 -3 -8 -8; 3 -6 4 10 9 0 -5; 0 -1 4 -4 -9 -3 "
      "-10; 2 1 4 5 -7 0 -8]"_mat};
    auto refAD{
      "[-38 -28 62 6 116 105 -138; -13 -22 -69 29 -10 -99 42; -1 54 91 45 "
      "-95 142 -36; -13 118 31 -91 78 8 151; 19 -74 15 26 153 31 -145; 86 "
      "-61 -18 -111 -22 -55 -135]"_mat};
    IntMatrix<> AD = A * D;
    expect(AD == refAD);
    IntMatrix<> E{
      "[-4 7 9 -4 2 9 -8; 3 -5 6 0 -1 8 7; -7 9 -1 1 -5 2 10; -3 10 -10 -3 6 "
      "5 5; -6 7 -4 -7 10 5 3; 9 -8 7 9 2 2 6]"_mat};
    IntMatrix<> m7EpAD = A * D - 7 * E;
    auto refADm7E{
      "[-10 -77 -1 34 102 42 -82; -34 13 -111 29 -3 -155 -7; 48 -9 98 38 -60 "
      "128 -106; 8 48 101 -70 36 -27 116; 61 -123 43 75 83 -4 -166; 23 -5 "
      "-67 -174 -36 -69 -177]"_mat};
    expect(m7EpAD == refADm7E);

    Vector<std::int64_t> a{-8};
    a.push_back(7);
    a.push_back(3);
    Vector<std::int64_t> b = a * 2;
    Vector<std::int64_t> c;
    c.push_back(-16);
    c.push_back(14);
    c.push_back(6);
    expect(b == c);
    c.resize(6);
    c << A[_, 1];
    Vector<std::int64_t> d(std::array<std::int64_t, 6>{-5, 6, -2, -7, -8, -8});
    expect(c == d);
    expect(b * c[_(0, 3)].t() == 152);
    IntMatrix<> dA1x1(DenseDims<>{row(1), col(1)}, 0);
    expect(dA1x1.isSquare());
    IntMatrix<> dA2x2(DenseDims<>{row(2), col(2)}, 0);
    dA1x1.antiDiag() << 1;
    expect((dA1x1[0, 0]) == 1);
    dA2x2.antiDiag() << 1;
    expect((dA2x2[0, 0]) == 0);
    expect((dA2x2[0, 1]) == 1);
    expect((dA2x2[1, 0]) == 1);
    expect((dA2x2[1, 1]) == 0);
    for (std::ptrdiff_t i = 1; i < 20; ++i) {
      IntMatrix<> F(DenseDims<>{row(i), col(i)});
      F << 0;
      F.antiDiag() << 1;
      for (std::ptrdiff_t j = 0; j < i; ++j)
        for (std::ptrdiff_t k = 0; k < i; ++k)
          expect(fatal(F[j, k] == (k + j == i - 1)));
    }
  };
  // NOLINTNEXTLINE(modernize-use-trailing-return-type)
  "ExpressionTemplateTest2 BasicAssertions"_test = [] {
    ManagedArray<double, DenseDims<>> W{{row(3), col(3)}, 0},
      X{{row(3), col(3)}, 0}, Y{{row(3), col(3)}, 0}, Z{{row(3), col(3)}, 0};
    W[0, 0] = 0.29483432115939806;
    W[0, 1] = 1.5777027461040212;
    W[0, 2] = 0.8171761007267028;
    W[1, 0] = 1.0463632179853855;
    W[1, 1] = 0.9503214631611095;
    W[1, 2] = -0.17890983978584624;
    W[2, 0] = 1.5853551451194254;
    W[2, 1] = -0.784875301203305;
    W[2, 2] = 1.7033024094365752;
    X[0, 0] = -1.1175097244313117;
    X[0, 1] = -0.21769215316295054;
    X[0, 2] = -0.7340630927749082;
    X[1, 0] = -0.5750426169922397;
    X[1, 1] = 0.27174064995044767;
    X[1, 2] = -1.0669896577273217;
    X[2, 0] = 0.9302424251181362;
    X[2, 1] = -1.3157431480603476;
    X[2, 2] = 1.546836705770486;
    Y[0, 0] = 1.1701212478097331;
    Y[0, 1] = 0.7747688878004019;
    Y[0, 2] = -0.926815554991563;
    Y[1, 0] = -1.4441713498640656;
    Y[1, 1] = -1.3615487160168993;
    Y[1, 2] = 0.7908183008408143;
    Y[2, 0] = -0.7626497248468547;
    Y[2, 1] = -0.21682371102755368;
    Y[2, 2] = -0.07604892144743511;
    Z[0, 0] = 3.3759933640164708e16;
    Z[0, 1] = 9.176788687153845e14;
    Z[0, 2] = -1.1081818546676994e15;
    Z[1, 0] = -1.7207794047001762e15;
    Z[1, 1] = 3.0768637505289172e16;
    Z[1, 2] = 9.277082601207064e14;
    Z[2, 0] = -8.956589651911538e14;
    Z[2, 1] = -2.7136623944168075e14;
    Z[2, 2] = 3.2308470074953084e16;
    ManagedArray<double, DenseDims<>> A{
      W * (W + 16380 * X + 40840800 * Y) +
      (33522128640 * W + 10559470521600 * X + 1187353796428800 * Y) +
      32382376266240000 * I};

    expect(A == Z);
  };

  // NOLINTNEXTLINE(modernize-use-trailing-return-type)
  "OffsetEnd BasicAssertions"_test = [] {
    auto A{"[3 3 3 3; 2 2 2 2; 1 1 1 1; 0 0 0 0]"_mat};
    auto B = IntMatrix<>{DenseDims<>{row(4), col(4)}};
    for (std::ptrdiff_t i = 0; i < 4; ++i) B[last - i, _] << i;
    expect(A == B);
  };
  "SetToIdentity"_test = [] {
    auto A{"[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]"_mat};
    auto B = IntMatrix<>{DenseDims<>{row(4), col(4)}};
    B << I;
    expect(A == B);
  };
  "SquareMatrixTest BasicAssertions"_test = [] {
    SquareMatrix<std::int64_t> A{SquareDims<>{row(4)}};
    for (std::ptrdiff_t i = 0; i < 4; ++i)
      for (std::ptrdiff_t j = 0; j < 4; ++j) A[i, j] = 4 * i + j;
    DenseMatrix<std::int64_t> B{DenseDims<>{row(4), col(2)}};
    B << A[_(end - 2, end), _].t();
    for (std::ptrdiff_t j = 0; j < 4; ++j)
      for (std::ptrdiff_t i = 0; i < 2; ++i)
        expect((B[j, i]) == 4 * (i + 2) + j);
  };
  "VectorTest BasicAssertions"_test = [] {
    alloc::OwningArena<> alloc;
    ResizeableView<std::int64_t, Length<>> x;
    for (std::size_t i = 0; i < 100; ++i) {
      if (x.getCapacity() <= x.size())
        x.reserve(&alloc, std::max<std::ptrdiff_t>(8, 2 * x.size()));
      x.emplace_back_within_capacity(i);
    }
    expect(x.size() == 100);
    expect(x.sum() == 100 * 99 / 2);
  };
  "SVectorTest BasicAssertions"_test = [] {
    SVector<std::int64_t, 3> x{1, 2, 3};
    // static_assert(simd::VecLen<3, std::int64_t> == 2);
    static_assert(utils::Compressible<SVector<std::int64_t, 3>>);
    static_assert(std::tuple_size_v<decltype(x)> == 3);
    static_assert(
      std::same_as<std::tuple_element_t<2, decltype(x)>, std::int64_t>);
    SVector<std::int64_t, 3> y{10, 20, 30};
    SVector<std::int64_t, 3, true> ycompress{10, 20, 30};
    y = ycompress;
    SVector<std::int64_t, 3> z{11, 22, 33};
    SVector<std::int64_t, 3> w = x + y;
    expect(w == z);
    expect(w == z);
    constexpr auto const_cmp = [](auto const &a, auto const &b) {
      return containers::tuple(a == b, std::is_constant_evaluated());
    };
    Vector<std::int64_t> v{::math::length(3)};
    v << _(1, 4);
    expect(const_cmp(v.size(), 3)._0);
    expect(!const_cmp(v.size(), 3)._1);
    expect(const_cmp(x.size(), 3)._0);
    expect(std::distance(
             v.begin(),
             std::ranges::find_if(v[_(1, end)],
                                  std::bind_front(std::equal_to<>{}, 3))) == 2);
    // expect(constCmp(v.size())._0);
    // expect(!(constCmp(v.size())._1);
    // expect(constCmp(x.size())._0);
    static_assert(const_cmp(decltype(x)::size(), unsigned(3))._1);
    static_assert(const_cmp(decltype(x)::size(), decltype(y)::size())._1);
    // expect(constCmp(x.size())._1);
    // expect(constCmp(x.size(), unsigned(3))._1);
    // expect(constCmp(x.size(), y.size())._1);
    auto [a, b, c] = w;
    expect(a == 11);
    expect(b == 22);
    expect(c == 33);
  };
  "TinyVectorTest BasicAssertions"_test = [] {
    {
      containers::TinyVector<int, 5> v{};
      static_assert(std::same_as<utils::eltype_t<decltype(v)>, int>);
      expect(v.empty());
      expect(v.size() == 0);
      v.resize(3);
      expect(!(v.empty()));
      expect(v.size() == 3);
      expect(v.back() == 0);
      v.push_back(2);
      expect(v.size() == 4);
      expect(v.back() == 2);
      expect(v.pop_back_val() == 2);
      expect(v.front() == 0);
      expect(v.back() == 0);
      expect(v.size() == 3);
      v.pop_back();
      expect(v.size() == 2);
      v.pop_back();
      expect(!(v.empty()));
      expect(v.size() == 1);
      v.pop_back();
      expect(v.empty());
      expect(v.size() == 0);
      int &y = v.emplace_back(2);
      y += 3;
      expect(v.front() == 5);
      expect(v.back() == 5);
      v.push_back(2);
      expect(v.front() == 5);
      expect(v.back() == 2);
      v.push_back(21);
      expect(v.back() == 21);
      int s = 0;
      for (auto x : v) s += x;
      expect(s == 28);
    }
    {
      containers::TinyVector<std::int8_t, 5, std::int8_t> v{};
      static_assert(std::same_as<utils::eltype_t<decltype(v)>, std::int8_t>);
      expect(v.empty());
      expect(v.size() == 0);
      v.resize(3);
      expect(!(v.empty()));
      expect(v.size() == 3);
      expect(v.back() == 0);
      v.push_back(2);
      expect(v.size() == 4);
      expect(v.back() == 2);
      expect(v.pop_back_val() == 2);
      expect(v.front() == 0);
      expect(v.back() == 0);
      expect(v.size() == 3);
      v.pop_back();
      expect(v.size() == 2);
      v.pop_back();
      expect(!(v.empty()));
      expect(v.size() == 1);
      v.pop_back();
      expect(v.empty());
      expect(v.size() == 0);
      std::int8_t &y = v.emplace_back(2);
      y += 3;
      expect(v.front() == 5);
      expect(v.back() == 5);
      v.push_back(2);
      expect(v.front() == 5);
      expect(v.back() == 2);
      v.push_back(21);
      expect(v.back() == 21);
      std::int8_t s = 0;
      for (auto x : v) s = std::int8_t(s + x);
      expect(s == 28);
    }
  };

  "NonTriviallyDestructible BasicAssertions"_test = [] {
    // 2 + 2*100, 3 + 2*100, 4 + 2*100
    Vector<std::int64_t> y{std::array<std::int64_t, 3>{204, 205, 206}};
    for (std::ptrdiff_t i = 0; i < 4; ++i) {
      auto [a, b] = y.split(i);
      expect(a == y[_(0, i)]);
      expect(b == y[_(i, end)]);
    }
    {
      auto [a, b] = y.popFront();
      expect(a == y[0]);
      expect(b == y[_(1, end)]);
    }
    Vector<std::int64_t> z{std::array<std::int64_t, 3>{0, 1, 2}};
    Vector<Vector<std::int64_t, 0>, 0> x{::math::length(5)};
    for (std::ptrdiff_t i = 0; i < 10; i += 2)
      x[i / 2] = std::array<std::int64_t, 3>{2 + 2 * i, 3 + 2 * i, 4 + 2 * i};
    for (std::ptrdiff_t i = 10; i < 102; i += 2)
      x.emplace_back(
        std::array<std::int64_t, 3>{2 + 2 * i, 3 + 2 * i, 4 + 2 * i});
    for (std::ptrdiff_t i = 1; i < 102; i += 2)
      x.insert(x.begin() + i,
               std::array<std::int64_t, 3>{2 + 2 * i, 3 + 2 * i, 4 + 2 * i});
    for (std::ptrdiff_t i = 0; i < x.size(); ++i)
      for (std::ptrdiff_t j = 0; j < 3; ++j) expect(x[i][j] == 2 * (i + 1) + j);
    expect(x.pop_back_val() == y);
    x.truncate(55);
    x.resize(45);
    z += 2 * x.size();
    expect(x.pop_back_val() == z);
    x.resizeForOverwrite(23);
    z -= 2 * (45 - x.size());
    expect(x.pop_back_val() == z);
    x.zero();
    expect(allZero(x));
    y.resize(2);
    y << 1;
    expect(allZero(y - 1));
  };

  "StringMat1x1 BasicAssertions"_test = [] {
    auto B = "[-5]"_mat;
    PtrMatrix<std::int64_t> Bp = B;
    IntMatrix<> A = "[-5]"_mat;
    ::math::DensePtrMatrix<std::int64_t> Dp = B;
    ::math::ManagedArray<std::int64_t, ::math::DenseDims<>> D = "[-5]"_mat;
    containers::Tuple{Bp, A, Dp, D}.apply([](const auto &x) {
      expect(std::ptrdiff_t(x.numCol()) == 1);
      expect(std::ptrdiff_t(x.numRow()) == 1);
      expect((x[0, 0]) == -5);
    });
  };

  "StringVector BasicAssertions"_test = [] {
    static_assert(!utils::Compressible<std::int64_t>);
    auto a = "[-5 3 7]"_mat;
    auto along = "[-5 3 7 -15 17 -5 -4 -3 -2 1 0 0 1 2 0 3 4 5 6 7]"_mat;
    static_assert(
      std::convertible_to<::math::StaticDims<std::int64_t, 1, 3, false>,
                          Length<>>);
    ::math::Array<std::int64_t, ::math::StaticDims<std::int64_t, 1, 3, false>>
      aps = a;
    PtrVector<std::int64_t> ap = a;
    Vector<std::int64_t> b = a;
    expect("[-5 3]"_mat == a[_(0, 2)]);
    const auto &ca = along;
    expect(a == ca[_(0, 3)]);
    containers::Tuple{aps, ap, b}.apply([](const auto &x) {
      expect(std::ptrdiff_t(x.size()) == 3);
      expect(x[0] == -5);
      expect(x[1] == 3);
      expect(x[2] == 7);
    });
    expect(anyNEZero(a));
    expect(anyGTZero(a));
    expect(anyLTZero(a));
  };

  "DimensionConversions BasicAssertions"_test = [] {
    // Test StridedDims conversions
    {
      // StridedDims<> -> StridedDims<4, 3, 4> (runtime to compile-time)
      StridedDims<> runtime_strided{row(4), col(3), stride(4)};
      StridedDims<4, 3, 4> static_strided = runtime_strided;
      expect(row(static_strided) == 4);
      expect(col(static_strided) == 3);
      expect(stride(static_strided) == 4);

      // StridedDims<4, 3, 4> -> StridedDims<> (compile-time to runtime)
      StridedDims<> back_to_runtime = static_strided;
      expect(row(back_to_runtime) == 4);
      expect(col(back_to_runtime) == 3);
      expect(stride(back_to_runtime) == 4);

      // Test with Arrays
      ManagedArray<std::int64_t, StridedDims<>> array_runtime{runtime_strided};
      ManagedArray<std::int64_t, StridedDims<4, 3, 4>> array_static{
        static_strided};
      for (std::ptrdiff_t i = 0; i < 4; ++i) {
        for (std::ptrdiff_t j = 0; j < 3; ++j) {
          array_runtime[i, j] = i * 10 + j;
          array_static[i, j] = i * 10 + j;
        }
      }
      expect((array_runtime[2, 1]) == 21);
      expect((array_static[2, 1]) == 21);
    }

    // Test DenseDims conversions
    {
      // DenseDims<> -> DenseDims<3, 4> (runtime to compile-time)
      DenseDims<> runtime_dense{row(3), col(4)};
      DenseDims<3, 4> static_dense = runtime_dense;
      expect(row(static_dense) == 3);
      expect(col(static_dense) == 4);

      // DenseDims<3, 4> -> DenseDims<> (compile-time to runtime)
      DenseDims<> back_to_runtime = static_dense;
      expect(row(back_to_runtime) == 3);
      expect(col(back_to_runtime) == 4);

      // Test with Arrays
      ManagedArray<std::int64_t, DenseDims<>> array_runtime{runtime_dense};
      ManagedArray<std::int64_t, DenseDims<3, 4>> array_static{static_dense};
      for (std::ptrdiff_t i = 0; i < 3; ++i) {
        for (std::ptrdiff_t j = 0; j < 4; ++j) {
          array_runtime[i, j] = i * 100 + j;
          array_static[i, j] = i * 100 + j;
        }
      }
      expect((array_runtime[1, 2]) == 102);
      expect((array_static[1, 2]) == 102);
    }

    // Test SquareDims conversions
    {
      // SquareDims<> -> SquareDims<5> (runtime to compile-time)
      SquareDims<> runtime_square{row(5)};
      SquareDims<5> static_square = runtime_square;
      expect(row(static_square) == 5);
      expect(col(static_square) == 5);

      // SquareDims<5> -> SquareDims<> (compile-time to runtime)
      SquareDims<> back_to_runtime = static_square;
      expect(row(back_to_runtime) == 5);
      expect(col(back_to_runtime) == 5);

      // Test with Arrays
      ManagedArray<std::int64_t, SquareDims<>> array_runtime{runtime_square};
      ManagedArray<std::int64_t, SquareDims<5>> array_static{static_square};
      for (std::ptrdiff_t i = 0; i < 5; ++i) {
        for (std::ptrdiff_t j = 0; j < 5; ++j) {
          array_runtime[i, j] = i * j;
          array_static[i, j] = i * j;
        }
      }
      expect((array_runtime[3, 4]) == 12);
      expect((array_static[3, 4]) == 12);
    }

    // Test cross-type conversions
    {
      // SquareDims -> DenseDims
      SquareDims<3> square_3x3;
      DenseDims<3, 3> dense_from_square = square_3x3;
      DenseDims<> dense_runtime_from_square = square_3x3;
      expect(row(dense_from_square) == 3);
      expect(col(dense_from_square) == 3);
      expect(row(dense_runtime_from_square) == 3);
      expect(col(dense_runtime_from_square) == 3);

      // DenseDims -> StridedDims
      DenseDims<2, 4> dense_2x4;
      StridedDims<2, 4, 4> strided_from_dense = dense_2x4;
      StridedDims<> strided_runtime_from_dense = dense_2x4;
      expect(row(strided_from_dense) == 2);
      expect(col(strided_from_dense) == 4);
      expect(stride(strided_from_dense) == 4);
      expect(row(strided_runtime_from_dense) == 2);
      expect(col(strided_runtime_from_dense) == 4);
      expect(stride(strided_runtime_from_dense) == 4);

      // SquareDims -> StridedDims
      SquareDims<6> square_6x6;
      StridedDims<6, 6, 6> strided_from_square = square_6x6;
      StridedDims<> strided_runtime_from_square = square_6x6;
      expect(row(strided_from_square) == 6);
      expect(col(strided_from_square) == 6);
      expect(stride(strided_from_square) == 6);
      expect(row(strided_runtime_from_square) == 6);
      expect(col(strided_runtime_from_square) == 6);
      expect(stride(strided_runtime_from_square) == 6);
    }

    // Test partial conversions (some known, some unknown parameters)
    {
      // StridedDims<-1, 3, -1> -> StridedDims<4, 3, 8>
      StridedDims<-1, 3, -1> partial_strided{row(4), col(3), stride(8)};
      StridedDims<4, 3, 8> full_strided = partial_strided;
      expect(row(full_strided) == 4);
      expect(col(full_strided) == 3);
      expect(stride(full_strided) == 8);

      // DenseDims<2, -1> -> DenseDims<2, 7>
      DenseDims<2, -1> partial_dense{row(2), col(7)};
      DenseDims<2, 7> full_dense = partial_dense;
      expect(row(full_dense) == 2);
      expect(col(full_dense) == 7);

      // Test with Arrays
      ManagedArray<std::int64_t, DenseDims<2, -1>> array_partial{partial_dense};
      ManagedArray<std::int64_t, DenseDims<2, 7>> array_full{full_dense};
      for (std::ptrdiff_t i = 0; i < 2; ++i) {
        for (std::ptrdiff_t j = 0; j < 7; ++j) {
          array_partial[i, j] = i + j;
          array_full[i, j] = i + j;
        }
      }
      expect((array_partial[1, 5]) == 6);
      expect((array_full[1, 5]) == 6);
    }
  };

  "MutPtrMatrixBool BasicAssertions"_test = [] {
    // Test MutPtrMatrix<bool> functionality
    ManagedArray<bool, DenseDims<>> boolArray{DenseDims<>{row(3), col(4)},
                                              false};
    MutPtrMatrix<bool> boolMatrix = boolArray;

    // Test assignment and access
    boolMatrix[0, 0] = true;
    boolMatrix[1, 2] = true;
    boolMatrix[2, 3] = true;

    expect(boolMatrix[0, 0] == true);
    expect(boolMatrix[1, 2] == true);
    expect(boolMatrix[2, 3] == true);
    expect(boolMatrix[0, 1] == false);
    expect(boolMatrix[1, 1] == false);

    // Test row and column access
    expect(boolMatrix.numRow() == 3);
    expect(boolMatrix.numCol() == 4);

    // Test setting entire row
    boolMatrix[2, _] << true;
    for (std::ptrdiff_t j = 0; j < 4; ++j) expect(boolMatrix[2, j] == true);

    // Test copying from another matrix
    ManagedArray<bool, DenseDims<>> boolArray2{DenseDims<>{row(3), col(4)},
                                               true};
    MutPtrMatrix<bool> boolMatrix2 = boolArray2;
    boolMatrix2[1, 1] = false;

    expect(boolMatrix2[0, 0] == true);
    expect(boolMatrix2[1, 1] == false);
    expect(boolMatrix2[2, 2] == true);
  };
  "Diag Test"_test = [] {
    DenseMatrix<std::int64_t> A{DenseDims<>{row(15), col(15)}},
      B{DenseDims<>{row(15), col(15)}}, C{DenseDims<>{row(15), col(15)}},
      D{DenseDims<>{row(15), col(15)}}, E{DenseDims<>{row(15), col(15)}},
      F{DenseDims<>{row(15), col(15)}};
    A.zero(), B.zero(), C.zero(), D.zero(), E.zero(), F.zero();
    using containers::Tuple;
    Tuple(A.diag(), B.diag(), C.diag(), D.diag(), E.diag(), F.diag())
      << Tuple(1, -2, 3, -4, 5, -6);
    // 1 - 2 + 6 - 4 + 5 - 6 == 0
    expect(allZero(A + B + 2 * C + D + E + F));
    // 1 + 2 + 3 - 4 + 5 - 6 == 1
    expect(A - B + C + D + E + F == I);
  };
  return 0;
}

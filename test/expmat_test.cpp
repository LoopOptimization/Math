import boost.ut;
import Array;
import ArrayConcepts;
import Dual;
import LinearAlgebra;
import ManagedArray;
import MatDim;
import Reductions;
import StaticArray;
import std;
import TinyVector;
import Tuple;
import BaseUtils;
import UniformScaling;

namespace {
using namespace ::math;
using utils::eltype_t;
using namespace boost::ut;

template <typename T>
constexpr void evalpoly(MutSquarePtrMatrix<T> B, MutSquarePtrMatrix<T> A,
                        SquarePtrMatrix<T> C, const auto &p) {
  std::ptrdiff_t N = p.size();
  invariant(N > 0);
  invariant(std::ptrdiff_t(B.numRow()), std::ptrdiff_t(C.numRow()));
  if (N & 1) std::swap(A, B);
  B << p[0] * C + p[1] * I;
  for (std::ptrdiff_t i = 2; i < N; ++i) {
    std::swap(A, B);
    B << A * C + p[i] * I;
  }
}

template <AbstractMatrix T> constexpr auto opnorm1(const T &A) {
  using S = decltype(extractvalue(std::declval<eltype_t<T>>()));
  auto [M, N] = shape(A);
  invariant(M > 0);
  invariant(N > 0);
  S a{};
  for (std::ptrdiff_t n = 0; n < N; ++n) {
    S s{};
    for (std::ptrdiff_t m = 0; m < M; ++m) s += std::abs(extractvalue(A[m, n]));
    a = std::max(a, s);
  }
  return a;
}

/// computes ceil(log2(x)) for x >= 1
constexpr auto log2ceil(double x) -> unsigned {
  invariant(x >= 1);
  std::uint64_t u = std::bit_cast<std::uint64_t>(x) - 1;
  return (u >> 52) - 1022;
}

template <typename T> constexpr void expmimpl(MutSquarePtrMatrix<T> A) {
  std::ptrdiff_t n = std::ptrdiff_t(A.numRow()), s = 0;
  SquareMatrix<T> A2{SquareDims<>{row(n)}}, U_{SquareDims<>{row(n)}};
  MutSquarePtrMatrix<T> U{U_};
  if (double nA = opnorm1(A); nA <= 0.015) {
    A2 << A * A;
    U << A * (A2 + 60.0 * I);
    A << 12.0 * A2 + 120.0 * I;
  } else {
    SquareMatrix<T> B{SquareDims<>{row(n)}};
    if (nA <= 2.1) {
      A2 << A * A;
      containers::TinyVector<double, 5> p0, p1;
      if (nA > 0.95) {
        p0 = {1.0, 3960.0, 2162160.0, 302702400.0, 8821612800.0};
        p1 = {90.0, 110880.0, 3.027024e7, 2.0756736e9, 1.76432256e10};
      } else if (nA > 0.25) {
        p0 = {1.0, 1512.0, 277200.0, 8.64864e6};
        p1 = {56.0, 25200.0, 1.99584e6, 1.729728e7};
      } else {
        p0 = {1.0, 420.0, 15120.0};
        p1 = {30.0, 3360.0, 30240.0};
      }
      evalpoly(B, U, A2, p0);
      U << A * B;
      evalpoly(A, B, A2, p1);
    } else {
      // s = std::max(unsigned(std::ceil(std::log2(nA / 5.4))), 0);
      s = nA > 5.4 ? log2ceil(nA / 5.4) : 0;
      if (s & 1) {       // we'll swap `U` and `A` an odd number of times
        std::swap(A, U); // so let them switch places
        A << U * exp2(-s);
      } else if (s > 0) A *= exp2(-s);
      A2 << A * A;
      // here we take an estrin (instead of horner) approach to cut down flops
      SquareMatrix<T> A4{A2 * A2}, A6{A2 * A4};
      B << A6 * (A6 + 16380 * A4 + 40840800 * A2) +
             (33522128640 * A6 + 10559470521600 * A4 + 1187353796428800 * A2) +
             32382376266240000 * I;
      U << A * B;
      A << A6 * (182 * A6 + 960960 * A4 + 1323241920 * A2) +
             (670442572800 * A6 + 129060195264000 * A4 +
              7771770303897600 * A2) +
             64764752532480000 * I;
    }
  }
  containers::tie(A, U) << containers::Tuple(A + U, A - U);
  LU::ldiv(U, MutPtrMatrix<T>(A));
  for (; s--; std::swap(A, U)) U << A * A;
}

template <typename T> constexpr auto expm(SquarePtrMatrix<T> A) {
  SquareMatrix<T> V{SquareDims{A.numRow()}};
  V << A;
  expmimpl(V);
  return V;
}
constexpr auto dualDeltaCmp(double x, double y) -> bool { return x < y; }
template <typename T, std::ptrdiff_t N>
constexpr auto dualDeltaCmp(Dual<T, N> x, double y) -> bool {
  if (!dualDeltaCmp(x.value(), y)) return false;
  for (std::ptrdiff_t i = 0; i < N; ++i)
    if (!dualDeltaCmp(x.gradient()[i], y)) return false;
  return true;
}
} // namespace
auto main() -> int {
  "ExpMatTest BasicAssertions"_test = [] -> void {
    static constexpr double adata[16]{
      0.13809508135032297, -0.10597225613986219, -0.5623996136438215,
      1.099556072129511,   0.7571409301354933,   -0.0725924122459707,
      -0.5732592019339723, 0.4216707913809331,   0.9551223749499392,
      1.0628072168060698,  -1.0919065664748313,  0.9125836172498181,
      -1.4826804140146677, 0.6550780207685463,   -0.6227535845719466,
      0.2280514374580733};

    static constexpr double bdata[16] = {
      0.2051199361909877,  -0.049831094437687434, -0.3980657896416266,
      0.6706580677244947,  0.14286173961464693,   0.7798855526203928,
      -0.468198822464851,  0.3841802990849566,    0.04397695538798724,
      0.7186280674937524,  -0.10256382048628668,  0.8710856559160713,
      -1.2481676162786608, 0.4472989810307132,    -0.11106692926404803,
      0.3930685232252409};

    static constexpr double addata[48] = {
      0.13809508135032297,  0.23145585885555967,  0.6736099502056541,
      -0.10597225613986219, -1.697508191315446,   -1.0754726189887889,
      -0.5623996136438215,  -0.1193943011160536,  -0.26291486453222596,
      1.099556072129511,    1.387881995203853,    3.1737580685186204,
      0.7571409301354933,   -1.7465876925755974,  -3.0307435092366,
      -0.0725924122459707,  0.3308073172514302,   0.07339531064172199,
      -0.5732592019339723,  -0.6947698845697236,  -1.0672845574608438,
      0.4216707913809331,   -0.37974214526732475, 0.23640427332714647,
      0.9551223749499392,   0.6139406820054037,   1.5323159670356596,
      1.0628072168060698,   -1.9568771179939957,  -0.4461533904276679,
      -1.0919065664748313,  1.1814457776729768,   -1.8787170755674358,
      0.9125836172498181,   1.2012117345451894,   -2.031045206544551,
      -1.4826804140146677,  1.8971407850439679,   0.05411248868441574,
      0.6550780207685463,   -1.121131948533708,   0.3543625073973555,
      -0.6227535845719466,  1.195008337302074,    -1.4233657256047536,
      0.2280514374580733,   -1.2001994532706792,  0.03274459682369542};
    static constexpr double bddata[48] = {
      0.20511993619098767,  0.09648410552837837, -2.2538795735050865,
      -0.04983109443768741, -0.9642251558357876, 0.19359255059179328,
      -0.39806578964162675, 0.26655194154076056, -0.44550440724595763,
      0.6706580677244948,   0.1094503631997486,  2.088335208353239,
      0.1428617396146469,   -0.5809206484565365, -3.1420312577464244,
      0.7798855526203933,   -0.1381973048177799, 0.08427369552387146,
      -0.46819882246485106, 0.0698482347780988,  0.24309332281112794,
      0.3841802990849566,   -1.4946422537617146, -0.5928606061328459,
      0.043976955387987425, -0.1226602268238547, -0.5881911245335547,
      0.7186280674937524,   -0.9224382403836392, -1.117883803285957,
      -0.1025638204862866,  0.5621700920192116,  -0.7417763253026113,
      0.8710856559160713,   0.15338641454990598, -0.9911471735749906,
      -1.2481676162786608,  2.1160878512755934,  -0.30475289132078964,
      0.4472989810307133,   0.43406926847671146, 0.2778024646943355,
      -0.11106692926404803, 0.18787583104988448, -0.09522217722177315,
      0.3930685232252409,   -0.9491225558721068, -2.5673833776578996};

    SquareMatrix<double> A(SquareDims<>{::math::row(4)});
    std::memcpy(A.data(), adata, 16 * 8UL);
    SquareMatrix<double> B(SquareDims<>{::math::row(4)});
    std::memcpy(B.data(), bdata, 16 * 8UL);
    expect(norm2(B - expm(A)) <= 1e-10);

    static_assert(utils::Compressible<Dual<double, 2>>);
    SquareMatrix<Dual<double, 2>> Ad(SquareDims<>{::math::row(4)});
    std::memcpy(Ad.data(), addata, 48 * 8UL);
    SquareMatrix<Dual<double, 2>> Bd(SquareDims<>{::math::row(4)});
    std::memcpy(Bd.data(), bddata, 48 * 8UL);

    expect(dualDeltaCmp(norm2(Bd - expm(Ad)), 1e-10));
    {
      auto x = Bd[3, 3] - 0.35;
      auto y = -1 * (0.35 - Bd[3, 3]);
      expect(approx(x.value(), y.value(), 1e-14));
      expect(approx(x.gradient()[0], y.gradient()[0], 1e-14));
      expect(approx(x.gradient()[1], y.gradient()[1], 1e-14));
    }
    expect(::math::smax(Bd[3, 3], 0.35) == ::math::smax(0.35, Bd[3, 3]));
  };

  return 0;
}

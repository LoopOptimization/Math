import boost.ut;
import Arena;
import Array;
import ArrayConcepts;
import CorePrint;
import Dual;
import ExprTemplates;
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

using namespace ::math;
using ::math::transpose;
using namespace boost::ut;

auto main() -> int {
  "DualTest BasicAssertions"_test = [] -> void {
    alloc::OwningArena arena;

    std::mt19937 gen(0);
    std::uniform_real_distribution<double> dist(-1, 1);
    SquareMatrix<double> A(SquareDims<>{::math::row(15)});
    Vector<double> x(length(15));
    for (auto &a : A) a = dist(gen);
    for (auto &xx : x) xx = dist(gen);
    SquareMatrix<double> B = A + A.t();
    const auto halfquadform = [&](const auto &y) {
      return 0.5 * ((y * B) * transpose(y));
    };
    Vector<double> g = x * B;
    auto f = halfquadform(x);

    auto [fx, gx] = gradient(&arena, x, halfquadform);
    auto [fxx, gxx, hxx] = hessian(&arena, x, halfquadform);
    expect(std::abs(fx - f) < 1e-10);
    expect(std::abs(fxx - f) < 1e-10);
    expect(norm2(g - gx) < 1e-10);
    expect(norm2(g - gxx) < 1e-10);
    utils::print("g = ");
    g.print();
    utils::print("\ngxx = ");
    gxx.print();
    utils::print('\n');
    utils::print("B = ");
    B.print();
    utils::print("\nhxx = ");
    hxx.print();
    utils::print('\n');
    for (std::ptrdiff_t i = 0; i < hxx.numRow(); ++i)
      for (std::ptrdiff_t j = i + 1; j < hxx.numCol(); ++j)
        hxx[i, j] = hxx[j, i];
    utils::print("hxx = ");
    hxx.print();
    utils::print('\n');
    expect(norm2(B - hxx) < 1e-10);
  };
  // NOLINTNEXTLINE(modernize-use-trailing-return-type)
  "IntDivDualTest BasicAssertions"_test = [] {
    Dual<double, 2> x{
      -0.1025638204862866,
      SVector<double, 2>{0.5621700920192116, -0.7417763253026113}};
    std::int32_t a = 4, b = 3;
    Dual<double, 2> y = a / x + b;
    expect(approx(y.value(), -36.00010726038452, 1e-14));
    expect(approx(y.gradient()[0], -213.76635331423668, 1e-14));
    expect(approx(y.gradient()[1], 282.061999181122, 1e-14));

    Dual<Dual<double, 8>, 2> w{
      Dual<double, 8>{
        0.474029877977747,
        SVector<double, 8>{0.3086698530598352, 0.2473835557435392,
                           0.3323869313428053, 0.5171334596793957,
                           0.38702317404526887, 0.00047786444101627357,
                           0.5631545736256198, 0.43922995203510906}},
      SVector<Dual<double, 8>, 2>{
        Dual<double, 8>{
          0.364700717714898,
          SVector<double, 8>{0.49478829029794746, 0.8045632402683945,
                             0.14233875934752005, 0.6248974091625752,
                             0.4100454368750559, 0.36093334233891017,
                             0.7299307580759404, 0.9599831981794166}},
        Dual<double, 8>{
          0.8935518055743592,
          SVector<double, 8>{0.1802316112143557, 0.7998299168331432,
                             0.6885713578218868, 0.16063226225861582,
                             0.8882724638859483, 0.3121292626973291,
                             0.8389228620948505, 0.8982533827333361}}}};
    Dual<Dual<double, 8>, 2> z = a / w + b;
    expect(approx(z.value().value(), 11.438286668900178, 1e-14));
    expect(approx(z.value().gradient()[0], -5.494684675315883, 1e-14));
    expect(approx(z.value().gradient()[1], -4.403716848906783, 1e-14));
    expect(approx(z.value().gradient()[2], -5.916876429038721, 1e-14));
    expect(approx(z.value().gradient()[3], -9.205580874924776, 1e-14));
    expect(approx(z.value().gradient()[4], -6.88946549958806, 1e-14));
    expect(approx(z.value().gradient()[5], -0.008506546379252392, 1e-14));
    expect(approx(z.value().gradient()[6], -10.024810569806137, 1e-14));
    expect(approx(z.value().gradient()[7], -7.818807254621021, 1e-14));
    expect(approx(z.gradient()[0].value(), -6.4920996489938965, 1e-14));
    expect(approx(z.gradient()[0].gradient()[0], -0.353004214108692, 1e-14));
    expect(approx(z.gradient()[0].gradient()[1], -7.546059942645349, 1e-14));
    expect(approx(z.gradient()[0].gradient()[2], 6.570646809055814, 1e-14));
    expect(approx(z.gradient()[0].gradient()[3], 3.040948459025321, 1e-14));
    expect(approx(z.gradient()[0].gradient()[4], 3.301701335355337, 1e-14));
    expect(approx(z.gradient()[0].gradient()[5], -6.4119467254828155, 1e-14));
    expect(approx(z.gradient()[0].gradient()[6], 2.431800795666499, 1e-14));
    expect(approx(z.gradient()[0].gradient()[7], -5.057833482826592, 1e-14));
    expect(approx(z.gradient()[1].value(), -15.906268020733835, 1e-14));
    expect(approx(z.gradient()[1].gradient()[0], 17.50675476102715, 1e-14));
    expect(approx(z.gradient()[1].gradient()[1], 2.3642057402333787, 1e-14));
    expect(approx(z.gradient()[1].gradient()[2], 10.049384954558796, 1e-14));
    expect(approx(z.gradient()[1].gradient()[3], 31.8458106722888, 1e-14));
    expect(approx(z.gradient()[1].gradient()[4], 10.161154827174284, 1e-14));
    expect(approx(z.gradient()[1].gradient()[5], -5.524196339292094, 1e-14));
    expect(approx(z.gradient()[1].gradient()[6], 22.859958982249076, 1e-14));
    expect(approx(z.gradient()[1].gradient()[7], 13.487122714859662, 1e-14));
  };

  return 0;
}

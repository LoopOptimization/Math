import boost.ut;
import Saturated;
import std;

using namespace boost::ut;

void testBasicAssertions() {

  constexpr int imax = std::numeric_limits<int>::max(),
                imin = std::numeric_limits<int>::min();
  expect(::math::add_sat(imax - 10, imax / 2) == imax);
  expect(::math::add_sat(imin + 10, imin / 2) == imin);
  expect(::math::sub_sat(imax - 10, imin / 2) == imax);
  expect(::math::sub_sat(imin + 10, imax / 2) == imin);
  expect(::math::mul_sat(imax - 10, imax / 2) == imax);
  expect(::math::mul_sat(imin + 10, imax / 2) == imin);

  constexpr auto iminlong = static_cast<long long>(imin),
                 imaxlong = static_cast<long long>(imax);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> idistrib(imin, imax);
  for (int i = 0; i < 1000; ++i) {
    int a = idistrib(gen), b = idistrib(gen);
    long long a64 = static_cast<long long>(a) + b,
              s64 = static_cast<long long>(a) - b,
              m64 = static_cast<long long>(a) * b;
    expect(::math::add_sat(a, b) == std::clamp(a64, iminlong, imaxlong));
    expect(::math::sub_sat(a, b) == std::clamp(s64, iminlong, imaxlong));
    expect(::math::mul_sat(a, b) == std::clamp(m64, iminlong, imaxlong));
  }

  constexpr unsigned int umax = std::numeric_limits<unsigned int>::max();
  expect(::math::add_sat(umax - 10, umax / 2) == umax);
  expect(::math::sub_sat(umax / 2, umax - 1) == 0);
  expect(::math::mul_sat(umax - 10, umax / 2) == umax);

  constexpr auto umaxlong = static_cast<long long>(umax);
  std::uniform_int_distribution<unsigned> udistrib(0, umax);
  for (int i = 0; i < 1000; ++i) {
    unsigned a = udistrib(gen), b = udistrib(gen);
    long long a64 = static_cast<long long>(a) + b,
              s64 = static_cast<long long>(a) - b;
    unsigned long long m64 = static_cast<unsigned long long>(a) * b;
    expect(::math::add_sat(a, b) == std::clamp(a64, (long long)0, umaxlong));
    expect(::math::sub_sat(a, b) == std::clamp(s64, (long long)0, umaxlong));
    expect(::math::mul_sat(a, b) == std::clamp(m64, (unsigned long long)0,
                                              (unsigned long long)umaxlong));
  }
}

int main() {
  "SaturatedArithmetic BasicAssertions"_test = [] {
    testBasicAssertions();
  };
  return 0;
}

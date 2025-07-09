import boost.ut;

import CorePrint;
import Int8;
import std;

using numbers::i8, numbers::u8, numbers::Flag8;
using namespace boost::ut;
static_assert(utils::Printable<numbers::i8>);
static_assert(utils::Printable<numbers::u8>);

int main() {
  "Int8Test BasicAssertions"_test = [] {
    for (std::uint8_t x = 0; x < std::numeric_limits<std::uint8_t>::max();
         ++x) {
      auto fx = static_cast<Flag8>(x);
      if (x) expect(bool(fx));
      else expect(!bool(fx));
      auto ux = static_cast<u8>(x);
      for (std::uint8_t y = 0; y < std::numeric_limits<std::uint8_t>::max();
           ++y) {
        auto uy = static_cast<u8>(y);
        utils::print(uy);
        expect((ux <=> uy) == (x <=> y));
        {
          u8 z = static_cast<u8>(x + y);
          expect(z == z);
          expect(!(z != z));
        }
        bool b = (ux > uy) == (x > y);
        expect(b);
        expect((ux == uy) == (x == y));
        expect((ux == y) == (x == y));
        expect((x == uy) == (x == y));
        expect((ux != uy) == (x != y));
        expect((ux != y) == (x != y));
        expect((x != uy) == (x != y));
        expect((ux > uy) == (x > y));
        expect((ux > y) == (x > y));
        expect((x > uy) == (x > y));
        expect((ux < uy) == (x < y));
        expect((ux < y) == (x < y));
        expect((x < uy) == (x < y));
        expect((ux >= uy) == (x >= y));
        expect((ux >= y) == (x >= y));
        expect((x >= uy) == (x >= y));
        expect((ux <= uy) == (x <= y));
        expect((ux <= y) == (x <= y));
        expect((x <= uy) == (x <= y));

        expect((ux > uy) == (x > y));

        {
          int z = x + y;
          if ((z >= 0) && (z <= std::numeric_limits<std::uint8_t>::max())) {
            u8 uz = ux;
            uz += uy;
            expect(uz == z);
            expect((ux + uy) == z);
            expect((ux + y) == z);
            expect((x + uy) == z);
          }
        }
        {
          int z = x * y;
          if ((z >= 0) && (z <= std::numeric_limits<std::uint8_t>::max())) {
            u8 uz = ux;
            uz *= uy;
            expect(uz == z);
            expect((ux * uy) == z);
            expect((ux * y) == z);
            expect((x * uy) == z);
          }
        }
        {
          int z = x - y;
          if ((z >= 0) && (z <= std::numeric_limits<std::uint8_t>::max())) {
            u8 uz = ux;
            uz -= uy;
            expect(uz == z);
            expect((ux - uy) == z);
            expect((ux - y) == z);
            expect((x - uy) == z);
          }
        }
        if (y) {
          int z = x / y;
          if ((z >= 0) && (z <= std::numeric_limits<std::uint8_t>::max())) {
            u8 uz = ux;
            uz /= uy;
            expect(uz == z);
            expect((ux / uy) == z);
            expect((ux / y) == z);
            expect((x / uy) == z);
          }
        }
        auto fy = static_cast<Flag8>(y);
        expect((fx & fy) == static_cast<Flag8>(x & y));
        expect((fx & y) == static_cast<Flag8>(x & y));
        expect((x & fy) == static_cast<Flag8>(x & y));
        expect((fx | fy) == static_cast<Flag8>(x | y));
        expect((fx | y) == static_cast<Flag8>(x | y));
        expect((x | fy) == static_cast<Flag8>(x | y));
        expect((fx ^ fy) == static_cast<Flag8>(x ^ y));
        expect((fx ^ y) == static_cast<Flag8>(x ^ y));
        expect((x ^ fy) == static_cast<Flag8>(x ^ y));

        expect((fx & fy) == (x & y));
        expect((fx & y) == (x & y));
        expect((x & fy) == (x & y));
        expect((fx | fy) == (x | y));
        expect((fx | y) == (x | y));
        expect((x | fy) == (x | y));
        expect((fx ^ fy) == (x ^ y));
        expect((fx ^ y) == (x ^ y));
        expect((x ^ fy) == (x ^ y));

        expect((x & y) == (fx & fy));
        expect((x & y) == (fx & y));
        expect((x & y) == (x & fy));
        expect((x | y) == (fx | fy));
        expect((x | y) == (fx | y));
        expect((x | y) == (x | fy));
        expect((x ^ y) == (fx ^ fy));
        expect((x ^ y) == (fx ^ y));
        expect((x ^ y) == (x ^ fy));
        expect((fx >> 1) == (x >> 1));
        if (y) {
          expect(static_cast<bool>(fy));
          expect(!(!fy));
        } else {
          expect(!bool(fy));
          expect(!fy);
        }
      }
    }
    for (std::int8_t x = std::numeric_limits<std::int8_t>::min();
         x < std::numeric_limits<std::int8_t>::max(); ++x) {
      for (std::int8_t y = std::numeric_limits<std::int8_t>::min();
           y < std::numeric_limits<std::int8_t>::max(); ++y) {
        i8 ix = static_cast<i8>(x);
        i8 iy = static_cast<i8>(y);

        {
          int z = x + y;
          if ((z >= std::numeric_limits<std::int8_t>::min()) &&
              (z <= std::numeric_limits<std::int8_t>::max())) {
            i8 iz = ix;
            iz += iy;
            expect(iz == z);
            expect((ix + iy) == z);
            expect((ix + y) == z);
            expect((x + iy) == z);
          }
        }
        {
          int z = x * y;
          if ((z >= std::numeric_limits<std::int8_t>::min()) &&
              (z <= std::numeric_limits<std::int8_t>::max())) {
            i8 iz = ix;
            iz *= iy;
            expect(iz == z);
            expect((ix * iy) == z);
            expect((ix * y) == z);
            expect((x * iy) == z);
          }
        }
        {
          int z = x - y;
          if ((z >= std::numeric_limits<std::int8_t>::min()) &&
              (z <= std::numeric_limits<std::int8_t>::max())) {
            i8 iz = ix;
            iz -= iy;
            expect(iz == z);
            expect((ix - iy) == z);
            expect((ix - y) == z);
            expect((x - iy) == z);
          }
        }
        if (y) {
          int z = x / y;
          if ((z >= std::numeric_limits<std::int8_t>::min()) &&
              (z <= std::numeric_limits<std::int8_t>::max())) {
            i8 iz = ix;
            iz /= iy;
            expect(iz == z);
            expect((ix / iy) == z);
            expect((ix / y) == z);
            expect((x / iy) == z);
          }
        }
      }
    }
  };
  return 0;
}

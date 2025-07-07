import boost.ut;
import BaseUtils;
import std;

using namespace boost::ut;

int main() {
  "BitTest BasicAssertions"_test = [] {
    for (int i = 0; i < 63; ++i) {
      auto e2 = double(std::uint64_t(1) << i);
      expect(bit::exp2unchecked(i) == e2);
      expect(bit::exp2unchecked(-i) == 1.0 / e2);
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<double> dist(15.6);
    for (int i = 0; i < 10000; ++i) {
      double x = dist(gen);
      expect(bit::next_pow2(x) == std::exp2(std::ceil(std::log2(x))));
    }
  };
  return 0;
}

import Testing;
import CorePrint;
import Permutation;
import std;
import TinyVector;

using namespace testing;
using Order = utils::Order<>;

auto main() -> int {
  "PermutationTest.BasicAssertions"_test = [] -> void {
    int count = 0;
    std::array<Order, 6> res{
      Order::create({0, 1, 2}), Order::create({1, 0, 2}),
      Order::create({2, 0, 1}), Order::create({0, 2, 1}),
      Order::create({1, 2, 0}), Order::create({2, 1, 0})};
    for (auto p : utils::Permutations(3)) {
      utils::print("Perm: ", p, '\n');
      expect(p == res[count++]);
    };

    expect(eq(count, 6));
  };
  "LoopPermutationTest.BasicAssertions"_test = [] -> void {
    using utils::LoopSet;
    containers::TinyVector<LoopSet, 14, std::int16_t> cmpts;
    // cmpts are [1], [0,3], [2, 4]
    cmpts.push_back(LoopSet::fromMask(0x02)); // 0x00000010
    cmpts.push_back(LoopSet::fromMask(0x09)); // 0x00001001
    cmpts.push_back(LoopSet::fromMask(0x14)); // 0x00010100
    int count = 0;
    std::array<Order, 6> res{
      Order::create({1, 3, 0, 4, 2}), Order::create({1, 0, 3, 4, 2}),
      Order::create({1, 0, 3, 2, 4}), Order::create({1, 3, 0, 2, 4})};
    for (auto p : utils::LoopPermutations(cmpts, 5)) {
      utils::print("Perm: ", p, '\n');
      expect(res[count++] == p);
    };

    expect(eq(count, 4));
  };
  return 0;
}

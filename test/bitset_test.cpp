#include <algorithm>

import boost.ut;

import Arena;
import ArrayParse;
import BitSet;
import CorePrint;
import ManagedArray;
import std;

#define STRINGIZE_DETAIL(x) #x
#define STRINGIZE(x) STRINGIZE_DETAIL(x)
using containers::BitSet, containers::BitSets, ::math::Vector,
  utils::operator""_mat;

auto main() -> int {
  using namespace boost::ut;

  "BitSet BasicAssertions"_test = [] -> void {
    BitSet bs(1000);
    std::ptrdiff_t count = 0;
    for (std::ptrdiff_t _ : bs) ++count;
    expect(fatal(count == 0));
    bs[4] = true;
    bs[10] = true;
    bs[200] = true;
    bs[117] = true;
    bs[87] = true;
    bs[991] = true;
    expect(fatal(bs.findFirstZero() == 0));
    bs[0] = true;
    expect(fatal(bs.findFirstZero() == 1));
    bs.print();
    utils::print('\n');
    expect(bs.isSubSet(bs));
    // expect(std::ranges::begin(bs), bs.begin());
    // expect(std::ranges::end(bs), bs.end());
    std::array bsc{0, 4, 10, 87, 117, 200, 991};
    std::ptrdiff_t j = 0;
    for (std::ptrdiff_t J : bs) {
      expect(eq(bs.getIterationIdx(J), j));
      expect(eq(J, bsc[j++]));
      expect(bs[J]);
      utils::println("We get: ", J);
    }
    j = 0;
    for (std::ptrdiff_t i : bs) {
      expect(i == bsc[j++]);
      expect(bs[i]);
      utils::println("We get: ", i);
    }
    expect(j == std::ssize(bsc));
    expect(j == std::ptrdiff_t(bs.size()));
    utils::println("About to create empty!");
    BitSet empty;
    std::ptrdiff_t c = 0, d = 0;
    utils::println("About to iterate empty!");
    for (std::ptrdiff_t b : empty) {
      ++c;
      d += b;
    }
    utils::println("Iterated empty!");
    expect(!c);
    expect(!d);
    utils::println("we made it to the end!");
  };
  "BitSet Insert"_test = [] -> void {
    BitSet<std::array<std::uint64_t, 2>> bs;
    std::ptrdiff_t count = 0;
    for (std::ptrdiff_t _ : bs) ++count;
    expect(fatal(count == 0));
    bs.insert(1);
    bs.insert(5);
    bs.insert(6);
    bs.insert(8);
    expect(bs.data_[0] == 354);
    expect(bs.data_[1] == 0);
    bs.insert(5);
    expect(bs.data_[0] == 354);

    std::array args{144, 146};
    containers::BitSet<::math::Vector<std::uint64_t, 8>> aug_vars{};
    for (int x : args) aug_vars.uncheckedInsert(x);
    expect(eq(aug_vars.front(), std::min(args[0], args[1])));
    for (std::ptrdiff_t i = 0; std::ptrdiff_t x : aug_vars)
      expect(eq(x, args[i++]));
  };
  "BitSet DynSize"_test = [] -> void {
    BitSet bs, bsd{BitSet<>::dense(11)};
    std::ptrdiff_t count = 0;
    for (auto chunk : bsd.chunks()) {
      expect(eq(chunk.begin(), 0));
      expect(eq(chunk.end(), 11));
      ++count;
    }
    expect(eq(count, 1));
    count = 0;
    for (std::ptrdiff_t _ : bs) ++count;
    expect(fatal(count == 0));
    for (std::ptrdiff_t _ : bsd) ++count;
    expect(fatal(count == 11));
    expect(bs.data_.size() == 0);
    expect(bs.isSubSet(bsd));
    bs[4] = true;
    bs[10] = true;
    expect(bsd.isSuperSet(bs));
    expect(!bsd.isSubSet(bs));
    expect(bs.data_.size() == 1);
    expect(bs.data_.front() == 1040);
    for (std::ptrdiff_t i = 0; i < 11; ++i)
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-else"
#endif
      if (!bs.contains(i)) expect(bsd.remove(i));
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    expect(bs == bsd);
    Vector<std::size_t> sv;
    for (auto i : bs) sv.push_back(i);
    expect(sv.size() == 2);
    expect(sv[0] == 4);
    expect(sv[1] == 10);
  };
  "BitSet FixedSize"_test = [] -> void {
    BitSet<std::array<std::uint64_t, 2>> bs;
    bs[4] = true;
    bs[10] = true;
    expect(bs.data_[0] == 0b10000010000);
    expect(bs.data_[1] == 0);
    expect(eq(bs.maxValue(), 10));
    Vector<std::size_t> sv;
    for (auto i : bs) sv.push_back(i);
    expect(sv.size() == 2);
    expect(sv[0] == 4);
    expect(sv[1] == 10);
  };
  "BitSet FixedSizeSmall"_test = [] -> void {
    using SB = BitSet<std::array<std::uint16_t, 1>>;
    static_assert(sizeof(SB) == 2);
    SB bs;
    std::ptrdiff_t count = 0;
    for (std::ptrdiff_t _ : bs) ++count;
    expect(fatal(count == 0));
    bs[4] = true;
    bs[10] = true;
    bs[7] = true;
    bs.insert(5);
    expect(bs.data_[0] == 1200);
    Vector<std::size_t> sv;
    for (auto i : bs) sv.push_back(i);
    expect(sv.size() == 4);
    expect(sv[0] == 4);
    expect(sv[1] == 5);
    expect(sv[2] == 7);
    expect(sv[3] == 10);
    expect(SB::fromMask(1200).data_[0] == 1200);
  };
  "BitSet IterTest"_test = [] -> void {
    ::math::Vector<std::uint64_t, 1> data{};
    data.push_back(9223372036854775808U);
    data.push_back(1732766553700568065U);
    data.push_back(1891655728U);
    BitSet<> bs{data};
    std::ptrdiff_t sz = bs.size(), i = 0;
    for (std::ptrdiff_t a : bs) {
      expect(eq(bs.getIterationIdx(a), i++));
      utils::print(a, " ");
    }
    utils::print('\n');
    expect(fatal(eq(i, sz)));
  };
  "BitSet EmptyIntersection"_test = [] -> void {
    BitSet bs1, bs2;
    expect(bs1.emptyIntersection(bs2));

    bs1[4] = true;
    bs1[10] = true;
    expect(bs1.emptyIntersection(bs2));
    expect(bs2.emptyIntersection(bs1));

    bs2[5] = true;
    bs2[11] = true;
    expect(bs1.emptyIntersection(bs2));
    expect(bs2.emptyIntersection(bs1));

    bs2[4] = true;
    expect(!bs1.emptyIntersection(bs2));
    expect(!bs2.emptyIntersection(bs1));

    BitSet<std::array<std::uint64_t, 4>> bs3, bs4;
    bs3[100] = true;
    bs4[200] = true;
    expect(bs3.emptyIntersection(bs4));

    bs4[100] = true;
    expect(!bs3.emptyIntersection(bs4));
  };

  "BitSet ExpressionTemplate"_test = [] -> void {
    // Test basic expression template assignment with PtrVector
    auto v{"[5 0 -3 7 0 -1 8 0 2]"_mat};

    BitSet<std::array<std::uint16_t, 1>> bs{};

    // Test bs << (x != 0) pattern
    bs << (v != 0);

    // Verify bits are set correctly for non-zero elements
    expect(bs[0]);  // 5 != 0
    expect(!bs[1]); // 0 == 0
    expect(bs[2]);  // -3 != 0
    expect(bs[3]);  // 7 != 0
    expect(!bs[4]); // 0 == 0
    expect(bs[5]);  // -1 != 0
    expect(bs[6]);  // 8 != 0
    expect(!bs[7]); // 0 == 0
    expect(bs[8]);  // 2 != 0

    expect(bs.size() == 6); // 6 non-zero elements
  };

  "BitSet ExpressionTemplate GreaterThan"_test = [] -> void {
    // Test with different comparison operators
    auto data{"[5 -2 3 -7 0 1 8 -1 2]"_mat};
    ::math::PtrVector<std::int64_t> vec{data};

    BitSet<std::array<std::uint16_t, 1>> bs{};

    // Test bs << (x > 0) pattern
    bs << (vec > 0);

    // Verify bits are set correctly for positive elements
    expect(bs[0]);  // 5 > 0
    expect(!bs[1]); // -2 <= 0
    expect(bs[2]);  // 3 > 0
    expect(!bs[3]); // -7 <= 0
    expect(!bs[4]); // 0 <= 0
    expect(bs[5]);  // 1 > 0
    expect(bs[6]);  // 8 > 0
    expect(!bs[7]); // -1 <= 0
    expect(bs[8]);  // 2 > 0

    expect(bs.size() == 5); // 5 positive elements
  };

  "BitSet ExpressionTemplate Complex"_test = [] -> void {
    // Test with more complex expressions
    auto data{"[10 5 -3 15 2 8 -1 20 0]"_mat};
    ::math::PtrVector<std::int64_t> vec{data};

    BitSet<std::array<std::uint16_t, 1>> bs1{}, bs2{};

    // Test bs << (x >= 5)
    bs1 << (vec >= 5);
    expect(bs1[0]);  // 10 >= 5
    expect(bs1[1]);  // 5 >= 5
    expect(!bs1[2]); // -3 < 5
    expect(bs1[3]);  // 15 >= 5
    expect(!bs1[4]); // 2 < 5
    expect(bs1[5]);  // 8 >= 5
    expect(!bs1[6]); // -1 < 5
    expect(bs1[7]);  // 20 >= 5
    expect(!bs1[8]); // 0 < 5

    // Test bs << (x < 0)
    bs2 << (vec < 0);
    expect(!bs2[0]); // 10 >= 0
    expect(!bs2[1]); // 5 >= 0
    expect(bs2[2]);  // -3 < 0
    expect(!bs2[3]); // 15 >= 0
    expect(!bs2[4]); // 2 >= 0
    expect(!bs2[5]); // 8 >= 0
    expect(bs2[6]);  // -1 < 0
    expect(!bs2[7]); // 20 >= 0
    expect(!bs2[8]); // 0 >= 0

    expect(!bs1.isSubSet(bs2));
    expect(!bs1.isSuperSet(bs2));

    expect(bs1.size() == 5); // 5 elements >= 5
    expect(bs2.size() == 2); // 2 elements < 0
  };

  "BitSet ExpressionTemplate DynamicSize"_test = [] -> void {
    // Test with dynamic size BitSet
    auto vec{"[1 0 3 0 5 6 0 8]"_mat};

    BitSet<> bs{}; // Dynamic size BitSet

    // Test assignment with dynamic resize
    bs << (vec != 0);

    expect(bs[0]);  // 1 != 0
    expect(!bs[1]); // 0 == 0
    expect(bs[2]);  // 3 != 0
    expect(!bs[3]); // 0 == 0
    expect(bs[4]);  // 5 != 0
    expect(bs[5]);  // 6 != 0
    expect(!bs[6]); // 0 == 0
    expect(bs[7]);  // 8 != 0

    expect(bs.size() == 5); // 5 non-zero elements
  };

  "BitSet CompoundAssignment OR"_test = [] -> void {
    // Test |= operator with expression templates
    auto vec1{"[1 0 1 0 1 0]"_mat};
    auto vec2{"[0 1 1 1 0 0]"_mat};

    BitSet<std::array<std::uint16_t, 1>> bs{};

    // Start with vec1 != 0: [1,0,1,0,1,0] -> bits 0,2,4 set
    bs << (vec1 != 0);
    expect(bs[0]);
    expect(!bs[1]);
    expect(bs[2]);
    expect(!bs[3]);
    expect(bs[4]);
    expect(!bs[5]);
    expect(bs.size() == 3);

    auto bs2{bs};
    // OR with vec2 != 0: [0,1,1,1,0,0] -> bits 1,2,3 set
    // Result should be: bits 0,1,2,3,4 set
    bs |= (vec2 != 0);
    expect(bs[0]);
    expect(bs[1]);
    expect(bs[2]);
    expect(bs[3]);
    expect(bs[4]);
    expect(!bs[5]);
    expect(bs.size() == 5);
    expect(!bs.isSubSet(bs2));
    expect(bs2.isSubSet(bs));
    expect(bs.isSuperSet(bs2));
    expect(!bs2.isSuperSet(bs));
  };

  "BitSet CompoundAssignment AND"_test = [] -> void {
    // Test &= operator with expression templates
    auto vec1{"[1 1 1 1 1 0]"_mat};
    auto vec2{"[1 0 1 0 1 1]"_mat};

    BitSet<std::array<std::uint16_t, 1>> bs{};

    // Start with vec1 != 0: [1,1,1,1,1,0] -> bits 0,1,2,3,4 set
    bs << (vec1 != 0);
    expect(bs[0]);
    expect(bs[1]);
    expect(bs[2]);
    expect(bs[3]);
    expect(bs[4]);
    expect(!bs[5]);
    expect(bs.size() == 5);

    // AND with vec2 != 0: [1,0,1,0,1,1] -> bits 0,2,4,5 set
    // Result should be: bits 0,2,4 set (intersection)
    bs &= (vec2 != 0);
    expect(bs[0]);
    expect(!bs[1]);
    expect(bs[2]);
    expect(!bs[3]);
    expect(bs[4]);
    expect(!bs[5]);
    expect(bs.size() == 3);
  };

  "BitSet CompoundAssignment XOR"_test = [] -> void {
    // Test ^= operator with expression templates
    auto vec1{"[1 1 0 1 0 1]"_mat};
    auto vec2{"[1 0 1 1 1 0]"_mat};

    BitSet<std::array<std::uint16_t, 1>> bs{};

    // Start with vec1 != 0: [1,1,0,1,0,1] -> bits 0,1,3,5 set
    bs << (vec1 != 0);
    expect(bs[0]);
    expect(bs[1]);
    expect(!bs[2]);
    expect(bs[3]);
    expect(!bs[4]);
    expect(bs[5]);
    expect(bs.size() == 4);

    // XOR with vec2 != 0: [1,0,1,1,1,0] -> bits 0,2,3,4 set
    // Result should be: bits 1,2,4,5 set (symmetric difference)
    bs ^= (vec2 != 0);
    expect(!bs[0]);
    expect(bs[1]);
    expect(bs[2]);
    expect(!bs[3]);
    expect(bs[4]);
    expect(bs[5]);
    expect(bs.size() == 4);
  };

  "BitSet CompoundAssignment Complex"_test = [] -> void {
    // Test compound assignment with more complex expressions
    auto data{"[5 -2 3 -7 0 1 8 -1]"_mat};
    ::math::PtrVector<std::int64_t> vec{data};

    BitSet<std::array<std::uint16_t, 1>> bs{};

    // Start with positive elements: vec > 0 -> bits 0,2,5,6 set
    bs << (vec > 0);
    expect(bs[0]);
    expect(!bs[1]);
    expect(bs[2]);
    expect(!bs[3]);
    expect(!bs[4]);
    expect(bs[5]);
    expect(bs[6]);
    expect(!bs[7]);
    expect(bs.size() == 4);

    // OR with elements >= 3: [5>=3, -2>=3, 3>=3, -7>=3, 0>=3, 1>=3, 8>=3,
    // -1>=3]
    // -> bits 0,2,6 set
    // Result should be: bits 0,2,5,6 set (union)
    bs |= (vec >= 3);
    expect(bs[0]);
    expect(!bs[1]);
    expect(bs[2]);
    expect(!bs[3]);
    expect(!bs[4]);
    expect(bs[5]);
    expect(bs[6]);
    expect(!bs[7]);
    expect(bs.size() == 4);

    // AND with elements != -2: all except bit 1 -> bits 0,2,5,6 set (no change)
    bs &= (vec != -2);
    expect(bs[0]);
    expect(!bs[1]);
    expect(bs[2]);
    expect(!bs[3]);
    expect(!bs[4]);
    expect(bs[5]);
    expect(bs[6]);
    expect(!bs[7]);
    expect(bs.size() == 4);

    // XOR with elements < 0: bits 1,3,7 set
    // Result should be: bits 0,1,2,3,5,6,7 set
    bs ^= (vec < 0);
    expect(bs[0]);
    expect(bs[1]);
    expect(bs[2]);
    expect(bs[3]);
    expect(!bs[4]);
    expect(bs[5]);
    expect(bs[6]);
    expect(bs[7]);
    expect(bs.size() == 7);
  };

  "BitSet CompoundAssignment DynamicSize"_test = [] -> void {
    // Test compound assignment with dynamic size BitSet
    auto vec1{"[1 0 1 0 1]"_mat};
    auto vec2{"[0 1 1 1 0]"_mat};

    BitSet<> bs{}; // Dynamic size BitSet

    // Start with vec1
    bs << (vec1 != 0);
    expect(bs.size() == 3);

    // OR with vec2
    bs |= (vec2 != 0);
    expect(bs[0]);
    expect(bs[1]);
    expect(bs[2]);
    expect(bs[3]);
    expect(bs[4]);
    expect(bs.size() == 5);

    // AND to keep only common elements with vec1 again
    bs &= (vec1 != 0);
    expect(bs[0]);
    expect(!bs[1]);
    expect(bs[2]);
    expect(!bs[3]);
    expect(bs[4]);
    expect(bs.size() == 3);
  };

  "BitSets BasicConstruction"_test = [] -> void {
    alloc::OwningArena<> alloc;

    // Test empty factory method with small collection
    BitSets bs_empty = BitSets::empty(&alloc, 3, 100);
    expect(bs_empty.size() == 3);
    expect(bs_empty.maximumSize() == 128);    // Rounded up to 64-bit boundaries
    expect(bs_empty.maximumNumChunks() == 2); // (100 + 63) / 64 = 2

    // Test undef factory method with larger collection
    BitSets bs_undef = BitSets::undef(&alloc, 5, 200);
    expect(bs_undef.size() == 5);
    expect(bs_undef.maximumSize() == 256);    // Rounded up to 64-bit boundaries
    expect(bs_undef.maximumNumChunks() == 4); // (200 + 63) / 64 = 4

    // Test accessing individual BitSets
    auto bs0 = bs_empty[0];
    auto bs1 = bs_empty[1];
    auto bs2 = bs_empty[2];

    // All should be empty initially
    expect(bs0.empty());
    expect(bs1.empty());
    expect(bs2.empty());
    expect(bs0.size() == 0);
    expect(bs1.size() == 0);
    expect(bs2.size() == 0);
  };

  "BitSets IndividualOperations"_test = [] -> void {
    alloc::OwningArena<> alloc;
    BitSets bs_coll = BitSets::empty(&alloc, 4, 64);

    // Get individual BitSets
    auto bs0 = bs_coll[0];
    auto bs1 = bs_coll[1];
    auto bs2 = bs_coll[2];
    auto bs3 = bs_coll[3];

    // Test setting bits in different BitSets
    bs0[5] = true;
    bs0[10] = true;
    bs1[3] = true;
    bs1[15] = true;
    bs2[0] = true;
    bs2[63] = true;

    // Verify each BitSet has correct bits set
    expect(bs0[5] && bs0[10]);
    expect(!bs0[3] && !bs0[15] && !bs0[0] && !bs0[63]);
    expect(eq(bs0.size(), 2));

    expect(bs1[3] && bs1[15]);
    expect(!bs1[5] && !bs1[10] && !bs1[0] && !bs1[63]);
    expect(eq(bs1.size(), 2));

    expect(bs2[0] && bs2[63]);
    expect(!bs2[3] && !bs2[5] && !bs2[10] && !bs2[15]);
    expect(eq(bs2.size(), 2));

    expect(bs3.empty());
    expect(eq(bs3.size(), 0));

    // Test contains method
    expect(bs0.contains(5));
    expect(bs0.contains(10));
    expect(!bs0.contains(3));
    expect(!bs1.contains(5));
    expect(bs1.contains(3));

    // Test insert method
    expect(!bs3.insert(42)); // Should return false (wasn't there)
    expect(bs3.insert(42));  // Should return true (now it's there)
    expect(bs3.contains(42));
    expect(eq(bs3.size(), 1));

    // Test remove method
    expect(bs0.remove(5));  // Should return true (was there)
    expect(!bs0.remove(5)); // Should return false (not there anymore)
    expect(!bs0.contains(5));
    expect(eq(bs0.size(), 1));
  };

  "BitSets MemoryLayout"_test = [] -> void {
    alloc::OwningArena<> alloc;

    // Test edge case: 1 element per set
    BitSets bs_tiny = BitSets::empty(&alloc, 2, 1);
    expect(eq(bs_tiny.size(), 2));
    expect(bs_tiny.maximumSize() == 64); // Minimum 1 chunk = 64 bits
    expect(bs_tiny.maximumNumChunks() == 1);

    // Test exactly 64 elements (1 chunk boundary)
    BitSets bs_64 = BitSets::empty(&alloc, 3, 64);
    expect(bs_64.maximumSize() == 64);
    expect(bs_64.maximumNumChunks() == 1);

    // Test 65 elements (crosses chunk boundary)
    BitSets bs_65 = BitSets::empty(&alloc, 3, 65);
    expect(bs_65.maximumSize() == 128);
    expect(bs_65.maximumNumChunks() == 2);

    // Test large collection
    BitSets bs_large = BitSets::empty(&alloc, 10, 1000);
    expect(eq(bs_large.size(), 10));
    expect(bs_large.maximumSize() == 1024); // (1000 + 63) / 64 * 64 = 16 * 64
    expect(bs_large.maximumNumChunks() == 16);

    // Verify we can access all sets
    for (std::ptrdiff_t i = 0; i < 10; ++i) {
      auto bs = bs_large[i];
      expect(bs.empty());
      bs[999] = true; // Set the highest possible bit
      expect(bs.contains(999));
      expect(bs.size() == 1);
    }
  };

  "BitSets Independence"_test = [] -> void {
    alloc::OwningArena<> alloc;
    BitSets bs_coll = BitSets::empty(&alloc, 3, 100);

    auto bs0 = bs_coll[0];
    auto bs1 = bs_coll[1];
    auto bs2 = bs_coll[2];

    // Set the same element indices in different BitSets
    bs0[25] = true;
    bs1[25] = true;
    bs2[25] = true;

    bs0[50] = true;
    bs1[75] = true;

    // Verify independence - each BitSet should only have its own bits
    expect(bs0.contains(25) && bs0.contains(50));
    expect(!bs0.contains(75));
    expect(bs0.size() == 2);

    expect(bs1.contains(25) && bs1.contains(75));
    expect(!bs1.contains(50));
    expect(bs1.size() == 2);

    expect(bs2.contains(25));
    expect(!bs2.contains(50) && !bs2.contains(75));
    expect(bs2.size() == 1);

    // Test operations on one don't affect others
    bs0.clear();
    expect(bs0.empty());
    expect(bs0.size() == 0);

    // Other BitSets should be unaffected
    expect(bs1.size() == 2);
    expect(bs2.size() == 1);
    expect(bs1.contains(25) && bs1.contains(75));
    expect(bs2.contains(25));

    // Test BitSet operations
    BitSet<> temp_bs;
    temp_bs[25] = true;
    temp_bs[99] = true;

    bs1 |= temp_bs;            // Union with temp_bs
    expect(eq(bs1.size(), 3)); // Should now have 25, 75, 99
    expect(bs1.contains(25) && bs1.contains(75) && bs1.contains(99));

    // bs2 should still be unaffected
    expect(bs2.size() == 1);
    expect(bs2.contains(25) && !bs2.contains(99));
  };

  "BitSets ExpressionTemplates"_test = [] -> void {
    alloc::OwningArena<> alloc;
    BitSets bs_coll = BitSets::empty(&alloc, 2, 50);

    auto bs0 = bs_coll[0];
    auto bs1 = bs_coll[1];

    // Test expression template assignment with array data
    auto data{"[5 0 -3 7 0 -1 8 0 2]"_mat};

    // Use expression templates on individual BitSets
    bs0 << (data != 0);
    bs1 << (data > 0);

    // Verify bs0 has bits set for non-zero elements
    expect(bs0[0] && !bs0[1] && bs0[2] && bs0[3]);
    expect(!bs0[4] && bs0[5] && bs0[6] && !bs0[7] && bs0[8]);
    expect(eq(bs0.size(), 6));

    // Verify bs1 has bits set for positive elements only
    expect(bs1[0] && !bs1[1] && !bs1[2] && bs1[3]);
    expect(!bs1[4] && !bs1[5] && bs1[6] && !bs1[7] && bs1[8]);
    expect(eq(bs1.size(), 4));

    // Test compound assignment operators
    auto data2{"[1 1 0 0 1 1 0 1 0]"_mat};

    bs0 &= (data2 != 0); // Intersection
    // bs0 had: [1,0,1,1,0,1,1,0,1] (data != 0)
    // data2:   [1,1,0,0,1,1,0,1,0] (data2 != 0)
    // Result:  [1,0,0,0,0,1,0,0,0]
    expect(bs0[0] && !bs0[1] && !bs0[2] && !bs0[3]);
    expect(!bs0[4] && bs0[5] && !bs0[6] && !bs0[7] && !bs0[8]);
    expect(bs0.size() == 2);
  };

  "BitSets Iterator Basic"_test = [] -> void {
    alloc::OwningArena<> alloc;

    // Test with non-empty collection
    BitSets bs_coll = BitSets::empty(&alloc, 4, 50);

    auto begin_it = bs_coll.begin();
    auto end_it = bs_coll.end();

    // Test basic iterator properties
    expect(begin_it != end_it);
    expect(!(begin_it == end_it));

    // Test iterator dereferencing
    auto first_bitset = *begin_it;
    expect(first_bitset.empty());
    expect(eq(first_bitset.size(), 0));

    // Test that we can modify through dereferenced iterator
    first_bitset[10] = true;
    expect(first_bitset.contains(10));
    expect(eq(first_bitset.size(), 1));

    // Test iterator increment
    auto second_it = begin_it;
    ++second_it;
    expect(second_it != begin_it);
    expect(second_it != end_it);

    auto second_bitset = *second_it;
    expect(second_bitset.empty()); // Should be different from first_bitset

    // Test postfix increment
    auto third_it = second_it++;
    expect(third_it != second_it);
    expect(third_it == begin_it + 1); // third_it should be the old position

    // Test iterator comparison and ordering
    expect(begin_it < second_it);
    expect(second_it < end_it);
    expect(begin_it <= second_it);
    expect(second_it >= begin_it);

    // Test empty collection
    BitSets empty_coll = BitSets::empty(&alloc, 0, 10);
    expect(empty_coll.begin() == empty_coll.end());
  };

  "BitSets Iterator RangeFor"_test = [] -> void {
    alloc::OwningArena<> alloc;
    BitSets bs_coll = BitSets::empty(&alloc, 5, 100);

    // Test range-based for loop
    std::ptrdiff_t count = 0;
    for (auto bitset : bs_coll) {
      expect(bitset.empty());
      expect(eq(bitset.size(), 0));

      // Set different bits in each BitSet
      bitset[count * 10] = true;
      bitset[count * 10 + 5] = true;

      expect(eq(bitset.size(), 2));
      expect(bitset.contains(count * 10));
      expect(bitset.contains(count * 10 + 5));

      ++count;
    }

    expect(eq(count, 5)); // Should have iterated over all 5 BitSets

    // Verify each BitSet maintains its data independently
    count = 0;
    for (auto bitset : bs_coll) {
      expect(eq(bitset.size(), 2));
      expect(bitset.contains(count * 10));
      expect(bitset.contains(count * 10 + 5));
      ++count;
    }
  };

  "BitSets Iterator STL"_test = [] -> void {
    alloc::OwningArena<> alloc;
    BitSets bs_coll = BitSets::empty(&alloc, 7, 64);

    // Test std::distance
    auto distance = std::distance(bs_coll.begin(), bs_coll.end());
    expect(eq(distance, 7));

    // Test with std::for_each (or ranges equivalent)
    std::ptrdiff_t counter = 0;
    std::ranges::for_each(bs_coll, [&counter](auto bitset) {
      bitset[counter * 5] = true;
      expect(bitset.contains(counter * 5));
      ++counter;
    });
    expect(eq(counter, 7));

    // Test iterator satisfies bidirectional iterator requirements
    auto it = bs_coll.begin();
    static_assert(std::bidirectional_iterator<decltype(it)>);

    // Test with std::advance
    std::advance(it, 3);
    expect(it != bs_coll.begin());
    expect(it != bs_coll.end());

    auto advanced_bitset = *it;
    expect(!advanced_bitset.contains(63));
    advanced_bitset[63] = true;
    expect(advanced_bitset.contains(63));
  };

  "BitSets Iterator Bidirectional"_test = [] -> void {
    alloc::OwningArena<> alloc;
    BitSets bs_coll = BitSets::empty(&alloc, 6, 128);

    // Test forward iteration
    auto it = bs_coll.begin();
    std::ptrdiff_t forward_count = 0;
    while (it != bs_coll.end()) {
      auto bitset = *it;
      bitset[forward_count] = true;
      ++forward_count;
      ++it;
    }
    expect(eq(forward_count, 6));

    // Test backward iteration
    it = bs_coll.end();
    std::ptrdiff_t backward_count = 0;
    while (it-- != bs_coll.begin()) {
      auto bitset = *it;
      expect(bitset.contains(5 - backward_count));
      ++backward_count;
    }
    expect(eq(backward_count, 6));

    // Test postfix decrement
    it = bs_coll.end();
    --it; // Move to last valid position
    auto last_it = it--;
    expect(last_it != it);
    auto bitset_at_it = *it;
    auto bitset_at_last = *last_it;

    // They should be different BitSets
    bitset_at_it[100] = true;
    bitset_at_last[101] = true;
    expect(bitset_at_it.contains(100) && !bitset_at_it.contains(101));
    expect(bitset_at_last.contains(101) && !bitset_at_last.contains(100));

    // Test mixed forward/backward
    it = bs_coll.begin();
    ++it;
    ++it; // Move forward 2
    --it; // Move back 1
    ++it;
    ++it; // Move forward 2 more

    std::ptrdiff_t expected_distance = 3;
    std::ptrdiff_t actual_distance = std::distance(bs_coll.begin(), it);
    expect(eq(actual_distance, expected_distance));
  };

  "BitSets Iterator Sizes"_test = [] -> void {
    alloc::OwningArena<> alloc;

    // Test empty collection
    BitSets empty_coll = BitSets::empty(&alloc, 0, 10);
    expect(empty_coll.begin() == empty_coll.end());
    expect(eq(std::distance(empty_coll.begin(), empty_coll.end()), 0));

    std::ptrdiff_t empty_count = 0;
    for ([[maybe_unused]] auto bs : empty_coll) ++empty_count;
    expect(eq(empty_count, 0));

    // Test single BitSet collection
    BitSets single_coll = BitSets::empty(&alloc, 1, 32);
    expect(single_coll.begin() != single_coll.end());
    expect(eq(std::distance(single_coll.begin(), single_coll.end()), 1));

    std::ptrdiff_t single_count = 0;
    for (auto bs : single_coll) {
      bs[15] = true;
      expect(bs.contains(15));
      ++single_count;
    }
    expect(eq(single_count, 1));

    // Test large collection
    BitSets large_coll = BitSets::empty(&alloc, 15, 200);
    expect(eq(std::distance(large_coll.begin(), large_coll.end()), 15));

    std::ptrdiff_t large_count = 0;
    for (auto bs : large_coll) {
      bs[large_count] = true;
      expect(bs.contains(large_count));
      ++large_count;
    }
    expect(eq(large_count, 15));

    // Test different num_elts values
    BitSets small_elts = BitSets::empty(&alloc, 3, 1);   // 1 element per BitSet
    BitSets med_elts = BitSets::empty(&alloc, 3, 64);    // Exactly 1 chunk
    BitSets large_elts = BitSets::empty(&alloc, 3, 129); // Multiple chunks

    for (auto &coll : {small_elts, med_elts, large_elts}) {
      expect(eq(std::distance(coll.begin(), coll.end()), 3));
      std::ptrdiff_t iter_count = 0;
      for (auto bs : coll) {
        bs[0] = true;
        expect(bs.contains(0));
        ++iter_count;
      }
      expect(eq(iter_count, 3));
    }
  };

  "BitSets Iterator Independence"_test = [] -> void {
    alloc::OwningArena<> alloc;
    BitSets bs_coll = BitSets::empty(&alloc, 4, 80);

    // Create multiple iterators
    auto it1 = bs_coll.begin();
    auto it2 = bs_coll.begin();
    auto it3 = bs_coll.begin();

    // Advance them to different positions
    ++it2;
    ++it3;
    ++it3;

    // Get BitSets through different iterators
    auto bs1 = *it1;
    auto bs2 = *it2;
    auto bs3 = *it3;

    // Modify each BitSet differently
    bs1[10] = true;
    bs2[20] = true;
    bs3[30] = true;

    // Verify independence
    expect(bs1.contains(10) && !bs1.contains(20) && !bs1.contains(30));
    expect(!bs2.contains(10) && bs2.contains(20) && !bs2.contains(30));
    expect(!bs3.contains(10) && !bs3.contains(20) && bs3.contains(30));

    // Test that advancing one iterator doesn't affect others
    ++it1;
    auto bs1_new = *it1;
    bs1_new[15] = true;

    // Original BitSet references should still be independent
    expect(bs1.contains(10) && !bs1.contains(15));
    expect(bs1_new.contains(15) && !bs1_new.contains(10));

    // Verify iterator positions remain correct
    expect(it2 == bs_coll.begin() + 1);
    expect(it3 == bs_coll.begin() + 2);
    expect(it1 == bs_coll.begin() + 1);

    // Test iterator stability during BitSet operations
    auto stable_it = bs_coll.begin();
    auto stable_bs = *stable_it;

    // Perform operations on the BitSet
    stable_bs[40] = true;
    stable_bs[50] = true;
    stable_bs.remove(40);

    // Iterator should still be valid and point to same position
    auto stable_bs_again = *stable_it;
    expect(stable_bs_again.contains(10));
    expect(!stable_bs_again.contains(40));
    expect(stable_bs_again.contains(50));
    expect(eq(stable_bs_again.size(), 2));
  };

  "BitSetIterator Plus Operator"_test = [] -> void {
    // Test operator+(std::ptrdiff_t) on BitSetIterator
    BitSet bs(200);

    // Set some scattered bits: 5, 12, 20, 35, 64, 100, 150
    bs[5] = true;
    bs[12] = true;
    bs[20] = true;
    bs[35] = true;
    bs[64] = true;
    bs[100] = true;
    bs[150] = true;

    auto it = bs.begin();
    expect(eq(*it, 5)); // First element

    // Test advancing by 0 (should return same position)
    auto same_it = it + 0;
    expect(eq(*same_it, 5));

    // Test advancing by 1
    auto next_it = it + 1;
    expect(eq(*next_it, 12));

    // Test advancing by 3
    auto third_it = it + 3;
    expect(eq(*third_it, 35));

    // Test advancing by 3
    auto fourth_it = it + 4;
    expect(eq(*fourth_it, 64));

    // Test advancing by all elements
    auto all_it = it + 7;
    expect(all_it == std::default_sentinel);

    // Test advancing beyond end
    auto beyond_it = it + 10;
    expect(beyond_it == std::default_sentinel);

    // Test from middle position
    auto mid_it = bs.begin() + 2; // Should be at position 20
    expect(eq(*mid_it, 20));

    auto from_mid = mid_it + 3; // Should advance to 100
    expect(eq(*from_mid, 100));

    // Test with negative value (should return same iterator)
    auto neg_it = it + (-5);
    expect(eq(*neg_it, 5)); // Should be same as original
  };

  "BitSetIterator Minus Operator"_test = [] -> void {
    // Test operator-(const BitSetIterator &other)
    BitSet bs(100);

    // Set bits at positions: 10, 25, 40, 55, 70, 85
    bs[10] = true;
    bs[25] = true;
    bs[40] = true;
    bs[55] = true;
    bs[70] = true;
    bs[85] = true;

    auto it_begin = bs.begin();
    expect(eq(*it_begin, 10));

    auto it_2nd = it_begin + 1;
    expect(eq(*it_2nd, 25));

    auto it_4th = it_begin + 3;
    expect(eq(*it_4th, 55));

    auto it_end = bs.begin() + 6;
    expect(it_end == std::default_sentinel);

    // Test difference calculation
    expect(eq(it_2nd - it_begin, 1));
    expect(eq(it_4th - it_begin, 3));
    expect(eq(it_4th - it_2nd, 2));
    expect(eq(it_begin - it_2nd, -1));
    expect(eq(it_2nd - it_4th, -2));

    // Test same iterator
    // NOLINTNEXTLINE(misc-redundant-expression)
    expect(eq(it_begin - it_begin, 0));
    // NOLINTNEXTLINE(misc-redundant-expression)
    expect(eq(it_2nd - it_2nd, 0));
  };

  "BitSetIterator Large Jump Performance"_test = [] -> void {
    // Test that operator+ efficiently handles large jumps using popcount
    BitSet bs(10000);

    // Set every 100th bit: 0, 100, 200, ..., 9900 (100 bits total)
    for (std::ptrdiff_t i = 0; i < 10000; i += 100) bs[i] = true;
    expect(eq(bs.size(), 100));

    auto it = bs.begin();
    expect(eq(*it, 0));

    // Test large jumps that should use popcount optimization
    auto jump_10 = it + 10;
    expect(eq(*jump_10,
              1000)); // 10th set bit (0-indexed) should be at position 1000

    auto jump_50 = it + 50;
    expect(eq(*jump_50,
              5000)); // 50th set bit (0-indexed) should be at position 5000

    auto jump_99 = it + 99;
    expect(eq(*jump_99,
              9900)); // 99th set bit (0-indexed) should be at position 9900

    auto jump_100 = it + 100;
    expect(jump_100 == std::default_sentinel); // Beyond end

    // Test from middle position
    auto mid_it = it + 25; // Should be at position 2500
    expect(eq(*mid_it, 2500));

    auto from_mid_jump = mid_it + 25; // Should jump to position 5000
    expect(eq(*from_mid_jump, 5000));
  };

  "BitSetIterator Dense BitSet"_test = [] -> void {
    // Test with a dense bitset to verify word-boundary handling
    BitSet bs = BitSet<>::dense(200);
    expect(eq(bs.size(), 200));

    auto it = bs.begin();
    expect(eq(*it, 0));

    // Test jumps that cross word boundaries (64-bit words)
    auto jump_63 = it + 63;
    expect(eq(*jump_63, 63));

    auto jump_64 = it + 64;
    expect(eq(*jump_64, 64));

    auto jump_65 = it + 65;
    expect(eq(*jump_65, 65));

    auto jump_128 = it + 128;
    expect(eq(*jump_128, 128));

    auto jump_199 = it + 199;
    expect(eq(*jump_199, 199));

    auto jump_200 = it + 200;
    expect(jump_200 == std::default_sentinel);

    // Test difference calculations with dense bitset
    auto it_0 = bs.begin();
    auto it_64 = it_0 + 64;
    auto it_128 = it_0 + 128;

    expect(eq(it_64 - it_0, 64));
    expect(eq(it_128 - it_64, 64));
    expect(eq(it_128 - it_0, 128));
  };

  "BitSetIterator Sparse BitSet"_test = [] -> void {
    // Test with very sparse bitset
    BitSet bs(1000);

    // Set only a few bits: 7, 8, 500, 999
    bs[7] = true;
    bs[8] = true;
    bs[250] = true;
    bs[500] = true;
    bs[749] = true;
    bs[999] = true;

    auto it = bs.begin();
    expect(*it == 7);

    expect(eq(*(it + 1), 8));
    expect(eq(*(it + 2), 250));
    expect(eq(*(it + 3), 500));
    expect(eq(*(it + 4), 749));
    expect(eq(*(it + 5), 999));
    expect((it + 6) == std::default_sentinel);
    expect((it + 7) == std::default_sentinel);
    expect((it + 777) == std::default_sentinel);

    // Test differences with sparse data
    for (std::ptrdiff_t i = 0; i < 5; ++i)
      for (std::ptrdiff_t j = 0; j < 5; ++j)
        expect(eq((it + i) - (it + j), i - j));
  };

  "BitSet ChunkedIterator"_test = [] -> void {
    // Test chunked iterator with example pattern: 0b100110111010
    // Expected chunks: _(1,2), _(4,6), _(7,10), _(11,12)
    BitSet bs(20);
    bs[1] = true;
    bs[4] = true;
    bs[5] = true;
    bs[7] = true;
    bs[8] = true;
    bs[9] = true;
    bs[11] = true;

    std::vector<std::pair<std::ptrdiff_t, std::ptrdiff_t>> expected_chunks = {
      {1, 2}, {4, 6}, {7, 10}, {11, 12}};

    std::ptrdiff_t chunk_idx = 0;
    for (auto chunk : bs.chunks()) {
      expect(chunk_idx < std::ssize(expected_chunks));
      expect(eq(chunk.begin(), expected_chunks[chunk_idx].first));
      expect(eq(chunk.end(), expected_chunks[chunk_idx].second));
      ++chunk_idx;
    }
    expect(eq(chunk_idx, expected_chunks.size()));
  };

  "BitSet ChunkedIterator Dense"_test = [] -> void {
    // Test with dense bitset - should produce one large chunk
    BitSet bs = BitSet<>::dense(100);

    std::ptrdiff_t chunk_count = 0;
    for (auto chunk : bs.chunks()) {
      expect(eq(chunk.begin(), chunk_count ? 64 : 0));
      expect(eq(chunk.end(), chunk_count ? 100 : 64));
      ++chunk_count;
    }
    expect(eq(chunk_count, 2));
  };

  "BitSet ChunkedIterator Empty"_test = [] -> void {
    // Test with empty bitset
    BitSet bs(100);

    std::ptrdiff_t chunk_count = 0;
    for (auto chunk : bs.chunks()) {
      ++chunk_count;
      utils::println("Unexpected chunk: ", chunk.begin(), "-", chunk.end());
    }
    expect(eq(chunk_count, 0));
    bs[4] = true;
    bs[5] = true;
    std::ptrdiff_t count = 0;
    for (auto chunk : bs.chunks()) {
      expect(eq(chunk.begin(), 4));
      expect(eq(chunk.end(), 6));
      ++count;
    }
    expect(eq(count, 1));
  };

  "BitSet ChunkedIterator SingleBits"_test = [] -> void {
    // Test with isolated single bits
    BitSet bs(100);
    bs[10] = true;
    bs[30] = true;
    bs[50] = true;
    bs[90] = true;

    std::vector<std::pair<std::ptrdiff_t, std::ptrdiff_t>> expected_chunks = {
      {10, 11}, {30, 31}, {50, 51}, {90, 91}};

    std::ptrdiff_t chunk_idx = 0;
    for (auto chunk : bs.chunks()) {
      expect(chunk_idx < std::ssize(expected_chunks));
      expect(eq(chunk.begin(), expected_chunks[chunk_idx].first));
      expect(eq(chunk.end(), expected_chunks[chunk_idx].second));
      ++chunk_idx;
    }
    expect(eq(chunk_idx, expected_chunks.size()));
  };

  "BitSet ChunkedIterator CrossingWords"_test = [] -> void {
    // Test chunks that cross 64-bit word boundaries
    BitSet bs(200);

    // Set bits around 64-bit boundaries
    for (std::ptrdiff_t i = 62; i <= 66; ++i) bs[i] = true; // 62,63,64,65,66
    for (std::ptrdiff_t i = 126; i <= 130; ++i)
      bs[i] = true; // 126,127,128,129,130

    std::vector<std::pair<std::ptrdiff_t, std::ptrdiff_t>> expected_chunks = {
      {62, 64}, {64, 67}, {126, 128}, {128, 131}};

    std::ptrdiff_t chunk_idx = 0;
    for (auto chunk : bs.chunks()) {
      expect(chunk_idx < std::ssize(expected_chunks));
      expect(eq(chunk.begin(), expected_chunks[chunk_idx].first));
      expect(eq(chunk.end(), expected_chunks[chunk_idx].second));
      ++chunk_idx;
    }
    expect(eq(chunk_idx, expected_chunks.size()));
  };

  return 0;
}

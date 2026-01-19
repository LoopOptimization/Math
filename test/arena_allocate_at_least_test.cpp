import boost.ut;
import Arena;
import std;

using namespace boost::ut;
using namespace alloc;

auto main() -> int {
  "ArenaAllocateAtLeast_Basic"_test = [] {
    OwningArena<> arena;

    auto [ptr, count] = arena.allocate_at_least<int>(1);
    expect(ptr != nullptr);
    expect(count >= 1);

    // Verify architecture-specific alignment
#ifdef __AVX512F__
    expect(count == 16); // 64 bytes / 4 bytes per int
#elifdef __AVX__
    expect(count == 8); // 32 bytes / 4 bytes per int
#else
    expect(count == 4); // 16 bytes / 4 bytes per int
#endif

    // Verify we can write to all allocated elements
    for (std::ptrdiff_t i = 0; i < count; ++i) ptr[i] = static_cast<int>(i);
    for (std::ptrdiff_t i = 0; i < count; ++i)
      expect(ptr[i] == static_cast<int>(i));
  };

  "ArenaAllocateAtLeast_Multiple"_test = [] {
    OwningArena<> arena;

    // Request more than one alignment unit
    auto [ptr, count] = arena.allocate_at_least<std::int64_t>(100);
    expect(ptr != nullptr);
    expect(count >= 100);

    // Verify proper alignment
    expect(reinterpret_cast<std::uintptr_t>(ptr) % alignof(std::int64_t) == 0);

    // Use the allocated memory
    for (std::ptrdiff_t i = 0; i < count; ++i) { ptr[i] = i * i; }
    // Verify written values
    for (std::ptrdiff_t i = 0; i < 100; ++i) { expect(ptr[i] == i * i); }
  };

  "ArenaAllocateAtLeast_ZeroRequest"_test = [] {
    OwningArena<> arena;

    auto [ptr, count] = arena.allocate_at_least<char>(0);
    expect(ptr != nullptr);

    // When requesting 0 elements, count can be 0 (aligned_up(0) = 0)
    // This is correct behavior - we get at least 0 elements
    expect(count >= 0);
  };

  "ArenaAllocateAtLeast_LargeType"_test = [] {
    OwningArena<> arena;

    struct Large {
      char data[128];
    };

    auto [ptr, count] = arena.allocate_at_least<Large>(1);
    expect(ptr != nullptr);
    expect(count >= 1);

    // Verify usability
    ptr[0].data[0] = 'A';
    expect(ptr[0].data[0] == 'A');
  };

  "ArenaAllocateAtLeast_CustomAlignment"_test = [] {
    OwningArena<> arena;

    auto [ptr, count] = arena.allocate_at_least<int, 128>(10);
    expect(ptr != nullptr);
    expect(count >= 10);

    // Verify 128-byte alignment
    expect(reinterpret_cast<std::uintptr_t>(ptr) % 128 == 0);

    // Use the memory
    for (std::ptrdiff_t i = 0; i < count; ++i) ptr[i] = static_cast<int>(i * 2);
  };

  "ArenaAllocateAtLeast_DifferentTypes"_test = [] {
    OwningArena<> arena;

    auto [c_ptr, c_count] = arena.allocate_at_least<char>(1);
    auto [s_ptr, s_count] = arena.allocate_at_least<short>(1);
    auto [i_ptr, i_count] = arena.allocate_at_least<int>(1);
    auto [l_ptr, l_count] = arena.allocate_at_least<long>(1);
    auto [d_ptr, d_count] = arena.allocate_at_least<double>(1);

    expect(c_ptr != nullptr && c_count >= 1);
    expect(s_ptr != nullptr && s_count >= 1);
    expect(i_ptr != nullptr && i_count >= 1);
    expect(l_ptr != nullptr && l_count >= 1);
    expect(d_ptr != nullptr && d_count >= 1);

    // Smaller types should get more elements for same alignment
    expect(c_count >= i_count);

    // Use the memory to verify it works
    c_ptr[0] = 'X';
    s_ptr[0] = 42;
    i_ptr[0] = 12345;
    l_ptr[0] = 999999;
    d_ptr[0] = 3.14159;

    expect(c_ptr[0] == 'X');
    expect(s_ptr[0] == 42);
    expect(i_ptr[0] == 12345);
    expect(l_ptr[0] == 999999);
    expect(d_ptr[0] == 3.14159);
  };

  "ArenaAllocateAtLeast_WithReset"_test = [] {
    OwningArena<> arena;

    auto [ptr1, count1] = arena.allocate_at_least<int>(5);
    for (std::ptrdiff_t i = 0; i < count1; ++i) ptr1[i] = 42;

    arena.reset();

    auto [ptr2, count2] = arena.allocate_at_least<int>(5);
    expect(ptr2 != nullptr);
    expect(count2 >= 5);

    // After reset, new allocation should work fine
    for (std::ptrdiff_t i = 0; i < count2; ++i) ptr2[i] = 99;
    expect(ptr2[0] == 99);
  };

  "ArenaAllocateAtLeast_ExactAlignmentMatch"_test = [] {
    OwningArena<> arena;

#ifdef __AVX512F__
    // Request exactly 16 ints (64 bytes with AVX512)
    auto [ptr, count] = arena.allocate_at_least<int>(16);
    expect(count == 16);
#elifdef __AVX__
    // Request exactly 8 ints (32 bytes with AVX)
    auto [ptr, count] = arena.allocate_at_least<int>(8);
    expect(count == 8);
#else
    // Request exactly 4 ints (16 bytes baseline)
    auto [ptr, count] = arena.allocate_at_least<int>(4);
    expect(count == 4);
#endif

    expect(ptr != nullptr);

    // Use all allocated memory
    for (std::ptrdiff_t i = 0; i < count; ++i) ptr[i] = static_cast<int>(i * 3);
  };

  "ArenaReallocateAtLeast_BasicGrowth"_test = [] {
    OwningArena<> arena;

    // Initial allocation
    auto [ptr1, count1] = arena.allocate_at_least<int>(5);
    expect(ptr1 != nullptr);
    expect(count1 >= 5);

    // Fill initial data
    for (std::ptrdiff_t i = 0; i < 5; ++i) {
      ptr1[i] = static_cast<int>(i * 10);
    }

    // Reallocate to larger size
    auto [ptr2, count2] = arena.reallocate_at_least(ptr1, count1, 20);
    expect(ptr2 != nullptr);
    expect(count2 >= 20);

    // Verify data was preserved
    for (std::ptrdiff_t i = 0; i < 5; ++i) {
      expect(ptr2[i] == static_cast<int>(i * 10));
    }

    // Use the new space
    for (std::ptrdiff_t i = 5; i < count2; ++i) {
      ptr2[i] = static_cast<int>(i * 2);
    }
  };

  "ArenaReallocateAtLeast_Shrinking"_test = [] {
    OwningArena<> arena;

    // Initial allocation
    auto [ptr1, count1] = arena.allocate_at_least<long>(50);
    expect(ptr1 != nullptr);
    expect(count1 >= 50);

    // Fill data
    for (std::ptrdiff_t i = 0; i < 50; ++i) ptr1[i] = i * 100L;

    // Reallocate to smaller size - should return same pointer
    auto [ptr2, count2] = arena.reallocate_at_least(ptr1, count1, 25);
    expect(ptr2 == ptr1); // Same pointer when shrinking
    expect(count2 >= 25);

    // Verify data is still intact
    for (std::ptrdiff_t i = 0; i < 25; ++i) expect(ptr2[i] == i * 100L);
  };

  "ArenaReallocateAtLeast_WithOldSize"_test = [] {
    OwningArena<> arena;

    // Allocate space for 20 elements but only use 10
    auto [ptr1, count1] = arena.allocate_at_least<double>(20);
    expect(count1 >= 20);

    // Fill only first 10 elements
    for (std::ptrdiff_t i = 0; i < 10; ++i)
      ptr1[i] = static_cast<double>(i) * 3.14;

    // Reallocate specifying oldSize = 10 (only copy 10 elements)
    auto [ptr2, count2] = arena.reallocate_at_least(ptr1, count1, 40, 10);
    expect(ptr2 != nullptr);
    expect(count2 >= 40);

    // Verify only the first 10 elements were copied
    for (std::ptrdiff_t i = 0; i < 10; ++i)
      expect(ptr2[i] == static_cast<double>(i) * 3.14);
  };

  "ArenaReallocateAtLeast_CustomAlignment"_test = [] {
    OwningArena<> arena;

    // Initial allocation with custom alignment
    auto [ptr1, count1] = arena.allocate_at_least<int, 128>(10);
    expect(ptr1 != nullptr);
    expect(count1 >= 10);
    expect(reinterpret_cast<std::uintptr_t>(ptr1) % 128 == 0);

    // Fill data
    for (std::ptrdiff_t i = 0; i < 10; ++i) ptr1[i] = static_cast<int>(i * 7);

    // Reallocate with same custom alignment
    auto [ptr2, count2] =
      arena.reallocate_at_least<false, int, 128>(ptr1, count1, 30);
    expect(ptr2 != nullptr);
    expect(count2 >= 30);
    expect(reinterpret_cast<std::uintptr_t>(ptr2) % 128 == 0);

    // Verify data preserved
    for (std::ptrdiff_t i = 0; i < 10; ++i)
      expect(ptr2[i] == static_cast<int>(i * 7));
  };

  "ArenaReallocateAtLeast_ForOverwrite"_test = [] {
    OwningArena<> arena;

    // Initial allocation
    auto [ptr1, count1] = arena.allocate_at_least<int>(8);
    expect(count1 >= 8);

    // Fill with data
    for (std::ptrdiff_t i = 0; i < 8; ++i) ptr1[i] = 42;

    // Reallocate with ForOverwrite=true (data may not be preserved)
    auto [ptr2, count2] = arena.reallocate_at_least<true>(ptr1, count1, 20);
    expect(ptr2 != nullptr);
    expect(count2 >= 20);

    // Just verify the allocation worked, don't check old data
    for (std::ptrdiff_t i = 0; i < count2; ++i) ptr2[i] = static_cast<int>(i);
  };

  "ArenaReallocateAtLeast_MultipleReallocations"_test = [] {
    OwningArena<> arena;

    // Start with small allocation
    auto [ptr, count] = arena.allocate_at_least<std::int64_t>(2);
    expect(count >= 2);
    ptr[0] = 100;
    ptr[1] = 200;
#ifdef __AVX512F__
    expect(count == 8);
#elifdef __AVX__
    expect(count == 4);
#else
    expect(count == 2);
#endif

    // Grow multiple times
    auto [ptr2, count2] = arena.reallocate_at_least(ptr, count, 10);
    expect(count2 >= 10);
    expect(ptr2[0] == 100);
    expect(ptr2[1] == 200);
#ifdef __AVX512F__
    expect(count2 == 16);
#elifdef __AVX__
    expect(count2 == 12);
#else
    expect(count2 == 10);
#endif

    auto [ptr3, count3] = arena.reallocate_at_least(ptr2, count2, 50);
    expect(count3 >= 50);
    expect(ptr3[0] == 100);
    expect(ptr3[1] == 200);
#ifdef __AVX512F__
    expect(count3 == 56);
#elifdef __AVX__
    expect(count3 == 52);
#else
    expect(count3 == 50);
#endif

    auto [ptr4, count4] = arena.reallocate_at_least(ptr3, count3, 100);
    expect(count4 >= 100);
    expect(ptr4[0] == 100);
    expect(ptr4[1] == 200);
#ifdef __AVX512F__
    expect(count4 == 104);
#elifdef __AVX__
    expect(count4 == 100);
#else
    expect(count4 == 100);
#endif
  };

  "ArenaReallocateAtLeast_ArchitectureSpecific"_test = [] {
    OwningArena<> arena;

    // Small initial allocation
    auto [ptr1, count1] = arena.allocate_at_least<char>(1);
    ptr1[0] = 'X';

    // Reallocate - should get at least one alignment unit
    auto [ptr2, count2] = arena.reallocate_at_least(ptr1, count1, 1);

    expect(count2 == count1);
#ifdef __AVX512F__
    expect(count2 == 64);
#elifdef __AVX__
    expect(count2 == 32);
#else
    expect(count2 == 16);
#endif

    expect(ptr2[0] == 'X'); // Data preserved
  };

  return 0;
}

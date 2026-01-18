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
    for (std::ptrdiff_t i = 0; i < count; ++i) {
      ptr[i] = static_cast<int>(i);
    }
    for (std::ptrdiff_t i = 0; i < count; ++i) {
      expect(ptr[i] == static_cast<int>(i));
    }
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
    for (std::ptrdiff_t i = 0; i < count; ++i) {
      ptr[i] = i * i;
    }
    // Verify written values
    for (std::ptrdiff_t i = 0; i < 100; ++i) {
      expect(ptr[i] == i * i);
    }
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
    for (std::ptrdiff_t i = 0; i < count; ++i) {
      ptr[i] = static_cast<int>(i * 2);
    }
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
    for (std::ptrdiff_t i = 0; i < count1; ++i) {
      ptr1[i] = 42;
    }

    arena.reset();

    auto [ptr2, count2] = arena.allocate_at_least<int>(5);
    expect(ptr2 != nullptr);
    expect(count2 >= 5);

    // After reset, new allocation should work fine
    for (std::ptrdiff_t i = 0; i < count2; ++i) {
      ptr2[i] = 99;
    }
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
    for (std::ptrdiff_t i = 0; i < count; ++i) {
      ptr[i] = static_cast<int>(i * 3);
    }
  };

  return 0;
}

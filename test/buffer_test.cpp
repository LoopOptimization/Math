import boost.ut;
import Buffer;
import std;

using namespace boost::ut;
using containers::Buffer;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
auto main() -> int {

  "Buffer DefaultConstruction"_test = [] {
    Buffer<int> buf;
    expect(buf.empty());
    expect(eq(buf.size(), 0u));
    expect(eq(buf.capacity(), 0u));
    expect(buf.data() == nullptr);
  };

  "Buffer ConstructWithSize"_test = [] {
    Buffer<int> buf(10);
    expect(!buf.empty());
    expect(eq(buf.size(), 10u));
    expect(buf.capacity() >= 10u);
    expect(buf.data() != nullptr);

    // Elements should be default initialized
    for (std::uint32_t i = 0; i < buf.size(); ++i) {
      buf[i] = static_cast<int>(i);
      expect(eq(buf[i], static_cast<int>(i)));
    }
  };

  "Buffer ConstructWithSizeAndValue"_test = [] {
    Buffer<int> buf(5, 42);
    expect(eq(buf.size(), 5u));
    expect(buf.capacity() >= 5u);

    for (std::uint32_t i = 0; i < buf.size(); ++i) expect(eq(buf[i], 42));
  };

  "Buffer InitializerList"_test = [] {
    Buffer<int> buf{1, 2, 3, 4, 5};
    expect(eq(buf.size(), 5u));
    expect(buf.capacity() >= 5u);

    for (std::uint32_t i = 0; i < buf.size(); ++i)
      expect(eq(buf[i], static_cast<int>(i + 1)));
  };

  "Buffer CopyConstruction"_test = [] {
    Buffer<int> buf1{10, 20, 30, 40, 50};
    Buffer<int> buf2(buf1);

    expect(eq(buf2.size(), buf1.size()));
    expect(buf2.data() != buf1.data()); // Different storage

    for (std::uint32_t i = 0; i < buf1.size(); ++i)
      expect(eq(buf2[i], buf1[i]));
  };

  "Buffer MoveConstruction"_test = [] {
    Buffer<int> buf1{1, 2, 3, 4, 5};
    auto *old_data = buf1.data();
    auto old_size = buf1.size();
    auto old_capacity = buf1.capacity();

    Buffer<int> buf2(std::move(buf1));

    expect(eq(buf2.size(), old_size));
    expect(eq(buf2.capacity(), old_capacity));
    expect(buf2.data() == old_data); // Moved data pointer

    expect(buf1.empty());
    expect(eq(buf1.size(), 0u));
    expect(eq(buf1.capacity(), 0u));
    expect(buf1.data() == nullptr);

    expect(eq(buf2[0], 1));
    expect(eq(buf2[4], 5));
  };

  "Buffer CopyAssignment"_test = [] {
    Buffer<int> buf1{1, 2, 3};
    Buffer<int> buf2{10, 20, 30, 40, 50};

    buf2 = buf1;
    expect(eq(buf2.size(), buf1.size()));

    for (std::uint32_t i = 0; i < buf1.size(); ++i)
      expect(eq(buf2[i], buf1[i]));
  };

  "Buffer CopyAssignmentReuseStorage"_test = [] {
    Buffer<int> buf1{1, 2, 3, 4, 5};
    Buffer<int> buf2;
    buf2.reserve(10);
    auto old_capacity = buf2.capacity();

    buf2 = buf1;
    expect(eq(buf2.size(), buf1.size()));
    expect(eq(buf2.capacity(), old_capacity)); // Reused storage

    for (std::uint32_t i = 0; i < buf1.size(); ++i)
      expect(eq(buf2[i], buf1[i]));
  };

  "Buffer MoveAssignment"_test = [] {
    Buffer<int> buf1{1, 2, 3, 4, 5};
    Buffer<int> buf2{10, 20};

    auto *old_data = buf1.data();
    auto old_size = buf1.size();

    buf2 = std::move(buf1);

    expect(eq(buf2.size(), old_size));
    expect(buf2.data() == old_data);

    expect(buf1.empty());
    expect(buf1.data() == nullptr);

    expect(eq(buf2[0], 1));
    expect(eq(buf2[4], 5));
  };

  "Buffer InitializerListAssignment"_test = [] {
    Buffer<int> buf{1, 2, 3};
    buf = {10, 20, 30, 40, 50};

    expect(eq(buf.size(), 5u));
    expect(eq(buf[0], 10));
    expect(eq(buf[4], 50));
  };

  "Buffer ElementAccess"_test = [] {
    Buffer<int> buf{10, 20, 30, 40, 50};

    expect(eq(buf[0], 10));
    expect(eq(buf[4], 50));
    expect(eq(buf.front(), 10));
    expect(eq(buf.back(), 50));

    buf[2] = 99;
    expect(eq(buf[2], 99));

    buf.front() = 11;
    buf.back() = 55;
    expect(eq(buf[0], 11));
    expect(eq(buf[4], 55));
  };

  "Buffer At"_test = [] {
    Buffer<int> buf{1, 2, 3, 4, 5};

    expect(eq(buf.at(0), 1));
    expect(eq(buf.at(4), 5));

    // bool caught = false;
    // try {
    //   [[maybe_unused]] auto val = buf.at(10);
    // } catch (const std::out_of_range &) { caught = true; }
    // expect(caught);
  };

  "Buffer Iterators"_test = [] {
    Buffer<int> buf{1, 2, 3, 4, 5};

    int sum = 0;
    for (auto it = buf.begin(); it != buf.end(); ++it) sum += *it;
    expect(eq(sum, 15));

    sum = 0;
    for (const auto &val : buf) sum += val;
    expect(eq(sum, 15));

    // Test reverse iterators
    sum = 0;
    for (auto it = buf.rbegin(); it != buf.rend(); ++it) sum += *it;
    expect(eq(sum, 15));
  };

  "Buffer PushBack"_test = [] {
    Buffer<int> buf;
    expect(buf.empty());

    buf.push_back(10);
    expect(eq(buf.size(), 1u));
    expect(eq(buf[0], 10));

    buf.push_back(20);
    expect(eq(buf.size(), 2u));
    expect(eq(buf[1], 20));

    // Test push_back with rvalue
    buf.push_back(30);
    expect(eq(buf.size(), 3u));
    expect(eq(buf[2], 30));
  };

  "Buffer EmplaceBack"_test = [] {
    Buffer<std::string> buf;

    auto &ref = buf.emplace_back("hello");
    expect(eq(buf.size(), 1u));
    expect(buf[0] == "hello");
    expect(ref == "hello");

    buf.emplace_back(5, 'x');
    expect(eq(buf.size(), 2u));
    expect(buf[1] == "xxxxx");
  };

  "Buffer PopBack"_test = [] {
    Buffer<int> buf{1, 2, 3, 4, 5};

    buf.pop_back();
    expect(eq(buf.size(), 4u));
    expect(eq(buf.back(), 4));

    auto val = buf.pop_back_val();
    expect(eq(val, 4));
    expect(eq(buf.size(), 3u));
    expect(eq(buf.back(), 3));
  };

  "Buffer Reserve"_test = [] {
    Buffer<int> buf;
    expect(eq(buf.capacity(), 0u));

    buf.reserve(100);
    expect(buf.capacity() >= 100u);
    expect(eq(buf.size(), 0u));

    // Reserve smaller size should not reduce capacity
    auto old_capacity = buf.capacity();
    buf.reserve(50);
    expect(eq(buf.capacity(), old_capacity));
  };

  "Buffer Resize"_test = [] {
    Buffer<int> buf{1, 2, 3};

    buf.resize(5);
    expect(eq(buf.size(), 5u));
    expect(eq(buf[0], 1));
    expect(eq(buf[1], 2));
    expect(eq(buf[2], 3));

    buf.resize(2);
    expect(eq(buf.size(), 2u));
    expect(eq(buf[0], 1));
    expect(eq(buf[1], 2));
  };

  "Buffer ResizeWithValue"_test = [] {
    Buffer<int> buf{1, 2, 3};

    buf.resize(6, 42);
    expect(eq(buf.size(), 6u));
    expect(eq(buf[0], 1));
    expect(eq(buf[2], 3));
    expect(eq(buf[3], 42));
    expect(eq(buf[5], 42));
  };

  "Buffer Clear"_test = [] {
    Buffer<int> buf{1, 2, 3, 4, 5};
    auto old_capacity = buf.capacity();

    buf.clear();
    expect(buf.empty());
    expect(eq(buf.size(), 0u));
    expect(eq(buf.capacity(), old_capacity)); // Capacity unchanged
  };

  "Buffer ShrinkToFit"_test = [] {
    Buffer<int> buf;
    buf.reserve(100);
    expect(buf.capacity() >= 100u);

    buf.push_back(1);
    buf.push_back(2);
    buf.push_back(3);

    buf.shrink_to_fit();
    expect(buf.capacity() >= buf.size());
    expect(buf.capacity() < 100u);
    expect(eq(buf.size(), 3u));
    expect(eq(buf[0], 1));
    expect(eq(buf[2], 3));
  };

  "Buffer Insert"_test = [] {
    Buffer<int> buf{1, 2, 4, 5};

    auto it = buf.insert(buf.begin() + 2, 3);
    expect(eq(buf.size(), 5u));
    expect(eq(*it, 3));
    expect(eq(buf[0], 1));
    expect(eq(buf[1], 2));
    expect(eq(buf[2], 3));
    expect(eq(buf[3], 4));
    expect(eq(buf[4], 5));

    // Insert at beginning
    it = buf.insert(buf.begin(), 0);
    expect(eq(buf.size(), 6u));
    expect(eq(*it, 0));
    expect(eq(buf[0], 0));
    expect(eq(buf[1], 1));

    // Insert at end
    it = buf.insert(buf.end(), 6);
    expect(eq(buf.size(), 7u));
    expect(eq(*it, 6));
    expect(eq(buf[6], 6));
  };

  "Buffer InsertMove"_test = [] {
    Buffer<std::string> buf;
    buf.push_back("hello");
    buf.push_back("world");

    std::string temp = "inserted";
    auto it = buf.insert(buf.begin() + 1, std::move(temp));
    expect(eq(buf.size(), 3u));
    expect(*it == "inserted");
    expect(buf[0] == "hello");
    expect(buf[1] == "inserted");
    expect(buf[2] == "world");
  };

  "Buffer Erase"_test = [] {
    Buffer<int> buf{1, 2, 3, 4, 5};

    auto it = buf.erase(buf.begin() + 2);
    expect(eq(buf.size(), 4u));
    expect(eq(*it, 4));
    expect(eq(buf[0], 1));
    expect(eq(buf[1], 2));
    expect(eq(buf[2], 4));
    expect(eq(buf[3], 5));

    // Erase at beginning
    it = buf.erase(buf.begin());
    expect(eq(buf.size(), 3u));
    expect(eq(*it, 2));
    expect(eq(buf[0], 2));

    // Erase at end
    it = buf.erase(buf.end() - 1);
    expect(eq(buf.size(), 2u));
    expect(it == buf.end());
  };

  "Buffer EraseRange"_test = [] {
    Buffer<int> buf{1, 2, 3, 4, 5, 6, 7, 8};

    auto it = buf.erase(buf.begin() + 2, buf.begin() + 5);
    expect(eq(buf.size(), 5u));
    expect(eq(*it, 6));
    expect(eq(buf[0], 1));
    expect(eq(buf[1], 2));
    expect(eq(buf[2], 6));
    expect(eq(buf[3], 7));
    expect(eq(buf[4], 8));

    // Erase all
    it = buf.erase(buf.begin(), buf.end());
    expect(buf.empty());
    expect(it == buf.end());
  };

  "Buffer Swap"_test = [] {
    Buffer<int> buf1{1, 2, 3};
    Buffer<int> buf2{10, 20, 30, 40, 50};

    auto *data1 = buf1.data();
    auto *data2 = buf2.data();
    auto size1 = buf1.size();
    auto size2 = buf2.size();

    buf1.swap(buf2);

    expect(buf1.data() == data2);
    expect(buf2.data() == data1);
    expect(eq(buf1.size(), size2));
    expect(eq(buf2.size(), size1));

    expect(eq(buf1[0], 10));
    expect(eq(buf2[0], 1));
  };

  "Buffer Comparison"_test = [] {
    Buffer<int> buf1{1, 2, 3};
    Buffer<int> buf2{1, 2, 3};
    Buffer<int> buf3{1, 2, 4};
    Buffer<int> buf4{1, 2};

    expect(buf1 == buf2);
    expect(!(buf1 != buf2));

    expect(buf1 != buf3);
    expect(buf1 < buf3);
    expect(buf3 > buf1);
    expect(buf1 <= buf3);
    expect(buf3 >= buf1);

    expect(buf4 < buf1);
    expect(buf1 > buf4);
  };

  "Buffer Growth"_test = [] {
    Buffer<int> buf;

    // Test that capacity grows appropriately
    for (std::uint32_t i = 0; i < 100; ++i) buf.push_back(static_cast<int>(i));

    expect(eq(buf.size(), 100u));
    for (std::uint32_t i = 0; i < 100; ++i)
      expect(eq(buf[i], static_cast<int>(i)));
  };

  "Buffer NonTrivialType"_test = [] {
    Buffer<std::string> buf;

    buf.push_back("one");
    buf.push_back("two");
    buf.push_back("three");

    expect(eq(buf.size(), 3u));
    expect(buf[0] == "one");
    expect(buf[1] == "two");
    expect(buf[2] == "three");

    buf.pop_back();
    expect(eq(buf.size(), 2u));

    buf.clear();
    expect(buf.empty());
  };

  "Buffer MaxSize"_test = [] {
    Buffer<int> buf;
    expect(eq(buf.max_size(), std::numeric_limits<std::uint32_t>::max()));
  };

  "Buffer LargeAllocation"_test = [] {
    // Test that allocate_at_least is being used
    Buffer<int> buf;
    buf.reserve(10);

    // Capacity should be >= 10, potentially larger due to allocate_at_least
    auto cap = buf.capacity();
    expect(cap >= 10u);

    // Fill to capacity without reallocation
    for (std::uint32_t i = 0; i < cap; ++i) buf.push_back(static_cast<int>(i));

    expect(eq(buf.size(), cap));
  };

  "Buffer EmptyOperations"_test = [] {
    Buffer<int> buf;

    // Test iterators on empty buffer
    expect(buf.begin() == buf.end());
    expect(buf.cbegin() == buf.cend());
    expect(buf.rbegin() == buf.rend());

    // Test clear on empty buffer
    buf.clear();
    expect(buf.empty());

    // Test shrink_to_fit on empty buffer
    buf.shrink_to_fit();
    expect(buf.empty());
  };

  "Buffer ConstCorrectness"_test = [] {
    const Buffer<int> buf{1, 2, 3, 4, 5};

    expect(eq(buf.size(), 5u));
    expect(eq(buf[0], 1));
    expect(eq(buf.front(), 1));
    expect(eq(buf.back(), 5));
    expect(buf.data() != nullptr);

    int sum = 0;
    for (const auto &val : buf) sum += val;
    expect(eq(sum, 15));

    auto it = buf.cbegin();
    expect(eq(*it, 1));
  };

  "Buffer InsertWithReallocation"_test = [] {
    Buffer<int> buf;
    buf.reserve(2);
    buf.push_back(1);
    buf.push_back(3);

    // This should trigger reallocation
    auto it = buf.insert(buf.begin() + 1, 2);
    expect(eq(buf.size(), 3u));
    expect(eq(*it, 2));
    expect(eq(buf[0], 1));
    expect(eq(buf[1], 2));
    expect(eq(buf[2], 3));
  };

  "Buffer StressTest"_test = [] {
    Buffer<int> buf;

    // Stress test: many insertions and deletions
    for (std::uint32_t i = 0; i < 1000; ++i) buf.push_back(static_cast<int>(i));

    expect(eq(buf.size(), 1000u));

    for (std::uint32_t i = 0; i < 500; ++i) buf.pop_back();

    expect(eq(buf.size(), 500u));

    for (std::uint32_t i = 0; i < buf.size(); ++i)
      expect(eq(buf[i], static_cast<int>(i)));

    buf.clear();
    expect(buf.empty());
  };

  "Buffer ReserveGrow"_test = [] {
    Buffer<int> buf;
    expect(eq(buf.capacity(), 0u));

    // First reserveGrow should allocate at least the requested capacity
    buf.reserveGrow(10);
    auto first_cap = buf.capacity();
    expect(first_cap >= 10u);
    expect(eq(buf.size(), 0u));

    // reserveGrow should grow by 2x (capacity + capacity) when current
    // capacity is non-zero
    buf.reserveGrow(first_cap + 1);
    auto second_cap = buf.capacity();
    expect(second_cap > first_cap);
    // Should be at least 2x the first capacity due to growth strategy
    expect(second_cap >= first_cap * 2);

    // reserveGrow with smaller size should not reduce capacity
    auto old_capacity = buf.capacity();
    buf.reserveGrow(5);
    expect(eq(buf.capacity(), old_capacity));
  };

  "Buffer ResizeGrow"_test = [] {
    Buffer<int> buf{1, 2, 3};
    auto initial_cap = buf.capacity();

    // resizeGrow should increase size and use growth strategy for capacity
    buf.resizeGrow(initial_cap + 5);
    expect(eq(buf.size(), initial_cap + 5));
    // Should grow by 2x, so capacity should be at least 2x initial capacity
    expect(buf.capacity() >= initial_cap * 2);
    // Original elements should be preserved
    expect(eq(buf[0], 1));
    expect(eq(buf[1], 2));
    expect(eq(buf[2], 3));

    // Test with value type that requires initialization
    Buffer<int> buf2;
    buf2.resizeGrow(10);
    expect(eq(buf2.size(), 10u));
    expect(buf2.capacity() >= 10u);
  };

  "Buffer PushBackReturnReference"_test = [] {
    Buffer<int> buf;

    // Test that push_back returns a reference to the inserted element
    auto &ref1 = buf.push_back(10);
    expect(eq(ref1, 10));
    expect(eq(buf.size(), 1u));
    expect(eq(buf[0], 10));

    // Verify the reference points to the actual element in the buffer
    ref1 = 20;
    expect(eq(buf[0], 20));

    // Test with rvalue
    auto &ref2 = buf.push_back(30);
    expect(eq(ref2, 30));
    expect(eq(buf.size(), 2u));
    expect(eq(buf[1], 30));

    // Test with non-trivial type
    Buffer<std::string> str_buf;
    auto &str_ref = str_buf.push_back("hello");
    expect(str_ref == "hello");
    str_ref = "world";
    expect(str_buf[0] == "world");

    // Test push_back with move
    std::string temp = "moved";
    auto &moved_ref = str_buf.push_back(std::move(temp));
    expect(moved_ref == "moved");
    expect(str_buf[1] == "moved");
  };

  "Buffer ReserveGrowVsReserve"_test = [] {
    Buffer<int> buf1;
    Buffer<int> buf2;

    // Both start with reserve to get initial capacity
    buf1.reserve(10);
    buf2.reserve(10);

    auto cap1 = buf1.capacity();
    auto cap2 = buf2.capacity();
    expect(eq(cap1, cap2));

    // reserve with larger value
    buf1.reserve(cap1 + 1);
    auto new_cap1 = buf1.capacity();

    // reserveGrow with same larger value
    buf2.reserveGrow(cap1 + 1);
    auto new_cap2 = buf2.capacity();

    // reserveGrow should allocate more due to growth strategy
    expect(ge(new_cap2, new_cap1));
    // reserveGrow should be at least 2x the old capacity
    expect(ge(new_cap2, cap2 * 2));
  };

  return 0;
}

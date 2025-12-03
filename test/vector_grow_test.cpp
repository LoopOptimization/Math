import boost.ut;

using namespace boost::ut;
import ManagedArray;
import std;

// Test for reserveGrow and resizeGrow methods
auto main() -> int {
  "reserveGrow BasicTest"_test = [] -> void {
    ::math::Vector<int> v;
    expect(v.size() == 0);
    auto initial_cap = std::ptrdiff_t(v.getCapacity());
    expect(initial_cap > 0);

    // reserveGrow should grow capacity aggressively (not just to exact size)
    v.reserveGrow(initial_cap + 1);
    auto cap_after_ip1 = std::ptrdiff_t(v.getCapacity());
    expect(ge(cap_after_ip1, 2 * initial_cap + 2));
    expect(v.size() == 0); // size should not change

    // reserveGrow with smaller size should not shrink
    v.reserveGrow(5);
    expect(v.getCapacity() == cap_after_ip1);
    expect(v.size() == 0);

    // reserveGrow with larger size should grow aggressively
    v.reserveGrow(100);
    auto cap_after_100 = v.getCapacity();
    expect(cap_after_100 >= 100);
    expect(cap_after_100 > cap_after_ip1); // should have grown
    expect(v.size() == 0);
  };

  "reserveGrow vs reserve"_test = [] -> void {
    ::math::Vector<int> v1, v2;

    // reserveGrow should allocate more than regular reserve
    v1.reserve(10);
    v2.reserveGrow(10);

    // Both should have at least capacity 10
    expect(v1.getCapacity() >= 10);
    expect(v2.getCapacity() >= 10);

    // reserveGrow should typically allocate more aggressively
    // (though exact behavior depends on newCapacity implementation)
    expect(v2.getCapacity() >= v1.getCapacity());
  };

  "resizeGrow BasicTest"_test = [] -> void {
    ::math::Vector<int> v;
    expect(v.size() == 0);
    auto initial_cap = std::ptrdiff_t(v.getCapacity());
    expect(initial_cap > 0);

    // resizeGrow should both resize AND grow capacity aggressively
    v.resizeGrow(initial_cap + 2);
    auto cap_after_ip1 = std::ptrdiff_t(v.getCapacity());
    expect(ge(cap_after_ip1, 2 * initial_cap + 4));
    expect(eq(v.size(), initial_cap + 2)); // size SHOULD change

    // Verify elements are initialized
    for (std::ptrdiff_t i = 0; i < 10; ++i) v[i] = static_cast<int>(i);

    // resizeGrow with smaller size
    v.resizeGrow(5);
    expect(v.size() == 5);
    expect(v.getCapacity() == cap_after_ip1); // capacity should not shrink

    // Verify remaining elements are intact
    for (std::ptrdiff_t i = 0; i < 5; ++i) expect(v[i] == static_cast<int>(i));

    // resizeGrow with larger size should grow capacity aggressively
    v.resizeGrow(100);
    auto cap_after_100 = v.getCapacity();
    expect(cap_after_100 >= 100);
    expect(cap_after_100 > cap_after_ip1);
    expect(v.size() == 100);

    // Verify old elements are still intact
    for (std::ptrdiff_t i = 0; i < 5; ++i) expect(v[i] == static_cast<int>(i));
  };

  "resizeGrow with push_back"_test = [] -> void {
    ::math::Vector<int> v;

    // Start with resizeGrow to get initial capacity
    v.resizeGrow(5);
    auto cap_after_resize = v.getCapacity();
    expect(v.size() == 5);
    expect(cap_after_resize >= 5);

    // Fill with values
    for (std::ptrdiff_t i = 0; i < 5; ++i) v[i] = static_cast<int>(i * 10);

    // Now resize down
    v.resize(3);
    expect(v.size() == 3);
    expect(v.getCapacity() == cap_after_resize); // capacity unchanged

    // Use resizeGrow to grow again
    v.resizeGrow(20);
    expect(v.size() == 20);
    expect(v.getCapacity() >= 20);

    // Original elements should be preserved
    expect(v[0] == 0);
    expect(v[1] == 10);
    expect(v[2] == 20);
  };

  "reserveGrow reallocation"_test = [] -> void {
    ::math::Vector<double> v;

    // Add some elements
    for (int i = 0; i < 5; ++i) v.push_back(static_cast<double>(i) * 1.5);
    expect(v.size() == 5);

    // Use reserveGrow to ensure capacity for more
    v.reserveGrow(50);
    expect(v.getCapacity() >= 50);
    expect(v.size() == 5); // size unchanged

    // Verify data is intact after potential reallocation
    for (std::ptrdiff_t i = 0; i < 5; ++i)
      expect(v[i] == static_cast<double>(i) * 1.5);

    // Add more elements without reallocation
    for (int i = 5; i < 50; ++i) v.push_back(static_cast<double>(i) * 1.5);
    expect(v.size() == 50);

    // All elements should be correct
    for (std::ptrdiff_t i = 0; i < 50; ++i)
      expect(v[i] == static_cast<double>(i) * 1.5);
  };

  return 0;
}

import Testing;

using namespace testing;
import Array;
import AxisTypes;
import ManagedArray;
import Range;
import SOA;
import std;
import Tuple;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
int main() {
  "SOATest BasicAssertions"_test = [] -> void {
    containers::Tuple x{3, 2.0, 5.0F};
    // static_assert(math::CumSizeOf_v<0, decltype(x)> == 0);
    // static_assert(math::CumSizeOf_v<1, decltype(x)> == 4);
    // static_assert(math::CumSizeOf_v<2, decltype(x)> == 12);
    // static_assert(math::CumSizeOf_v<3, decltype(x)> == 16);
    // math::ManagedSOA soa{std::type_identity<decltype(x)>{}, 5};
    using T = decltype(x);
    static_assert(std::is_trivially_default_constructible_v<
                  containers::Tuple<int, double>>);
    static_assert(std::is_trivially_default_constructible_v<T>);
    static_assert(std::is_trivially_destructible_v<T>);
    ::math::ManagedSOA soa(std::type_identity<decltype(x)>{},
                           ::math::length(5z));
    expect(soa.capacity() == 8);
    static_assert(!std::is_standard_layout_v<decltype(soa)::Base>);
    soa[0] = x;
    soa[1] = {5, 2.25, 5.5F};
    soa.template get<0>(2) = 7;
    soa.template get<1>(2) = 2.5;
    soa.template get<2>(2) = 6.0F;
    soa.template get<0>()[3] = 9;
    soa.template get<1>()[3] = 2.75;
    soa.template get<2>()[3] = 6.5F;
    soa[4] = {11, 3.0, 7.0F};
    {
      std::ptrdiff_t j = 0;
      for (auto [i, d, f] : soa) {
        static_assert(std::same_as<decltype(i), int>);
        static_assert(std::same_as<decltype(d), double>);
        static_assert(std::same_as<decltype(f), float>);
        expect(i == (3 + 2 * j));
        expect(d == (2.0 + 0.25 * j));
        expect(f == (5.0F + 0.5F * j));
        expect(i == soa.get<0>(j));
        expect(d == soa.get<1>(j));
        expect(f == soa.get<2>(j));
        ++j;
      }
    }
    soa.resize(7);
    soa[5] = {13, 3.25, 7.5F};
    soa[6] = {15, 3.5, 8.0F};
    for (std::ptrdiff_t j = 0; j < 7; ++j) {
      decltype(x) y = *(soa.begin() + j);
      // decltype(x) y = soa[j];
      auto [i, d, f] = y;
      static_assert(std::same_as<decltype(i), int>);
      static_assert(std::same_as<decltype(d), double>);
      static_assert(std::same_as<decltype(f), float>);
      expect(i == (3 + 2 * j));
      expect(d == (2.0 + 0.25 * j));
      expect(f == (5.0F + 0.5F * j));
      expect(i == soa.get<0>(j));
      expect(d == soa.get<1>(j));
      expect(f == soa.get<2>(j));
    }
    for (int j = 7; j < 65; ++j) {
      int i = 3 + 2 * j;
      double d = 2.0 + 0.25 * j;
      float f = 5.0F + 0.5F * float(j);
      if (j & 1) soa.emplace_back(i, d, f);
      else soa.push_back({i, d, f});
    }
    for (std::ptrdiff_t j = 0; j < 65; ++j) {
      decltype(x) y = soa[j];
      auto [i, d, f] = y;
      static_assert(std::same_as<decltype(i), int>);
      static_assert(std::same_as<decltype(d), double>);
      static_assert(std::same_as<decltype(f), float>);
      expect(i == (3 + 2 * j));
      expect(d == (2.0 + 0.25 * j));
      expect(f == (5.0F + 0.5F * j));
      expect(i == soa.get<0>(j));
      expect(d == soa.get<1>(j));
      expect(f == soa.get<2>(j));
    }
    expect(soa.size() == 65);
  };

  "SOAPairTest BasicAssertions"_test = [] -> void {
    containers::Pair x{3, 2.0};
    using Tup = decltype(x);
    // static_assert(math::CumSizeOf_v<0, decltype(x)> == 0);
    // static_assert(math::CumSizeOf_v<1, decltype(x)> == 4);
    // static_assert(math::CumSizeOf_v<2, decltype(x)> == 12);
    // math::ManagedSOA soa{std::type_identity<decltype(x)>{}, 5};
    ::math::ManagedSOA<Tup> soa;
    // math::ManagedSOA soa(std::type_identity<decltype(x)>{});
    expect(soa.capacity() == 0);
    soa.push_back(x);
    soa[0] = x;
    expect(soa.capacity() == 8);
    soa.resize(5);
    soa[1] = {5, 2.25};
    soa.template get<0>(2) = 7;
    soa.template get<1>(2) = 2.5;
    soa.template get<0>()[3] = 9;
    soa.template get<1>()[3] = 2.75;
    soa[4] = {11, 3.0};
    auto soa_copy = soa;
    for (std::ptrdiff_t j = 0; j < 5; ++j) {
      Tup y = soa[j];
      auto [i, d] = y;
      static_assert(std::same_as<decltype(i), int>);
      static_assert(std::same_as<decltype(d), double>);
      expect(i == (3 + 2 * j));
      expect(d == (2.0 + 0.25 * j));
      expect(i == soa.get<0>(j));
      expect(d == soa.get<1>(j));
      expect(i == soa_copy.get<0>(j));
      expect(d == soa_copy.get<1>(j));
    }
    soa.resize(7);
    soa[5] = {13, 3.25};
    soa[6] = {15, 3.5};
    for (std::ptrdiff_t j = 0; j < 7; ++j) {
      decltype(x) y = soa[j];
      auto [i, d] = y;
      static_assert(std::same_as<decltype(i), int>);
      static_assert(std::same_as<decltype(d), double>);
      expect(i == (3 + 2 * j));
      expect(d == (2.0 + 0.25 * j));
      expect(i == soa.get<0>(j));
      expect(d == soa.get<1>(j));
    }
    for (int j = 7; j < 65; ++j) {
      int i = 3 + 2 * j;
      double d = 2.0 + 0.25 * j;
      if (j & 1) soa.emplace_back(i, d);
      else soa.push_back({i, d});
    }
    for (std::ptrdiff_t j = 0; j < 65; ++j) {
      decltype(x) y = soa[j];
      auto [i, d] = y;
      static_assert(std::same_as<decltype(i), int>);
      static_assert(std::same_as<decltype(d), double>);
      expect(i == (3 + 2 * j));
      expect(d == (2.0 + 0.25 * j));
      expect(i == soa.get<0>(j));
      expect(d == soa.get<1>(j));
    }
    expect(soa.size() == 65);
  };

  "VecOfSOATest BasicAssertions"_test = [] -> void {
    ::math::Vector<::math::ManagedSOA<containers::Tuple<int, double, float>>>
      vsoa;
    vsoa.emplace_back();
    vsoa.emplace_back();
    vsoa.pop_back();
    vsoa.emplace_back();
    vsoa.emplace_back();
    vsoa.emplace_back();
    expect(vsoa.size() == 4);
  };

  return 0;
}

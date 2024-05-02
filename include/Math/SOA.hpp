#pragma once
#include "Math/Array.hpp"
#include <cstddef>
#ifndef SOA_hpp_INCLUDED
#define SOA_hpp_INCLUDED

#include "Alloc/Mallocator.hpp"
#include "Containers/Pair.hpp"
#include "Containers/Tuple.hpp"
#include "Math/MatrixDimensions.hpp"
#include <memory>
#include <type_traits>
#include <utility>
namespace poly::math {

template <typename... T> struct Types {};

static_assert(std::tuple_size_v<containers::Tuple<int, double>> == 2);
static_assert(std::tuple_size_v<containers::Pair<int, double>> == 2);

template <typename T,
          typename I = std::make_index_sequence<std::tuple_size_v<T>>>
struct TupleTypes {
  using type = void;
};

template <typename T, size_t... I>
struct TupleTypes<T, std::index_sequence<I...>> {
  using type = Types<std::tuple_element_t<I, T>...>;
};

template <typename T> using TupleTypes_t = typename TupleTypes<T>::type;

static_assert(std::is_same_v<TupleTypes_t<containers::Tuple<int, double>>,
                             Types<int, double>>);

template <typename T, typename U = TupleTypes<T>::type>
struct CanConstructFromMembers {
  static constexpr bool value = false;
};
template <typename T, typename... Elts>
struct CanConstructFromMembers<T, Types<Elts...>> {
  static constexpr bool value = std::is_constructible_v<T, Elts...>;
};

static_assert(CanConstructFromMembers<containers::Tuple<int, double>>::value);

namespace CapacityCalculators {

struct Length {
  static constexpr auto operator()(auto sz) -> ptrdiff_t {
    return static_cast<ptrdiff_t>(sz);
  }
  constexpr Length() = default;
  constexpr Length(auto){};
  /// args are old nz, new nz
  /// returns `0` if not growing
  constexpr auto oldnewcap(auto osz, auto nsz) -> std::array<ptrdiff_t, 2> {
    return {(*this)(osz), (*this)(nsz)};
  }
  constexpr auto operator=(ptrdiff_t) -> Length & { return *this; }
};
constexpr auto nextpow2(ptrdiff_t x) -> ptrdiff_t {
  return static_cast<ptrdiff_t>(std::bit_ceil(static_cast<size_t>(x)));
}
struct NextPow2 {
  static constexpr auto operator()(auto sz) -> ptrdiff_t {
    ptrdiff_t i = ptrdiff_t(sz);
    return i ? std::max(8z, nextpow2(i)) : 0z;
  }
  constexpr NextPow2() = default;
  constexpr NextPow2(auto){};
  /// args are old nz, new nz
  /// returns `0` if not growing
  constexpr auto oldnewcap(auto osz, auto nsz) -> std::array<ptrdiff_t, 2> {
    return {(*this)(osz), (*this)(nsz)};
  }
  constexpr auto operator=(ptrdiff_t) -> NextPow2 & { return *this; }
};
struct Explicit {
  ptrdiff_t capacity{0z};
  constexpr Explicit() = default;
  constexpr Explicit(ptrdiff_t sz)
    : capacity{sz ? ((sz > 8z) ? nextpow2(sz) : 8z) : 0z} {}
  constexpr auto operator()(auto) const -> ptrdiff_t { return capacity; }
  constexpr auto operator=(ptrdiff_t cap) -> Explicit & {
    capacity = cap;
    return *this;
  }
  /// args are old nz, new nz
  /// returns `0` if not growing
  constexpr auto oldnewcap(auto, auto sz) -> std::array<ptrdiff_t, 2> {
    return {capacity, nextpow2(sz)};
  }
};

} // namespace CapacityCalculators

template <size_t I, typename T> struct CumSizeOf {
  static constexpr size_t value = 0;
};

template <size_t I, typename T, typename... Ts>
struct CumSizeOf<I, Types<T, Ts...>> {
  static constexpr size_t value =
    sizeof(T) + CumSizeOf<I - 1, Types<Ts...>>::value;
};
template <typename T, typename... Ts> struct CumSizeOf<0, Types<T, Ts...>> {
  static constexpr size_t value = 0;
};
template <> struct CumSizeOf<0, Types<>> {
  static constexpr size_t value = 0;
};
template <size_t I, typename T>
inline constexpr size_t CumSizeOf_v = CumSizeOf<I, TupleTypes_t<T>>::value;

// std::conditional_t<MatrixDimension<S>, CapacityCalculators::Length,
//                                CapacityCalculators::NextPow2>
/// requires 16 byte alignment of allocated pointer
template <typename T, typename S = ptrdiff_t,
          typename C = CapacityCalculators::Explicit,
          typename TT = TupleTypes_t<T>,
          typename II = std::make_index_sequence<std::tuple_size_v<T>>>
struct SOA {};
template <typename T, typename S, typename C, typename... Elts, size_t... II>
requires(CanConstructFromMembers<T>::value)
struct SOA<T, S, C, Types<Elts...>, std::index_sequence<II...>> {
  char *data;
  [[no_unique_address]] S sz;
  [[no_unique_address]] C capacity;
  static constexpr bool trivial =
    std::is_trivially_default_constructible_v<T> &&
    std::is_trivially_destructible_v<T>;
  static_assert(trivial);
  struct Reference {
    char *ptr;
    ptrdiff_t stride;
    ptrdiff_t i;
    operator T() const {
      char *p = std::assume_aligned<16>(ptr);
      return T(*reinterpret_cast<std::tuple_element_t<II, T> *>(
        p + CumSizeOf_v<II, T> * stride +
        sizeof(std::tuple_element_t<II, T>) * i)...);
    }
    template <size_t I> void assign(const T &x) {
      char *p = std::assume_aligned<16>(ptr);
      *reinterpret_cast<std::tuple_element_t<I, T> *>(
        p + CumSizeOf_v<I, T> * stride +
        sizeof(std::tuple_element_t<I, T>) * i) = x.template get<I>();
    }
    auto operator=(const T &x) -> Reference & {
      ((void)assign<II>(x), ...);
      // assign<II...>(x);
      return *this;
    }
    auto operator=(Reference x) -> Reference & {
      (*this) = T(x);
      return *this;
    }
    template <size_t I> auto get() -> std::tuple_element_t<I, T> & {
      invariant(i >= 0);
      return *reinterpret_cast<std::tuple_element_t<I, T> *>(
        ptr + CumSizeOf_v<I, T> * stride +
        sizeof(std::tuple_element_t<I, T>) * i);
    }
    template <size_t I> auto get() const -> const std::tuple_element_t<I, T> & {
      invariant(i >= 0);
      return *reinterpret_cast<const std::tuple_element_t<I, T> *>(
        ptr + CumSizeOf_v<I, T> * stride +
        sizeof(std::tuple_element_t<I, T>) * i);
    }
  };
  auto operator[](ptrdiff_t i) const -> T {
    char *p = std::assume_aligned<16>(data);
    ptrdiff_t stride = capacity(sz);
    return T(*reinterpret_cast<std::tuple_element_t<II, T> *>(
      p + CumSizeOf_v<II, T> * stride +
      sizeof(std::tuple_element_t<II, T>) * i)...);
  }
  auto operator[](ptrdiff_t i) -> Reference { return {data, capacity(sz), i}; }
  static constexpr auto totalSizePer() -> size_t {
    return CumSizeOf_v<sizeof...(II), T>;
  }
  [[nodiscard]] constexpr auto size() const -> ptrdiff_t { return sz; }
  template <size_t I> auto get(ptrdiff_t i) -> std::tuple_element_t<I, T> & {
    invariant(i >= 0);
    invariant(i < size());
    return *reinterpret_cast<std::tuple_element_t<I, T> *>(
      data + CumSizeOf_v<I, T> * capacity(sz) +
      sizeof(std::tuple_element_t<I, T>) * i);
  }
  template <size_t I>
  auto get(ptrdiff_t i) const -> const std::tuple_element_t<I, T> & {
    invariant(i >= 0);
    invariant(i < size());
    return *reinterpret_cast<const std::tuple_element_t<I, T> *>(
      data + CumSizeOf_v<I, T> * capacity(sz) +
      sizeof(std::tuple_element_t<I, T>) * i);
  }
  template <size_t I>
  auto get() -> math::MutArray<std::tuple_element_t<I, T>, S> {
    return {reinterpret_cast<std::tuple_element_t<I, T> *>(
              data + CumSizeOf_v<I, T> * capacity(sz)),
            sz};
  }
  template <size_t I>
  auto get() const -> math::Array<std::tuple_element_t<I, T>, S> {
    return {reinterpret_cast<const std::tuple_element_t<I, T> *>(
              data + CumSizeOf_v<I, T> * capacity(sz)),
            sz};
  }
};

template <typename T, typename S = ptrdiff_t,
          typename C = CapacityCalculators::Explicit,
          typename TT = TupleTypes_t<T>,
          typename II = std::make_index_sequence<std::tuple_size_v<T>>,
          class A = alloc::Mallocator<char>>
struct ManagedSOA : public SOA<T, S, C, TT, II> {
  using Base = SOA<T, S, C, TT, II>;
  /// uninitialized allocation
  constexpr ManagedSOA() = default;
  ManagedSOA(S nsz) {
    this->sz = nsz;
    this->capacity = C{nsz};
    ptrdiff_t stride = this->capacity(this->sz);
    this->data = alloc(stride);
    // this->data =
    //   A::allocate(stride * this->totalSizePer(), std::align_val_t{16});
  }
  constexpr ManagedSOA(std::type_identity<T>) : ManagedSOA() {}
  ManagedSOA(std::type_identity<T>, S nsz) : ManagedSOA(nsz) {}
  ManagedSOA(const ManagedSOA &other) {
    this->sz = other.sz;
    this->capacity = C{this->sz};
    this->data = alloc(this->capacity(this->sz));
  }

  void free(ptrdiff_t stride) {
    A::deallocate(this->data, stride * Base::totalSizePer());
  }
  static auto alloc(ptrdiff_t stride) -> char * {
    return A::allocate(stride * Base::totalSizePer());
  }
  ~ManagedSOA() {
    if (this->data) free(this->capacity(this->sz));
  }
  constexpr ManagedSOA(ManagedSOA &&other)
    : SOA<T, S, C, TT, II>{other.data, other.sz, other.capacity} {
    other.data = nullptr;
    other.sz = {};
    other.capacity = {};
  }
  constexpr auto operator=(ManagedSOA &&other) -> ManagedSOA & {
    std::swap(this->data, other.data);
    std::swap(this->sz, other.sz);
    std::swap(this->capacity, other.capacity);
    return *this;
  }
  constexpr auto operator=(const ManagedSOA &other) -> ManagedSOA & {
    if (this == &other) return *this;
    auto L = static_cast<ptrdiff_t>(other.sz);
    resizeForOverwrite(other.sz);
    ManagedSOA &self{*this};
    for (ptrdiff_t i = 0z; i < L; ++i) self[i] = other[i];
    return *this;
  }
  void resizeForOverwrite(S nsz) {
    auto [ocap, ncap] = this->capacity.oldnewcap(this->sz, nsz);
    this->sz = nsz;
    if (ocap >= ncap) return;
    if (this->data) free(ocap);
    this->data = A::allocate(ncap * this->totalSizePer());
    this->capacity = ncap;
  }
  void resize(S nsz) {
    auto [ocap, ncap] = this->capacity.oldnewcap(this->sz, nsz);
    if (ocap >= ncap) {
      this->sz = nsz;
      return;
    }
    ManagedSOA other{nsz};
    ManagedSOA &self{*this};
    std::swap(self.data, other.data);
    std::swap(self.sz, other.sz);
    std::swap(self.capacity, other.capacity);
    // FIXME: only accepts non-copyable
    for (ptrdiff_t i = 0, L = std::min(ptrdiff_t(self.sz), ptrdiff_t(other.sz));
         i < L; ++i)
      self[i] = other[i];
  }
  constexpr void clear() { this->sz = {}; }
  template <typename... Args> void emplace_back(Args &&...args) {
    push_back(T(args...));
  }
  void push_back(T arg) {
    S osz = this->sz;
    resize(osz + 1);
    (*this)[osz] = arg;
  }
};
template <typename T, typename S>
ManagedSOA(std::type_identity<T>, S) -> ManagedSOA<T, S>;

} // namespace poly::math
#endif // SOA_hpp_INCLUDED

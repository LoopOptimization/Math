module;
#include <cstddef>
#include <memory>

export module ArrayConstructors;

import Arena;
import Array;
import AxisTypes;
import MatDim;

using utils::eltype_t;
export namespace math {

template <alloc::Allocator A, typename T>
using rebound_alloc =
  typename std::allocator_traits<A>::template rebind_alloc<T>;

template <class T, alloc::FreeAllocator A>
constexpr auto vector(A a, ptrdiff_t M) {
  if constexpr (std::same_as<T, eltype_t<A>>) return vector(a, M);
  else return vector(rebound_alloc<A, T>{}, M);
}
template <class T>
constexpr auto vector(WArena<T> alloc,
                      ptrdiff_t M) -> ResizeableView<T, Length<>> {
  return {alloc.allocate(M), length(M), capacity(M)};
}
template <class T, size_t SlabSize, bool BumpUp>
constexpr auto vector(Arena<SlabSize, BumpUp> *alloc,
                      ptrdiff_t M) -> ResizeableView<T, Length<>> {
  return {alloc->template allocate<compressed_t<T>>(M), length(M), capacity(M)};
}

template <typename T, alloc::FreeAllocator A>
constexpr auto vector(A a, ptrdiff_t M, std::type_identity_t<T> x) {
  if constexpr (std::same_as<T, eltype_t<A>>) return vector(a, M, x);
  else return vector(rebound_alloc<A, T>{}, M, x);
}

template <class T>
constexpr auto vector(WArena<T> alloc, ptrdiff_t M,
                      T x) -> ResizeableView<T, Length<>> {
  ResizeableView<T, Length<>> a{alloc.allocate(M), length(M), capacity(M)};
  a.fill(x);
  return a;
}
template <class T, size_t SlabSize, bool BumpUp>
constexpr auto vector(Arena<SlabSize, BumpUp> *alloc, ptrdiff_t M,
                      T x) -> ResizeableView<T, Length<>> {
  ResizeableView<T, Length<>> a{alloc->template allocate<compressed_t<T>>(M),
                                length(M), capacity(M)};
  a.fill(x);
  return a;
}

template <class T, alloc::FreeAllocator A>
constexpr auto square_matrix(A a, ptrdiff_t M) {
  if constexpr (std::same_as<T, eltype_t<A>>) return square_matrix(a, M);
  else return square_matrix(rebound_alloc<A, T>{}, M);
}

template <class T>
constexpr auto square_matrix(WArena<T> alloc,
                             ptrdiff_t M) -> MutSquarePtrMatrix<T> {
  return {alloc.allocate(M * M), SquareDims<>{row(M)}};
}
template <class T, size_t SlabSize, bool BumpUp>
constexpr auto square_matrix(Arena<SlabSize, BumpUp> *alloc,
                             ptrdiff_t M) -> MutSquarePtrMatrix<T> {
  return {alloc->template allocate<compressed_t<T>>(M * M),
          SquareDims<>{row(M)}};
}
template <class T, alloc::FreeAllocator A>
constexpr auto square_matrix(A a, ptrdiff_t M, std::type_identity_t<T> x) {
  if constexpr (std::same_as<T, eltype_t<A>>) return square_matrix(a, M, x);
  else return square_matrix(rebound_alloc<A, T>{}, M, x);
}
template <class T>
constexpr auto square_matrix(WArena<T> alloc, ptrdiff_t M,
                             T x) -> MutSquarePtrMatrix<T> {
  MutSquarePtrMatrix<T> A{alloc.allocate(M * M), SquareDims<>{row(M)}};
  A.fill(x);
  return A;
}
template <class T, size_t SlabSize, bool BumpUp>
constexpr auto square_matrix(Arena<SlabSize, BumpUp> *alloc, ptrdiff_t M,
                             T x) -> MutSquarePtrMatrix<T> {
  MutSquarePtrMatrix<T> A{alloc->template allocate<compressed_t<T>>(M * M),
                          SquareDims<>{row(M)}};
  A.fill(x);
  return A;
}

template <class T, alloc::FreeAllocator A, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(A a, Row<R> M, Col<C> N) {
  if constexpr (std::same_as<T, eltype_t<A>>) return matrix(a, M, N);
  else return matrix(rebound_alloc<A, T>{}, M, N);
}
template <class T, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(WArena<T> alloc, Row<R> M,
                      Col<C> N) -> MutArray<T, DenseDims<R, C>> {
  auto memamt = size_t(ptrdiff_t(M) * ptrdiff_t(N));
  return {alloc.allocate(memamt), DenseDims{M, N}};
}
template <class T, size_t SlabSize, bool BumpUp, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(Arena<SlabSize, BumpUp> *alloc, Row<R> M,
                      Col<C> N) -> MutArray<T, DenseDims<R, C>> {
  auto memamt = size_t(ptrdiff_t(M) * ptrdiff_t(N));
  return {alloc->template allocate<compressed_t<T>>(memamt), M, N};
}
template <class T, size_t SlabSize, bool BumpUp, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(Arena<SlabSize, BumpUp> *alloc, CartesianIndex<R, C> dim)
  -> MutArray<T, DenseDims<R, C>> {
  Row M = Row(dim);
  Col N = Col(dim);
  auto memamt = size_t(ptrdiff_t(M) * ptrdiff_t(N));
  return {alloc->template allocate<compressed_t<T>>(memamt), M, N};
}

template <class T, alloc::FreeAllocator A, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(A a, Row<R> M, Col<C> N, std::type_identity_t<T> x) {
  if constexpr (std::same_as<T, eltype_t<A>>) return square_matrix(a, M, N, x);
  else return square_matrix(rebound_alloc<A, T>{}, M, N, x);
}
template <class T, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(WArena<T> alloc, Row<R> M, Col<C> N,
                      T x) -> MutArray<T, DenseDims<R, C>> {
  auto memamt = size_t(ptrdiff_t(M) * ptrdiff_t(N));
  MutArray<T, DenseDims<R, C>> A{alloc.allocate(memamt), DenseDims{M, N}};
  A.fill(x);
  return A;
}
template <class T, size_t SlabSize, bool BumpUp, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(Arena<SlabSize, BumpUp> *alloc, Row<R> M, Col<C> N,
                      T x) -> MutArray<T, DenseDims<R, C>> {
  auto memamt = size_t(ptrdiff_t(M) * ptrdiff_t(N));
  MutArray<T, DenseDims<R, C>> A{
    alloc->template allocate<compressed_t<T>>(memamt), DenseDims{M, N}};
  A.fill(x);
  return A;
}

template <class T, alloc::FreeAllocator A>
constexpr auto identity(A a, ptrdiff_t M) {
  if constexpr (std::same_as<T, eltype_t<A>>) return identity(a, M);
  else return identity(rebound_alloc<A, T>{}, M);
}
template <class T>
constexpr auto identity(WArena<T> alloc, ptrdiff_t M) -> MutSquarePtrMatrix<T> {
  MutSquarePtrMatrix<T> A{square_matrix(alloc, M, T{})};
  A.diag() << T{1};
  return A;
}
template <class T, size_t SlabSize, bool BumpUp>
constexpr auto identity(Arena<SlabSize, BumpUp> *alloc,
                        ptrdiff_t M) -> MutSquarePtrMatrix<T> {
  MutSquarePtrMatrix<T> A{square_matrix(alloc, M, T{})};
  A.diag() << T{1};
  return A;
}

template <typename T, typename I>
concept Alloc = requires(T t, ptrdiff_t M, Row<> r, Col<> c, I i) {
  { identity<I>(t, M) } -> std::convertible_to<MutSquarePtrMatrix<I>>;
  { square_matrix<I>(t, M) } -> std::convertible_to<MutSquarePtrMatrix<I>>;
  { square_matrix<I>(t, M, i) } -> std::convertible_to<MutSquarePtrMatrix<I>>;
  { matrix<I>(t, r, c) } -> std::convertible_to<MutDensePtrMatrix<I>>;
  { matrix(t, r, c, i) } -> std::convertible_to<MutDensePtrMatrix<I>>;
  { vector<I>(t, M) } -> std::convertible_to<MutPtrVector<I>>;
};
static_assert(Alloc<std::allocator<int64_t>, int64_t>);
static_assert(Alloc<alloc::Mallocator<int64_t>, int64_t>);
static_assert(Alloc<Arena<> *, int64_t>);
} // namespace math

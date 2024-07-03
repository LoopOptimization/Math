module;

#include <concepts>
#include <memory>

export module GenericArrayConstructors;

export import ArrayConstructors;
import ManagedArray;

export namespace math {

using alloc::Arena, alloc::WArena, alloc::OwningArena, utils::eltype_t;

template <alloc::Allocator A, typename T>
using rebound_alloc =
  typename std::allocator_traits<A>::template rebind_alloc<T>;

// template <alloc::Allocator A>
// using alloc_type = typename std::allocator_traits<A>::value_type;

template <alloc::FreeAllocator A>
constexpr auto vector(A a, ptrdiff_t M)
  -> ManagedArray<eltype_t<A>, Length<>,
                  containers::PreAllocStorage<eltype_t<A>, ptrdiff_t>(), A> {
  return {length(M), a};
}
template <alloc::FreeAllocator A>
constexpr auto vector(A a, ptrdiff_t M, eltype_t<A> x)
  -> ManagedArray<eltype_t<A>, Length<>,
                  containers::PreAllocStorage<eltype_t<A>, ptrdiff_t>(), A> {
  return {length(M), x, a};
}

template <class T, alloc::FreeAllocator A>
constexpr auto vector(A a, ptrdiff_t M) {
  if constexpr (std::same_as<T, eltype_t<A>>) return vector(a, M);
  else return vector(rebound_alloc<A, T>{}, M);
}
template <typename T, alloc::FreeAllocator A>
constexpr auto vector(A a, ptrdiff_t M, std::type_identity_t<T> x) {
  if constexpr (std::same_as<T, eltype_t<A>>) return vector<A>(a, M, x);
  else return vector<A>(rebound_alloc<A, T>{}, M, x);
}

template <alloc::FreeAllocator A, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(A a, Row<R> M, Col<C> N) -> DenseMatrixAlloc<A, R, C> {
  return {DenseDims{M, N}, a};
}
template <alloc::FreeAllocator A, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(A a, Row<R> M, Col<C> N,
                      eltype_t<A> x) -> DenseMatrixAlloc<A, R, C> {
  return {DenseDims{M, N}, x, a};
}

template <class T, alloc::FreeAllocator A, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(A a, Row<R> M, Col<C> N) -> DenseMatrixAlloc<A, R, C> {
  return {DenseDims{M, N}, rebound_alloc<A, T>{}};
}
template <class T, alloc::FreeAllocator A, ptrdiff_t R, ptrdiff_t C>
constexpr auto matrix(A a, Row<R> M, Col<C> N,
                      std::type_identity_t<T> x) -> DenseMatrixAlloc<A, R, C> {
  return {DenseDims{M, N}, x, rebound_alloc<A, T>{}};
}

template <alloc::FreeAllocator A>
constexpr auto square_matrix(A a, ptrdiff_t M) -> SquareMatrixAlloc<A> {
  return {SquareDims<>{row(M)}, A{}};
}
template <alloc::FreeAllocator A>
constexpr auto square_matrix(A a, ptrdiff_t M, std::type_identity_t<T> x)
  -> SquareMatrixAlloc<A> {
  return {SquareDims<>{row(M)}, x, A{}};
}
template <class T, alloc::FreeAllocator A>
constexpr auto
square_matrix(A a, ptrdiff_t M) -> SquareMatrixAlloc<rebound_alloc<A, T>> {
  return {SquareDims<>{row(M)}, rebound_alloc<A, T>{}};
}
template <class T, alloc::FreeAllocator A>
constexpr auto square_matrix(A a, ptrdiff_t M, std::type_identity_t<T> x)
  -> SquareMatrixAlloc<rebound_alloc<A, T>> {
  return {SquareDims<>{row(M)}, x, rebound_alloc<A, T>{}};
}

template <alloc::FreeAllocator A>
constexpr auto identity(A a, ptrdiff_t M) -> SquareMatrixAlloc<A::value_type> {
  SquareMatrixAlloc<A> B{SquareDims<>{row(M)}, A::value_type{}, a};
  B.diag() << eltype_t<A>{1};
  return B;
}
template <class T, alloc::FreeAllocator A>
constexpr auto
identity(A a, ptrdiff_t M) -> SquareMatrixAlloc<T, rebound_alloc<A, T>> {
  return identity<rebound_alloc<A, T>>(rebound_alloc<A, T>{}, M);
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
static_assert(Alloc<alloc::Arena<> *, int64_t>);

} // namespace math

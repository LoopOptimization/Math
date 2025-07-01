#ifdef USE_MODULE
module;
#else
#pragma once
#endif

#ifndef USE_MODULE
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "Math/Array.cxx"
#include "Math/AxisTypes.cxx"
#include "Math/ManagedArray.cxx"
#include "Math/MatrixDimensions.cxx"
#include "Utilities/ArrayPrint.cxx"
#else
export module SmallSparseMatrix;

import Array;
import AxisTypes;
import CorePrint;
import ManagedArray;
import MatDim;
import std;
#endif

#ifdef USE_MODULE
export namespace math {
#else
namespace math {
#endif
// this file is not used at the moment
template <typename T> class SmallSparseMatrix {
  // non-zeros
  [[no_unique_address]] Vector<T> nonZeros{};
  // masks, the upper 8 bits give the number of elements in previous rows
  // the remaining 24 bits are a mask indicating non-zeros within this row
  static constexpr std::ptrdiff_t maxElemPerRow = 24;
  [[no_unique_address]] Vector<std::uint32_t> rows;
  [[no_unique_address]] Col<> ncol;

public:
  [[nodiscard]] constexpr auto getNonZeros() const -> PtrVector<T> {
    return nonZeros;
  }
  [[nodiscard]] constexpr auto getRows() const -> PtrVector<std::uint32_t> {
    return rows;
  }

  [[nodiscard]] constexpr auto numRow() const -> Row<> {
    return row(rows.size());
  }
  [[nodiscard]] constexpr auto numCol() const -> Col<> { return ncol; }
  [[nodiscard]] constexpr auto size() const
    -> CartesianIndex<std::ptrdiff_t, std::ptrdiff_t> {
    return {numRow(), numCol()};
  }
  [[nodiscard]] constexpr auto dim() const -> DenseDims<> {
    return {numRow(), numCol()};
  }
  // [[nodiscard]] constexpr auto view() const -> auto & { return *this; };
  constexpr SmallSparseMatrix(Row<> numRows, Col<> numCols)
    : rows(length(std::ptrdiff_t(numRows)), 0), ncol{numCols} {
    invariant(std::ptrdiff_t(ncol) <= maxElemPerRow);
  }
  constexpr auto get(Row<> i, Col<> j) const -> T {
    invariant(j < ncol);
    std::uint32_t r(rows[std::ptrdiff_t(i)]);
    std::uint32_t jshift = std::uint32_t(1) << std::uint32_t(std::ptrdiff_t(j));
    if (!(r & jshift)) return T{};
    // offset from previous rows
    std::uint32_t prevRowOffset = r >> maxElemPerRow;
    std::uint32_t rowOffset = std::popcount(r & (jshift - 1));
    return nonZeros[rowOffset + prevRowOffset];
  }
  constexpr auto operator[](std::ptrdiff_t i, std::ptrdiff_t j) const -> T {
    return get(row(i), col(j));
  }
  constexpr void insert(T x, Row<> i, Col<> j) {
    invariant(j < ncol);
    std::uint32_t r{rows[std::ptrdiff_t(i)]};
    std::uint32_t jshift = std::uint32_t(1) << std::ptrdiff_t(j);
    // offset from previous rows
    std::uint32_t prevRowOffset = r >> maxElemPerRow;
    std::uint32_t rowOffset = std::popcount(r & (jshift - 1));
    std::ptrdiff_t k = rowOffset + prevRowOffset;
    if (r & jshift) {
      nonZeros[k] = std::move(x);
    } else {
      nonZeros.insert(nonZeros.begin() + k, std::move(x));
      rows[std::ptrdiff_t(i)] = r | jshift;
      for (std::ptrdiff_t l = std::ptrdiff_t(i) + 1; l < rows.size(); ++l)
        rows[l] += std::uint32_t(1) << maxElemPerRow;
    }
  }

  struct Reference {
    [[no_unique_address]] SmallSparseMatrix<T> *A;
    [[no_unique_address]] std::ptrdiff_t i, j;
    constexpr operator T() const { return A->get(row(i), col(j)); }
    constexpr auto operator=(T x) -> Reference & {
      A->insert(std::move(x), row(i), col(j));
      return *this;
    }
  };
  constexpr auto operator[](std::ptrdiff_t i, std::ptrdiff_t j) -> Reference {
    return Reference{this, i, j};
  }
  template <std::convertible_to<T> Y, MatrixDimension S, std::ptrdiff_t L,
            typename A>
  operator ManagedArray<Y, S, L, A>() const {
    ManagedArray<Y, S, L, A> B(dim(), 0);
    std::ptrdiff_t k = 0;
    for (std::ptrdiff_t i = 0; i < numRow(); ++i) {
      std::uint32_t m = getRows()[i] & 0x00ffffff;
      std::ptrdiff_t j = 0;
      while (m) {
        std::uint32_t tz = std::countr_zero(m);
        m >>= tz + 1;
        j += tz;
        B[i, j++] = T(getNonZeros()[k++]);
      }
    }
    invariant(k == getNonZeros().size());
    return B;
  }

  void print() const {
    std::ptrdiff_t k = 0;
    utils::print("[ ");
    for (std::ptrdiff_t i = 0; i < numRow(); ++i) {
      if (i) utils::print("  ");
      std::uint32_t m = rows[i] & 0x00ffffff;
      std::ptrdiff_t j = 0;
      while (m) {
        if (j) utils::print(" ");
        std::uint32_t tz = std::countr_zero(m);
        m >>= (tz + 1);
        j += (tz + 1);
        while (tz--) utils::print(" 0 ");
        const T &x = nonZeros[k++];
        if (x >= 0) utils::print(" ");
        utils::print(x);
      }
      for (; j < numCol(); ++j) utils::print("  0");
      utils::print('\n');
    }
    utils::print(" ]");
    invariant(k == nonZeros.size());
  }

private:
  template <std::convertible_to<T> Y, MatrixDimension S>
  [[gnu::flatten]] friend constexpr auto operator<<(MutArray<Y, S> A,
                                                    const SmallSparseMatrix &B)
    -> MutArray<Y, S> {
    std::ptrdiff_t M = std::ptrdiff_t(A.numRow()),
                   N = std::ptrdiff_t(A.numCol()),
                   X = std::ptrdiff_t(A.rowStride()), k = 0;
    invariant(M, std::ptrdiff_t(B.numRow()));
    invariant(N, std::ptrdiff_t(B.numCol()));
    T *mem = A.data();
    PtrVector<T> nz = B.getNonZeros();
    PtrVector<std::uint32_t> rws = B.getRows();
    for (std::ptrdiff_t i = 0; i < M; ++i) {
      std::uint32_t m = rws[i] & 0x00ffffff;
      std::ptrdiff_t j = 0, l = X * i;
      while (m) {
        std::uint32_t tz = std::countr_zero(m);
        m >>= tz + 1;
        for (; tz; --tz) mem[l + j++] = T{};
        mem[l + j++] = nz[k++];
      }
      for (; j < N; ++j) mem[l + j] = T{};
    }
    invariant(k == nz.size());
    return A;
  }
};

}; // namespace math

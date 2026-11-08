/* Sparse Linear Algebra
Version 0.1.0
https://github.com/FrancoisCarouge/SparseLinearAlgebra

SPDX-License-Identifier: Unlicense

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org> */

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_PATTERN_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_PATTERN_HPP

//! @file
//! @brief The sparsity patterns: the positions of the stored elements of a
//! matrix, their construction, and the named patterns.

#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <cstddef>
#include <limits>

namespace fcarouge {
namespace sparse_linear_algebra::internal {
//! @brief Called by a pattern's invalid construction: never defined, its call
//! ends the constant evaluation, and is the name the compiler's diagnostic
//! reports.
void position_out_of_range_or_duplicate();

//! @brief The predicate of every position.
inline constexpr auto every{
    [](std::size_t /*row*/, std::size_t /*column*/) -> bool { return true; }};

//! @brief The predicate of the positions of the `Lower` sub-diagonals, the main
//! diagonal, and the `Upper` super-diagonals.
template <std::size_t Lower, std::size_t Upper>
inline constexpr auto within_band{
    [](std::size_t row, std::size_t column) -> bool {
      return row <= column + Lower && column <= row + Upper;
    }};
} // namespace sparse_linear_algebra::internal

namespace sparse {
//! @name Patterns
//! @{

//! @brief The index of a position a pattern does not store.
inline constexpr std::size_t absent{std::numeric_limits<std::size_t>::max()};

//! @brief The row and column of a stored element.
struct position {
  //! @brief The row of the element.
  std::size_t row{0};

  //! @brief The column of the element.
  std::size_t column{0};

  //! @brief The origin position.
  constexpr position() = default;

  //! @brief The position of the given row and column.
  constexpr position(std::size_t row_index, std::size_t column_index)
      : row{row_index}, column{column_index} {}

  //! @brief Whether two positions are the same.
  [[nodiscard]] friend constexpr auto operator==(const position &lhs,
                                                 const position &rhs)
      -> bool = default;

  //! @brief The row-major order.
  [[nodiscard]] friend constexpr auto operator<(const position &lhs,
                                                const position &rhs) -> bool {
    return lhs.row < rhs.row || (lhs.row == rhs.row && lhs.column < rhs.column);
  }
};

//! @brief The positions of the stored elements of a matrix.
//!
//! @details A structural type, usable as a template argument, of strictly
//! ascending row-major positions: two equal patterns are one and the same type
//! of matrix, and a stored element's index is its rank among the positions.
//! Prefer `make_pattern`, `generate`, and the named patterns, which establish
//! that order.
//!
//! @tparam Rows The count of rows of the matrix.
//! @tparam Columns The count of columns of the matrix.
//! @tparam Count The count of stored elements of the matrix.
template <std::size_t Rows, std::size_t Columns, std::size_t Count>
struct pattern {
  //! @brief The count of rows of the matrix.
  std::size_t rows{Rows};

  //! @brief The count of columns of the matrix.
  std::size_t columns{Columns};

  //! @brief The count of stored elements of the matrix.
  std::size_t count{Count};

  //! @brief The stored positions, in row-major order.
  //!
  //! @details One unused position when none is stored: not every standard
  //! library's empty array is a structural type.
  std::array<position, (Count > 0 ? Count : 1)> positions{};

  //! @brief The index of the first stored position not before the given one.
  //!
  //! @complexity Logarithmic in the count of stored elements.
  [[nodiscard]] constexpr auto lower_bound(position value) const
      -> std::size_t {
    std::size_t first{0};
    std::size_t size{Count};

    while (size > 0) {
      const std::size_t half{size / 2};

      if (positions[first + half] < value) {
        first += half + 1;
        size -= half + 1;
      } else {
        size = half;
      }
    }

    return first;
  }

  //! @brief The index of the stored element at the given row and column, or
  //! `absent` for a structural zero.
  //!
  //! @complexity Logarithmic in the count of stored elements.
  [[nodiscard]] constexpr auto index_of(std::size_t row,
                                        std::size_t column) const
      -> std::size_t {
    const position value{row, column};
    const std::size_t index{lower_bound(value)};

    return index < Count && positions[index] == value ? index : absent;
  }

  //! @brief Whether the given row and column is a stored element.
  //!
  //! @complexity Logarithmic in the count of stored elements.
  [[nodiscard]] constexpr auto stores(std::size_t row, std::size_t column) const
      -> bool {
    return index_of(row, column) != absent;
  }
};

//! @brief Whether two patterns, of any sizes, store the same positions of the
//! same shape.
template <std::size_t LhsRows, std::size_t LhsColumns, std::size_t LhsCount,
          std::size_t RhsRows, std::size_t RhsColumns, std::size_t RhsCount>
[[nodiscard]] constexpr auto
operator==(const pattern<LhsRows, LhsColumns, LhsCount> &lhs,
           const pattern<RhsRows, RhsColumns, RhsCount> &rhs) -> bool {
  if constexpr (LhsRows != RhsRows || LhsColumns != RhsColumns ||
                LhsCount != RhsCount) {
    return false;
  } else {
    return lhs.positions == rhs.positions;
  }
}

//! @brief The pattern of a matrix without any stored element.
template <std::size_t Rows, std::size_t Columns>
[[nodiscard]] consteval auto make_pattern() -> pattern<Rows, Columns, 0> {
  return {};
}

//! @brief The pattern of the given positions, in any order.
//!
//! @details A position out of the matrix, or a duplicate position, does not
//! compile.
//!
//! @tparam Rows The count of rows of the matrix.
//! @tparam Columns The count of columns of the matrix.
//! @tparam Count The count of the given positions.
//!
//! @param positions The stored positions, in any order.
template <std::size_t Rows, std::size_t Columns, std::size_t Count>
[[nodiscard]] consteval auto make_pattern(const position (&positions)[Count])
    -> pattern<Rows, Columns, Count> {
  pattern<Rows, Columns, Count> result;

  for (std::size_t i{0}; i < Count; ++i) {
    const position value{positions[i]};
    std::size_t j{i};

    if (value.row >= Rows || value.column >= Columns) {
      sparse_linear_algebra::internal::position_out_of_range_or_duplicate();
    }

    for (; j > 0 && value < result.positions[j - 1]; --j) {
      result.positions[j] = result.positions[j - 1];
    }

    if (j > 0 && value == result.positions[j - 1]) {
      sparse_linear_algebra::internal::position_out_of_range_or_duplicate();
    }

    result.positions[j] = value;
  }

  return result;
}

//! @brief The pattern of the positions satisfying a predicate.
//!
//! @tparam Rows The count of rows of the matrix.
//! @tparam Columns The count of columns of the matrix.
//! @tparam Predicate A stateless predicate of a row and a column.
template <std::size_t Rows, std::size_t Columns, auto Predicate>
[[nodiscard]] consteval auto generate() {
  constexpr std::size_t count{[] -> std::size_t {
    std::size_t result{0};

    for (std::size_t row{0}; row < Rows; ++row) {
      for (std::size_t column{0}; column < Columns; ++column) {
        result += Predicate(row, column) ? std::size_t{1} : std::size_t{0};
      }
    }

    return result;
  }()};

  pattern<Rows, Columns, count> result;

  for (std::size_t row{0}, index{0}; row < Rows; ++row) {
    for (std::size_t column{0}; column < Columns; ++column) {
      if (Predicate(row, column)) {
        result.positions[index] = {row, column};
        ++index;
      }
    }
  }

  return result;
}

//! @brief Every position is stored: a matrix without structural zeros.
template <std::size_t Rows, std::size_t Columns>
inline constexpr auto dense{
    generate<Rows, Columns, sparse_linear_algebra::internal::every>()};

//! @brief No position is stored: the zero matrix.
template <std::size_t Rows, std::size_t Columns>
inline constexpr auto zero{make_pattern<Rows, Columns>()};

//! @brief The positions of the `Lower` sub-diagonals, the main diagonal, and
//! the `Upper` super-diagonals.
template <std::size_t Rows, std::size_t Columns, std::size_t Lower,
          std::size_t Upper>
inline constexpr auto banded{
    generate<Rows, Columns,
             sparse_linear_algebra::internal::within_band<Lower, Upper>>()};

//! @brief The main diagonal of a square matrix.
template <std::size_t Size>
inline constexpr auto diagonal{banded<Size, Size, 0, 0>};

//! @brief The three central diagonals of a square matrix.
template <std::size_t Size>
inline constexpr auto tridiagonal{banded<Size, Size, 1, 1>};

//! @brief The main diagonal and below of a square matrix.
template <std::size_t Size>
inline constexpr auto lower_triangular{banded<Size, Size, Size, 0>};

//! @brief The main diagonal and above of a square matrix.
template <std::size_t Size>
inline constexpr auto upper_triangular{banded<Size, Size, 0, Size>};

//! @}
} // namespace sparse
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_PATTERN_HPP

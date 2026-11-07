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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_PATTERN_ALGEBRA_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_PATTERN_ALGEBRA_HPP

//! @file
//! @brief The pattern algebra: the compile-time computations the operations
//! derive the patterns of their results with.

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra/pattern.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <cstddef>

namespace fcarouge::sparse {
//! @name Pattern Algebra
//! @{

//! @brief The positions stored by either pattern of the same shape: the pattern
//! of a sum or a difference.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
inline constexpr auto pattern_union{
    sla::to_pattern<Lhs.rows, Lhs.columns,
                    sla::merge_positions<Lhs, Rhs, true>()>()};

//! @brief The positions stored by both patterns of the same shape.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
inline constexpr auto pattern_intersection{
    sla::to_pattern<Lhs.rows, Lhs.columns,
                    sla::merge_positions<Lhs, Rhs, false>()>()};

//! @brief The transposed pattern.
//!
//! @details The positions are bucketed by column: ascending rows within each
//! column are already what the transposed row-major order requires.
template <sparsity_pattern auto Pattern>
inline constexpr auto pattern_transpose{[] {
  pattern<Pattern.columns, Pattern.rows, Pattern.count> result;
  std::size_t index{0};

  for (std::size_t column{0}; column < Pattern.columns; ++column) {
    for (std::size_t i{0}; i < Pattern.count; ++i) {
      if (Pattern.positions[i].column == column) {
        result.positions[index] = {column, Pattern.positions[i].row};
        ++index;
      }
    }
  }

  return result;
}()};

//! @brief The count of elementary products of a matrix product.
//!
//! @details The count of the scalar multiplications, not of the stored
//! positions of the result: several products may reach one position.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.columns == Rhs.rows)
inline constexpr std::size_t product_term_count{sla::term_count<Lhs, Rhs>};

//! @brief The structural pattern of a matrix product.
//!
//! @details A position is stored whenever a term reaches it, even if the values
//! of the terms may cancel each other out at run time.
template <sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.columns == Rhs.rows)
inline constexpr auto pattern_product{
    sla::to_pattern<Lhs.rows, Rhs.columns,
                    sla::product_positions<Lhs, Rhs>()>()};

//! @brief Whether a pattern stores every position of another, of the same
//! shape.
template <sparsity_pattern auto Inner, sparsity_pattern auto Outer>
inline constexpr bool is_subset_of{[] -> bool {
  if constexpr (Inner.rows != Outer.rows or Inner.columns != Outer.columns) {
    return false;
  } else {
    return pattern_intersection<Inner, Outer>.count == Inner.count;
  }
}()};

//! @}
} // namespace fcarouge::sparse

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_PATTERN_ALGEBRA_HPP

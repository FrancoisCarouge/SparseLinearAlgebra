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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_INTERNAL_COMBINE_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_INTERNAL_COMBINE_HPP

//! @file
//! @brief The element-wise combination of two matrices of the same shape.

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra/pattern_algebra.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <cstddef>

namespace fcarouge::sparse_linear_algebra::internal {
//! @brief The positions of the union of two patterns, grouped by the patterns
//! storing them.
template <auto Lhs, auto Rhs>
inline constexpr const auto &union_terms{
    stored_terms<Lhs, Rhs, sparse::pattern_union<Lhs, Rhs>>};

//! @brief Combine the elements stored by either matrix of the same shape, into
//! the union of their patterns.
//!
//! @param lhs The left matrix.
//! @param rhs The right matrix, of the same shape.
//! @param both The operation on the elements stored by both matrices.
//! @param rhs_only The operation on the elements stored by the rhs matrix only.
//! An element stored by the lhs matrix only is kept as-is.
template <typename Type, auto Lhs, auto Rhs>
[[nodiscard]] constexpr auto combine(const sparse_matrix<Type, Lhs> &lhs,
                                     const sparse_matrix<Type, Rhs> &rhs,
                                     const auto &both, const auto &rhs_only)
    -> sparse_matrix<Type, sparse::pattern_union<Lhs, Rhs>> {
  sparse_matrix<Type, sparse::pattern_union<Lhs, Rhs>> result;

  // The same patterns store the same positions at the same indexes: the
  // combination is a loop over the stored elements, which the compilers
  // vectorize, rather than a walk through a table of indexes.
  if constexpr (Lhs == Rhs) {
    for (std::size_t i{0}; i < Lhs.count; ++i) {
      result.values[i] = both(lhs.values[i], rhs.values[i]);
    }
  } else {
    constexpr const auto &table{union_terms<Lhs, Rhs>};

    for (const term &element : table.both) {
      result.values[element.result] =
          both(lhs.values[element.lhs], rhs.values[element.rhs]);
    }

    for (const term &element : table.lhs_only) {
      result.values[element.result] = lhs.values[element.lhs];
    }

    for (const term &element : table.rhs_only) {
      result.values[element.result] = rhs_only(rhs.values[element.rhs]);
    }
  }

  return result;
}
} // namespace fcarouge::sparse_linear_algebra::internal

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_INTERNAL_COMBINE_HPP

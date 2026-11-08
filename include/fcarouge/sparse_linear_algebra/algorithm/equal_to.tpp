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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_EQUAL_TO_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_EQUAL_TO_TPP

#include "fcarouge/sparse_linear_algebra/pattern_algebra.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <cstddef>

namespace fcarouge::sparse {
//! @brief Element-wise equality of two matrices of the same shape, of any
//! patterns.
//!
//! @details Compares values, not patterns: a structural zero equals a stored
//! zero.
//!
//! @complexity Linear in the count of stored elements of the union of the
//! patterns, each read logarithmic in the count of stored elements.
template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator==(const sparse_matrix<Type, Lhs> &lhs,
                                        const sparse_matrix<Type, Rhs> &rhs)
    -> bool {
  constexpr auto combined{pattern_union<Lhs, Rhs>};

  for (std::size_t i{0}; i < combined.count; ++i) {
    const auto [row, column]{combined.positions[i]};

    if (not(lhs[row, column] == rhs[row, column])) {
      return false;
    }
  }

  return true;
}
} // namespace fcarouge::sparse

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_EQUAL_TO_TPP

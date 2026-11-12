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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_MATRIX_ONE_NORM_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_MATRIX_ONE_NORM_TPP

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <cstddef>

namespace fcarouge {
//! @brief The one norm of a matrix: its largest absolute column sum.
//!
//! @details Of every shape, a vector included: the one norm of a row vector is
//! its largest absolute element, the one norm of a column vector is the sum of
//! its absolute elements.
//!
//! @see std::linalg::matrix_one_norm
//!
//! @complexity Linear in the count of stored elements and of columns.
template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
matrix_one_norm(const sparse_matrix<Type, Pattern> &value) -> Type {
  std::array<Type, Pattern.columns> sums{};
  Type result{};

  // Without a structural zero, the storage is row-major: the column of each
  // stored element is the arithmetic of its index, which the compilers
  // vectorize, rather than a read of its position.
  if constexpr (Pattern.count == Pattern.rows * Pattern.columns) {
    for (std::size_t row{0}; row < Pattern.rows; ++row) {
      for (std::size_t column{0}; column < Pattern.columns; ++column) {
        sums[column] += sparse_linear_algebra::internal::absolute(
            value.values[(row * Pattern.columns) + column]);
      }
    }
  } else {
    for (std::size_t i{0}; i < Pattern.count; ++i) {
      sums[Pattern.positions[i].column] +=
          sparse_linear_algebra::internal::absolute(value.values[i]);
    }
  }

  for (const Type &sum : sums) {
    result = result < sum ? sum : result;
  }

  return result;
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_MATRIX_ONE_NORM_TPP

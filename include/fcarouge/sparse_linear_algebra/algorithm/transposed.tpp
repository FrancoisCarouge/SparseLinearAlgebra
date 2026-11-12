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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_TRANSPOSED_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_TRANSPOSED_TPP

#include "fcarouge/sparse_linear_algebra/pattern_algebra.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <cstddef>

namespace fcarouge {
namespace sparse_linear_algebra::internal {
//! @brief For each stored element of the transposed pattern, the index of its
//! mirrored stored element.
template <auto Pattern>
inline constexpr auto transpose_indexes{[] {
  constexpr auto mirrored{sparse::pattern_transpose<Pattern>};
  std::array<std::size_t, Pattern.count> result{};

  for (std::size_t i{0}; i < Pattern.count; ++i) {
    result[i] = Pattern.index_of(mirrored.positions[i].column,
                                 mirrored.positions[i].row);
  }

  return result;
}()};
} // namespace sparse_linear_algebra::internal

//! @brief The transposed matrix, of the transposed pattern.
//!
//! @details Returns a matrix, by value, rather than a view.
//!
//! @see std::linalg::transposed
//!
//! @complexity Linear in the count of stored elements.
template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
transposed(const sparse_matrix<Type, Pattern> &value)
    -> sparse_matrix<Type, sparse::pattern_transpose<Pattern>> {
  sparse_matrix<Type, sparse::pattern_transpose<Pattern>> result;

  // Without a structural zero, the storage is row-major: the transposition is
  // the arithmetic of the indexes, in the order of the stored elements of the
  // result, rather than a walk through a table of indexes.
  if constexpr (Pattern.count == Pattern.rows * Pattern.columns) {
    for (std::size_t column{0}; column < Pattern.columns; ++column) {
      for (std::size_t row{0}; row < Pattern.rows; ++row) {
        result.values[(column * Pattern.rows) + row] =
            value.values[(row * Pattern.columns) + column];
      }
    }
  } else {
    constexpr const auto &indexes{
        sparse_linear_algebra::internal::transpose_indexes<Pattern>};

    for (std::size_t i{0}; i < Pattern.count; ++i) {
      result.values[i] = value.values[indexes[i]];
    }
  }

  return result;
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_TRANSPOSED_TPP

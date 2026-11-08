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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DOT_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DOT_TPP

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

namespace fcarouge {
//! @brief The dot product of two vectors of the same size, of any orientation.
//!
//! @details The sum of the products of the elements stored by both vectors.
//!
//! @see std::linalg::dot
//!
//! @complexity Linear in the count of elements stored by both vectors.
template <typename Type, sparse::sparsity_pattern auto Lhs,
          sparse::sparsity_pattern auto Rhs>
  requires((Lhs.rows == 1 or Lhs.columns == 1) and
           (Rhs.rows == 1 or Rhs.columns == 1) and
           Lhs.rows * Lhs.columns == Rhs.rows * Rhs.columns)
[[nodiscard]] constexpr auto dot(const sparse_matrix<Type, Lhs> &lhs,
                                 const sparse_matrix<Type, Rhs> &rhs) -> Type {
  Type result{};

  for (const sparse_linear_algebra::internal::term &element :
       sparse_linear_algebra::internal::dot_terms<Lhs, Rhs>) {
    result += lhs.values[element.lhs] * rhs.values[element.rhs];
  }

  return result;
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DOT_TPP

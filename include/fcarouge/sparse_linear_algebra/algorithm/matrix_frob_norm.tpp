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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_MATRIX_FROB_NORM_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_MATRIX_FROB_NORM_TPP

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <concepts>

namespace fcarouge {
//! @brief The Frobenius norm of a matrix.
//!
//! @details The square root of the sum of the squares of the stored elements.
//! Only for a floating-point element type: an integer has no exact square root.
//!
//! @note At compile time, and on a freestanding implementation, the square
//! root is the library's own, whose result is not guaranteed to be identical
//! to the last digit to that of `std::sqrt`.
//!
//! @see std::linalg::matrix_frob_norm
//!
//! @complexity Linear in the count of stored elements.
template <std::floating_point Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
matrix_frob_norm(const sparse_matrix<Type, Pattern> &value) -> Type {
  Type result{};

  for (const Type &element : value.values) {
    result += element * element;
  }

  return sparse_linear_algebra::internal::square_root(result);
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_MATRIX_FROB_NORM_TPP

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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DIVIDE_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DIVIDE_TPP

#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <type_traits>

namespace fcarouge {
//! @brief Scaling by the inverse of a scalar divisor, of the same pattern.
//!
//! @note A zero divisor leaves every structural zero an exact zero: see the
//! note on non-finite values.
//!
//! @todo The division by a matrix needs a solver, a factorization.
//!
//! @complexity Linear in the count of stored elements.
template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator/(const sparse_matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> sparse_matrix<Type, Pattern> {
  sparse_matrix<Type, Pattern> result{lhs};

  for (Type &element : result.values) {
    element = element / rhs;
  }

  return result;
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DIVIDE_TPP

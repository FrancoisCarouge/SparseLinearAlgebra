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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_ADD_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_ADD_TPP

#include "fcarouge/sparse_linear_algebra/internal/combine.hpp"

namespace fcarouge::sparse_linear_algebra::internal {
//! @brief The sum of the elements stored by both terms.
//!
//! @note The element operations are variable templates of the element type
//! only, not lambdas of the operator: those would be instantiated for every
//! pair of patterns, and never called for a pair without a common position, or
//! without a position stored by the right term only.
template <typename Type>
inline constexpr auto sum{
    [](const Type &left, const Type &right) -> Type { return left + right; }};

//! @brief The element stored by the right term only, copied.
template <typename Type>
inline constexpr auto identity{[](const Type &right) -> Type { return right; }};
} // namespace fcarouge::sparse_linear_algebra::internal

namespace fcarouge {
//! @brief Addition of two matrices of the same shape, of any patterns.
//!
//! @details The result stores the union of the terms' patterns. An element
//! stored by both terms costs one addition, an element stored by one term only
//! is copied.
//!
//! @complexity Linear in the count of stored elements of the result.
template <typename Type, sparse::sparsity_pattern auto Lhs,
          sparse::sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator+(const sparse_matrix<Type, Lhs> &lhs,
                                       const sparse_matrix<Type, Rhs> &rhs)
    -> sparse_matrix<Type, sparse::pattern_union<Lhs, Rhs>> {
  return sla::combine(lhs, rhs, sla::sum<Type>, sla::identity<Type>);
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_ADD_TPP

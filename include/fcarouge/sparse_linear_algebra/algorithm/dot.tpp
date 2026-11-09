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

#include <cstddef>
#include <utility>

namespace fcarouge {
namespace sparse_linear_algebra::internal {
//! @brief The count of elementary products of a dot product.
//!
//! @details Named, for CppCheck, which misparses a member access in a template
//! argument list.
template <auto Lhs, auto Rhs>
inline constexpr std::size_t dot_term_count{dot_terms<Lhs, Rhs>.size()};

//! @brief The dot product: the sum of its elementary products.
//!
//! @details Each elementary product is a separate expression of compile-time
//! indexes: the sum is one expression, kept in registers, rather than a loop
//! over a table of indexes. The sum starts from a zero element and adds the
//! products in the order of the stored elements, as the accumulation into a
//! zero result does.
template <auto Lhs, auto Rhs, typename Type, std::size_t... Terms>
[[nodiscard]] constexpr auto
dot_product(const sparse_matrix<Type, Lhs> &lhs,
            const sparse_matrix<Type, Rhs> &rhs,
            std::index_sequence<Terms...> /*terms*/) -> Type {
  constexpr const auto &terms{dot_terms<Lhs, Rhs>};

  return (Type{} + ... +
          (lhs.values[terms[Terms].lhs] * rhs.values[terms[Terms].rhs]));
}
} // namespace sparse_linear_algebra::internal

//! @brief The dot product of two vectors of the same size, of any orientation.
//!
//! @details The sum of the products of the elements stored by both vectors.
//!
//! @see std::linalg::dot
//!
//! @complexity Linear in the count of elements stored by both vectors.
template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires((Lhs.rows == 1 or Lhs.columns == 1) and
           (Rhs.rows == 1 or Rhs.columns == 1) and
           Lhs.rows * Lhs.columns == Rhs.rows * Rhs.columns)
[[nodiscard]] constexpr auto dot(const sparse_matrix<Type, Lhs> &lhs,
                                 const sparse_matrix<Type, Rhs> &rhs) -> Type {
  return sparse_linear_algebra::internal::dot_product(
      lhs, rhs,
      std::make_index_sequence<
          sparse_linear_algebra::internal::dot_term_count<Lhs, Rhs>>{});
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DOT_TPP

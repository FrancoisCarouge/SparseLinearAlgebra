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

//! @brief The count of fold expressions expanding the elementary products of
//! a dot product.
template <auto Lhs, auto Rhs>
inline constexpr std::size_t dot_fold_count{
    fold_count<dot_term_count<Lhs, Rhs>>};

//! @brief The count of elementary products of a fold expression of a dot
//! product.
template <auto Lhs, auto Rhs, std::size_t Fold>
inline constexpr std::size_t dot_fold_term_count{
    fold_term_count<dot_term_count<Lhs, Rhs>, Fold>};

//! @brief Add a run of elementary products of a dot product to a partial sum.
//!
//! @details Each elementary product is a separate expression of compile-time
//! indexes, rather than a loop over a table of indexes loaded from memory. The
//! products add to the partial sum in the order of the stored elements.
template <auto Lhs, auto Rhs, std::size_t Fold, typename Type,
          std::size_t... Terms>
[[nodiscard]] constexpr auto dot_fold(const Type &sum,
                                      const sparse_matrix<Type, Lhs> &lhs,
                                      const sparse_matrix<Type, Rhs> &rhs,
                                      std::index_sequence<Terms...> /*terms*/)
    -> Type {
  constexpr const auto &terms{dot_terms<Lhs, Rhs>};
  constexpr std::size_t first{Fold * fold_size};

  return (sum + ... +
          (lhs.values[terms[first + Terms].lhs] *
           rhs.values[terms[first + Terms].rhs]));
}

//! @brief A dot product, the sum of all its elementary products from zero.
template <auto Lhs, auto Rhs, typename Type, std::size_t... Folds>
[[nodiscard]] constexpr auto
dot_product(const sparse_matrix<Type, Lhs> &lhs,
            const sparse_matrix<Type, Rhs> &rhs,
            std::index_sequence<Folds...> /*folds*/) -> Type {
  Type result{};

  ((result = dot_fold<Lhs, Rhs, Folds>(
        result, lhs, rhs,
        std::make_index_sequence<dot_fold_term_count<Lhs, Rhs, Folds>>{})),
   ...);

  return result;
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
          sparse_linear_algebra::internal::dot_fold_count<Lhs, Rhs>>{});
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_DOT_TPP

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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_PRODUCT_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_PRODUCT_TPP

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra/pattern_algebra.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <cstddef>
#include <type_traits>

namespace fcarouge {
namespace sparse_linear_algebra::internal {
//! @brief The elementary products of a matrix product.
template <auto Lhs, auto Rhs>
inline constexpr auto product_terms{
    [] -> std::array<term, sparse::product_term_count<Lhs, Rhs>> {
      constexpr auto result_pattern{sparse::pattern_product<Lhs, Rhs>};
      std::array<term, sparse::product_term_count<Lhs, Rhs>> result{};
      std::size_t index{0};

      for (std::size_t i{0}; i < Lhs.count; ++i) {
        const auto [row, inner]{Lhs.positions[i]};
        const std::size_t inner_end{Rhs.lower_bound({inner + 1, 0})};

        for (std::size_t j{Rhs.lower_bound({inner, 0})}; j < inner_end; ++j) {
          result[index] = {
              .result = result_pattern.index_of(row, Rhs.positions[j].column),
              .lhs = i,
              .rhs = j};
          ++index;
        }
      }

      return result;
    }()};
} // namespace sparse_linear_algebra::internal

//! @brief Matrix multiplication, of matrices and vectors of compatible shapes,
//! of any patterns.
//!
//! @details The result stores the structural product of the factors' patterns:
//! a position reached by at least one elementary product. Only the elementary
//! products the patterns admit are computed.
//!
//! @complexity Linear in the count of stored elements of the result and of
//! elementary products, `product_term_count`.
template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.columns == Rhs.rows)
[[nodiscard]] constexpr auto operator*(const sparse_matrix<Type, Lhs> &lhs,
                                       const sparse_matrix<Type, Rhs> &rhs)
    -> sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> {
  sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> result;

  for (const sparse_linear_algebra::internal::term &element :
       sparse_linear_algebra::internal::product_terms<Lhs, Rhs>) {
    result.values[element.result] +=
        lhs.values[element.lhs] * rhs.values[element.rhs];
  }

  return result;
}

//! @brief Scaling by a scalar factor, on the right, of the same pattern.
//!
//! @note A structural zero stays an exact zero, whatever the factor: see the
//! note on non-finite values.
//!
//! @complexity Linear in the count of stored elements.
template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator*(const sparse_matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> sparse_matrix<Type, Pattern> {
  sparse_matrix<Type, Pattern> result{lhs};

  for (Type &element : result.values) {
    element = element * rhs;
  }

  return result;
}

//! @brief Scaling by a scalar factor, on the left, of the same pattern.
//!
//! @note A structural zero stays an exact zero, whatever the factor: see the
//! note on non-finite values.
//!
//! @complexity Linear in the count of stored elements.
template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator*(const std::type_identity_t<Type> &lhs,
                                       const sparse_matrix<Type, Pattern> &rhs)
    -> sparse_matrix<Type, Pattern> {
  sparse_matrix<Type, Pattern> result{rhs};

  for (Type &element : result.values) {
    element = lhs * element;
  }

  return result;
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_ALGORITHM_PRODUCT_TPP

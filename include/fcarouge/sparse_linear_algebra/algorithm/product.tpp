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
#include <utility>

namespace fcarouge {
namespace sparse_linear_algebra::internal {
//! @brief The elementary products of a matrix product, grouped by result
//! element.
//!
//! @details The terms of the result element `i` are
//! `terms[offsets[i]]` to `terms[offsets[i + 1]]`, in the order of the stored
//! elements of the left factor.
template <std::size_t TermCount, std::size_t ResultCount> struct product_table {
  std::array<term, TermCount> terms{};
  std::array<std::size_t, ResultCount + 1> offsets{};
};

//! @brief The count of stored elements of a matrix product.
//!
//! @details Named, for CppCheck, which misparses a member access in a template
//! argument list.
template <auto Lhs, auto Rhs>
inline constexpr std::size_t product_element_count{
    sparse::pattern_product<Lhs, Rhs>.count};

//! @brief The table of the elementary products of a matrix product.
template <auto Lhs, auto Rhs>
using product_table_of = product_table<sparse::product_term_count<Lhs, Rhs>,
                                       product_element_count<Lhs, Rhs>>;

//! @brief The elementary products of a matrix product.
template <auto Lhs, auto Rhs>
inline constexpr auto product_terms{[] -> product_table_of<Lhs, Rhs> {
  constexpr auto result_pattern{sparse::pattern_product<Lhs, Rhs>};
  std::array<term, sparse::product_term_count<Lhs, Rhs>> listed{};
  product_table_of<Lhs, Rhs> result{};
  std::size_t index{0};

  for (std::size_t i{0}; i < Lhs.count; ++i) {
    const auto [row, inner]{Lhs.positions[i]};
    const std::size_t inner_end{Rhs.lower_bound({inner + 1, 0})};

    for (std::size_t j{Rhs.lower_bound({inner, 0})}; j < inner_end; ++j) {
      listed[index] = {
          .result = result_pattern.index_of(row, Rhs.positions[j].column),
          .lhs = i,
          .rhs = j};
      ++result.offsets[listed[index].result + 1];
      ++index;
    }
  }

  for (std::size_t i{0}; i < result_pattern.count; ++i) {
    result.offsets[i + 1] += result.offsets[i];
  }

  // A stable counting sort by result element.
  std::array<std::size_t, result_pattern.count + 1> next{result.offsets};

  for (const term &element : listed) {
    result.terms[next[element.result]] = element;
    ++next[element.result];
  }

  return result;
}()};

//! @brief The count of elementary products of a stored element of a matrix
//! product.
//!
//! @details Named, for CppCheck, which misparses a member access in a template
//! argument list.
template <auto Lhs, auto Rhs, std::size_t Index>
inline constexpr std::size_t product_element_term_count{
    product_terms<Lhs, Rhs>.offsets[Index + 1] -
    product_terms<Lhs, Rhs>.offsets[Index]};

//! @brief A stored element of a matrix product: the sum of its elementary
//! products.
//!
//! @details Each elementary product is a separate expression of compile-time
//! indexes: the sum is one expression, kept in registers, rather than a loop
//! over a table of indexes accumulating in memory. The sum starts from a zero
//! element and adds the products in the order of the stored elements of the
//! left factor, as the accumulation into a zero-initialized result does.
template <auto Lhs, auto Rhs, std::size_t Index, typename Type,
          std::size_t... Terms>
[[nodiscard]] constexpr auto
product_element(const sparse_matrix<Type, Lhs> &lhs,
                const sparse_matrix<Type, Rhs> &rhs,
                std::index_sequence<Terms...> /*terms*/) -> Type {
  constexpr std::size_t first{product_terms<Lhs, Rhs>.offsets[Index]};
  constexpr const auto &terms{product_terms<Lhs, Rhs>.terms};

  return (Type{} + ... +
          (lhs.values[terms[first + Terms].lhs] *
           rhs.values[terms[first + Terms].rhs]));
}

//! @brief Every stored element of a matrix product.
template <auto Lhs, auto Rhs, typename Type, std::size_t... Indexes>
[[nodiscard]] constexpr auto
product(const sparse_matrix<Type, Lhs> &lhs,
        const sparse_matrix<Type, Rhs> &rhs,
        std::index_sequence<Indexes...> /*indexes*/)
    -> sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> {
  sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> result;

  ((result.values[Indexes] = product_element<Lhs, Rhs, Indexes>(
        lhs, rhs,
        std::make_index_sequence<
            product_element_term_count<Lhs, Rhs, Indexes>>{})),
   ...);

  return result;
}
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
  return sparse_linear_algebra::internal::product(
      lhs, rhs,
      std::make_index_sequence<
          sparse_linear_algebra::internal::product_element_count<Lhs, Rhs>>{});
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

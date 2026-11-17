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
//! element, in the order of the stored elements of the left factor within each
//! group.
template <auto Lhs, auto Rhs>
inline constexpr auto product_terms{
    [] -> std::array<term, sparse::product_term_count<Lhs, Rhs>> {
      constexpr auto result_pattern{sparse::pattern_product<Lhs, Rhs>};
      std::array<term, sparse::product_term_count<Lhs, Rhs>> listed{};
      std::array<std::size_t, result_pattern.count + 1> next{};
      std::size_t index{0};

      for (std::size_t i{0}; i < Lhs.count; ++i) {
        const auto [row, inner]{Lhs.positions[i]};
        const std::size_t inner_end{row_bounds<Rhs>[inner + 1]};

        for (std::size_t j{row_bounds<Rhs>[inner]}; j < inner_end; ++j) {
          listed[index] = {
              .result = result_pattern.index_of(row, Rhs.positions[j].column),
              .lhs = i,
              .rhs = j};
          ++next[listed[index].result + 1];
          ++index;
        }
      }

      // A stable counting sort by result element: the first position of each
      // group, then the terms in their listed order.
      for (std::size_t i{0}; i < result_pattern.count; ++i) {
        next[i + 1] += next[i];
      }

      std::array<term, sparse::product_term_count<Lhs, Rhs>> result{};

      for (const term &element : listed) {
        result[next[element.result]] = element;
        ++next[element.result];
      }

      return result;
    }()};

//! @brief The count of fold expressions expanding the elementary products of
//! a matrix product.
template <auto Lhs, auto Rhs>
inline constexpr std::size_t product_fold_count{
    fold_count<sparse::product_term_count<Lhs, Rhs>>};

//! @brief The count of elementary products of a fold expression of a matrix
//! product.
template <auto Lhs, auto Rhs, std::size_t Fold>
inline constexpr std::size_t product_fold_term_count{
    fold_term_count<sparse::product_term_count<Lhs, Rhs>, Fold>};

//! @brief Accumulate a run of elementary products into their result elements.
//!
//! @details Each elementary product is a separate expression of compile-time
//! indexes, rather than a loop over a table of indexes loaded from memory: the
//! compilers keep the result elements in registers. The products accumulate
//! into the zero-initialized result in the order of the table.
template <auto Lhs, auto Rhs, std::size_t Fold, typename Type,
          std::size_t... Terms>
constexpr void
product_fold(sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> &result,
             const sparse_matrix<Type, Lhs> &lhs,
             const sparse_matrix<Type, Rhs> &rhs,
             std::index_sequence<Terms...> terms) {
  static_cast<void>(terms);

  constexpr const auto &table{product_terms<Lhs, Rhs>};
  constexpr std::size_t first{Fold * fold_size};

  ((result.values[table[first + Terms].result] +=
    lhs.values[table[first + Terms].lhs] *
    rhs.values[table[first + Terms].rhs]),
   ...);
}

//! @brief The most elementary products a matrix product expands at compile
//! time, beyond which it accumulates them in a loop.
//!
//! @details The expanded products keep the small products in registers, where
//! a loop over the stored elements is several times slower. Each expanded
//! product is an expression of its own, and the compilers' time and memory
//! grow with them: beyond a few thousand, their cost bounds the size of the
//! products a program can compile, while the loop's run time catches up with
//! the expanded products'.
inline constexpr std::size_t unrolled_term_limit{4096};

//! @brief A matrix product, of all its elementary products, accumulated a row
//! of the result at a time.
//!
//! @details For each row of the left factor, each stored element scales the row
//! of the right factor it selects into an accumulator of a row of the result,
//! whose stored elements are then gathered, and the accumulator cleared. The
//! products accumulate into each result element in the order of the stored
//! elements of the left factor, the order of the expanded products: both
//! compute the same values. The compile-time tables are the row bounds of the
//! patterns, linear in their stored elements, rather than the elementary
//! products.
template <auto Lhs, auto Rhs, typename Type>
[[nodiscard]] constexpr auto row_product(const sparse_matrix<Type, Lhs> &lhs,
                                         const sparse_matrix<Type, Rhs> &rhs)
    -> sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> {
  constexpr auto result_pattern{sparse::pattern_product<Lhs, Rhs>};
  constexpr const auto &lhs_bounds{row_bounds<Lhs>};
  constexpr const auto &rhs_bounds{row_bounds<Rhs>};
  constexpr const auto &result_bounds{row_bounds<result_pattern>};

  sparse_matrix<Type, result_pattern> result;
  std::array<Type, column_count<Rhs>> accumulator{};

  for (std::size_t row{0}; row < Lhs.rows; ++row) {
    for (std::size_t i{lhs_bounds[row]}; i < lhs_bounds[row + 1]; ++i) {
      const std::size_t inner{Lhs.positions[i].column};

      for (std::size_t j{rhs_bounds[inner]}; j < rhs_bounds[inner + 1]; ++j) {
        accumulator[Rhs.positions[j].column] += lhs.values[i] * rhs.values[j];
      }
    }

    for (std::size_t k{result_bounds[row]}; k < result_bounds[row + 1]; ++k) {
      Type &element{accumulator[result_pattern.positions[k].column]};

      result.values[k] = element;
      element = Type{};
    }
  }

  return result;
}

//! @brief A matrix product, of all its elementary products.
template <auto Lhs, auto Rhs, typename Type, std::size_t... Folds>
[[nodiscard]] constexpr auto product(const sparse_matrix<Type, Lhs> &lhs,
                                     const sparse_matrix<Type, Rhs> &rhs,
                                     std::index_sequence<Folds...> folds)
    -> sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> {
  static_cast<void>(folds);

  sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>> result;

  (product_fold<Lhs, Rhs, Folds>(
       result, lhs, rhs,
       std::make_index_sequence<product_fold_term_count<Lhs, Rhs, Folds>>{}),
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
  if constexpr (sparse::product_term_count<Lhs, Rhs> <=
                sparse_linear_algebra::internal::unrolled_term_limit) {
    return sparse_linear_algebra::internal::product(
        lhs, rhs,
        std::make_index_sequence<
            sparse_linear_algebra::internal::product_fold_count<Lhs, Rhs>>{});
  } else {
    return sparse_linear_algebra::internal::row_product(lhs, rhs);
  }
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

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

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra/pattern_algebra.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace fcarouge {
namespace sparse_linear_algebra::internal {
//! @brief A copy of a stored element into the elimination buffer of a block,
//! or out of it.
struct transfer {
  std::size_t buffer{0};
  std::size_t element{0};
};

//! @brief The subtraction from an element of the elimination buffer of a block
//! of the elementary product of a stored element of the denominator and of an
//! element of the quotient, of an unknown of a block solved before.
struct update {
  std::size_t buffer{0};
  std::size_t denominator{0};
  std::size_t quotient{0};
};

//! @brief A diagonal block of a division: its count of unknowns, the width of
//! the rows of its elimination buffer, and the end of its transfers and of its
//! updates in the tables.
struct division_block {
  std::size_t size{0};
  std::size_t width{0};
  std::size_t denominator_end{0};
  std::size_t numerator_end{0};
  std::size_t update_end{0};
  std::size_t store_end{0};
};

//! @brief The counts of the entries of the tables of a division, and the size
//! of the largest elimination buffer of its blocks.
struct division_counts {
  std::size_t blocks{0};
  std::size_t denominator_loads{0};
  std::size_t numerator_loads{0};
  std::size_t updates{0};
  std::size_t stores{0};
  std::size_t buffer_size{0};

  constexpr void load_denominator(std::size_t /*buffer*/,
                                  std::size_t /*element*/) {
    ++denominator_loads;
  }

  constexpr void load_numerator(std::size_t /*buffer*/,
                                std::size_t /*element*/) {
    ++numerator_loads;
  }

  constexpr void subtract(std::size_t /*buffer*/, std::size_t /*denominator*/,
                          std::size_t /*quotient*/) {
    ++updates;
  }

  constexpr void store(std::size_t /*buffer*/, std::size_t /*element*/) {
    ++stores;
  }

  constexpr void close_block(std::size_t size, std::size_t width) {
    ++blocks;
    buffer_size = size * width > buffer_size ? size * width : buffer_size;
  }
};

//! @brief The tables of a division, block by block: the stored elements of the
//! denominator and of the numerator loaded into the elimination buffer, the
//! updates by the unknowns already solved, and the stores of the solution into
//! the quotient.
template <std::size_t Blocks, std::size_t DenominatorLoads,
          std::size_t NumeratorLoads, std::size_t Updates, std::size_t Stores>
struct division_tables {
  std::array<division_block, Blocks> blocks{};
  std::array<transfer, DenominatorLoads> denominator_loads{};
  std::array<transfer, NumeratorLoads> numerator_loads{};
  std::array<update, Updates> updates{};
  std::array<transfer, Stores> stores{};
  division_counts filled{};

  constexpr void load_denominator(std::size_t buffer, std::size_t element) {
    denominator_loads[filled.denominator_loads] = {.buffer = buffer,
                                                   .element = element};
    filled.load_denominator(buffer, element);
  }

  constexpr void load_numerator(std::size_t buffer, std::size_t element) {
    numerator_loads[filled.numerator_loads] = {.buffer = buffer,
                                               .element = element};
    filled.load_numerator(buffer, element);
  }

  constexpr void subtract(std::size_t buffer, std::size_t denominator,
                          std::size_t quotient) {
    updates[filled.updates] = {
        .buffer = buffer, .denominator = denominator, .quotient = quotient};
    filled.subtract(buffer, denominator, quotient);
  }

  constexpr void store(std::size_t buffer, std::size_t element) {
    stores[filled.stores] = {.buffer = buffer, .element = element};
    filled.store(buffer, element);
  }

  constexpr void close_block(std::size_t size, std::size_t width) {
    filled.close_block(size, width);
    blocks[filled.blocks - 1] = {.size = size,
                                 .width = width,
                                 .denominator_end = filled.denominator_loads,
                                 .numerator_end = filled.numerator_loads,
                                 .update_end = filled.updates,
                                 .store_end = filled.stores};
  }
};

//! @brief Plan the division `x * rhs = lhs`, block by block of the block
//! triangular form of the denominator, into counts or tables.
//!
//! @details The elimination buffer of a block holds a row per unknown of the
//! block, for its matched equation: the stored elements of the equation at the
//! unknowns of the block, then the right-hand side at the rows of the quotient
//! the block reaches, the stored elements of the numerator less the updates by
//! the unknowns of the blocks solved before. After the elimination, a row holds
//! the solution of its unknown.
template <auto Lhs, auto Rhs, typename Tables>
consteval void plan_division(Tables &tables) {
  constexpr auto form{block_triangular_form<Rhs>};
  constexpr auto reached{quotient_reach<Lhs, Rhs>};
  constexpr auto transposed{
      to_pattern<Rhs.columns, Rhs.rows, transpose_positions<Rhs>()>()};
  constexpr auto lhs_transposed{
      to_pattern<Lhs.columns, Lhs.rows, transpose_positions<Lhs>()>()};
  constexpr auto quotient{
      to_pattern<Lhs.rows, Rhs.columns, quotient_positions<Lhs, Rhs>()>()};

  for (std::size_t block{0}; block < form.blocks; ++block) {
    const std::size_t first{form.begin[block]};
    const std::size_t size{form.begin[block + 1] - first};
    std::array<std::size_t, Lhs.rows> rows{};
    std::array<std::size_t, Lhs.rows> column_of{};
    std::size_t reached_count{0};

    for (std::size_t row{0}; row < Lhs.rows; ++row) {
      if (reached[form.order[first]][row]) {
        rows[reached_count] = row;
        column_of[row] = size + reached_count;
        ++reached_count;
      }
    }

    const std::size_t width{size + reached_count};

    for (std::size_t k{0}; k < size; ++k) {
      const std::size_t unknown{form.order[first + k]};
      const std::size_t equation{form.equation[unknown]};
      const std::size_t end{transposed.lower_bound({equation + 1, 0})};
      const std::size_t lhs_end{lhs_transposed.lower_bound({equation + 1, 0})};

      for (std::size_t i{transposed.lower_bound({equation, 0})}; i < end; ++i) {
        const std::size_t dependency{transposed.positions[i].column};
        const std::size_t element{Rhs.index_of(dependency, equation)};

        if (form.block[dependency] == block) {
          tables.load_denominator((k * width) + form.rank[dependency] - first,
                                  element);
        } else {
          for (std::size_t j{0}; j < reached_count; ++j) {
            if (reached[dependency][rows[j]]) {
              tables.subtract((k * width) + size + j, element,
                              quotient.index_of(rows[j], dependency));
            }
          }
        }
      }

      for (std::size_t i{lhs_transposed.lower_bound({equation, 0})};
           i < lhs_end; ++i) {
        const std::size_t row{lhs_transposed.positions[i].column};

        tables.load_numerator((k * width) + column_of[row],
                              Lhs.index_of(row, equation));
      }

      for (std::size_t j{0}; j < reached_count; ++j) {
        tables.store((k * width) + size + j,
                     quotient.index_of(rows[j], unknown));
      }
    }

    tables.close_block(size, width);
  }
}

//! @brief The counts of the entries of the tables of a division.
template <auto Lhs, auto Rhs>
inline constexpr division_counts division_count{[] consteval {
  division_counts result{};

  plan_division<Lhs, Rhs>(result);

  return result;
}()};

//! @brief The tables of a division.
template <auto Lhs, auto Rhs>
inline constexpr auto division_plan{[] consteval {
  constexpr division_counts counts{division_count<Lhs, Rhs>};
  division_tables<counts.blocks, counts.denominator_loads,
                  counts.numerator_loads, counts.updates, counts.stores>
      result{};

  plan_division<Lhs, Rhs>(result);

  return result;
}()};
} // namespace sparse_linear_algebra::internal

//! @brief Matrix division, the quotient `x` of `x * rhs = lhs`, by a square
//! denominator of full structural rank, of any patterns.
//!
//! @details The system is solved block by block of the block triangular form of
//! the denominator: an irreducible diagonal block by a Gaussian elimination
//! with partial pivoting, after the subtraction of the terms of the blocks
//! solved before. A block of a single unknown, every block of a diagonal or of
//! a triangular denominator, is a division: the solution is a substitution,
//! without pivoting. A denominator of a structural rank below its size is
//! singular whatever its values, and is rejected. The result stores the
//! structural quotient of the patterns, `pattern_quotient`.
//!
//! @warning The denominator must be invertible: a singular denominator divides
//! by a zero pivot.
//!
//! @note Within an irreducible block of more than one unknown, the elimination
//! is dense: a structural zero of the block takes part in the elimination as a
//! zero, and a non-finite stored element may propagate through it. Elsewhere, a
//! structural zero is never computed: see the note on non-finite values.
//!
//! @complexity Linear in the count of stored elements of the operands and of
//! the result, and in the count of terms coupling the blocks, plus, for each
//! irreducible block, cubic in its count of unknowns and linear in the count
//! of rows of the quotient it reaches.
template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Rhs.rows == Rhs.columns and Lhs.columns == Rhs.rows and
           sparse::structural_rank<Rhs> == Rhs.rows)
[[nodiscard]] constexpr auto operator/(const sparse_matrix<Type, Lhs> &lhs,
                                       const sparse_matrix<Type, Rhs> &rhs)
    -> sparse_matrix<Type, sparse::pattern_quotient<Lhs, Rhs>> {
  constexpr auto &plan{
      sparse_linear_algebra::internal::division_plan<Lhs, Rhs>};
  sparse_matrix<Type, sparse::pattern_quotient<Lhs, Rhs>> result;
  std::array<Type, plan.filled.buffer_size> buffer{};
  std::size_t denominator_load{0};
  std::size_t numerator_load{0};
  std::size_t update{0};
  std::size_t store{0};

  for (const sparse_linear_algebra::internal::division_block &block :
       plan.blocks) {
    const std::size_t size{block.size};
    const std::size_t width{block.width};

    for (std::size_t i{0}; i < size * width; ++i) {
      buffer[i] = Type{};
    }

    for (; denominator_load < block.denominator_end; ++denominator_load) {
      const auto [index, element]{plan.denominator_loads[denominator_load]};

      buffer[index] = rhs.values[element];
    }

    for (; numerator_load < block.numerator_end; ++numerator_load) {
      const auto [index, element]{plan.numerator_loads[numerator_load]};

      buffer[index] = lhs.values[element];
    }

    for (; update < block.update_end; ++update) {
      const auto [index, denominator, quotient]{plan.updates[update]};

      buffer[index] =
          buffer[index] - (rhs.values[denominator] * result.values[quotient]);
    }

    for (std::size_t pivot{0}; pivot < size; ++pivot) {
      std::size_t largest{pivot};

      for (std::size_t row{pivot + 1}; row < size; ++row) {
        if (sparse_linear_algebra::internal::absolute(
                buffer[(largest * width) + pivot]) <
            sparse_linear_algebra::internal::absolute(
                buffer[(row * width) + pivot])) {
          largest = row;
        }
      }

      if (largest != pivot) {
        for (std::size_t column{pivot}; column < width; ++column) {
          std::swap(buffer[(pivot * width) + column],
                    buffer[(largest * width) + column]);
        }
      }

      for (std::size_t row{pivot + 1}; row < size; ++row) {
        const Type factor{buffer[(row * width) + pivot] /
                          buffer[(pivot * width) + pivot]};

        for (std::size_t column{pivot + 1}; column < width; ++column) {
          buffer[(row * width) + column] =
              buffer[(row * width) + column] -
              (factor * buffer[(pivot * width) + column]);
        }
      }
    }

    for (std::size_t i{0}; i < size; ++i) {
      const std::size_t row{size - 1 - i};

      for (std::size_t column{size}; column < width; ++column) {
        Type value{buffer[(row * width) + column]};

        for (std::size_t inner{row + 1}; inner < size; ++inner) {
          value = value - (buffer[(row * width) + inner] *
                           buffer[(inner * width) + column]);
        }

        buffer[(row * width) + column] = value / buffer[(row * width) + row];
      }
    }

    for (; store < block.store_end; ++store) {
      const auto [index, element]{plan.stores[store]};

      result.values[element] = buffer[index];
    }
  }

  return result;
}

//! @brief Scaling by the inverse of a scalar divisor, of the same pattern.
//!
//! @note A zero divisor leaves every structural zero an exact zero: see the
//! note on non-finite values.
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

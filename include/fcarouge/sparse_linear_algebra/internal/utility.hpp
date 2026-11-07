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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP

//! @file
//! @brief The internal implementation details: the pattern computations the
//! operations derive their results and their arithmetic terms with.

#include "fcarouge/sparse_linear_algebra/pattern.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>

#if __STDC_HOSTED__

#include <cmath>

#endif

namespace fcarouge::sparse_linear_algebra::internal {
//! @brief Whether a type is a sparse matrix.
template <typename Type> inline constexpr bool is_matrix{false};

template <typename Type, sparse::sparsity_pattern auto Pattern>
inline constexpr bool is_matrix<sparse::matrix<Type, Pattern>>{true};

//! @brief The element type of a sparse matrix.
template <typename Matrix> struct element {};

template <typename Type, sparse::sparsity_pattern auto Pattern>
struct element<sparse::matrix<Type, Pattern>> {
  using type = Type;
};

//! @brief The element type of a sparse matrix.
template <typename Matrix> using element_t = element<Matrix>::type;

//! @brief Whether a pattern stores every position: a pattern without a
//! structural zero.
template <auto Pattern>
inline constexpr bool is_dense{Pattern.count == Pattern.rows * Pattern.columns};

//! @brief A growing list of positions, of a known maximum size, from which the
//! pattern of an operation's result is sized.
template <std::size_t Capacity> struct position_list {
  std::array<sparse::position, (Capacity > 0 ? Capacity : 1)> items{};
  std::size_t size{0};

  constexpr void push(sparse::position value) {
    items[size] = value;
    ++size;
  }
};

//! @brief The pattern of the positions of a list.
template <std::size_t Rows, std::size_t Columns, auto List>
[[nodiscard]] consteval auto to_pattern()
    -> sparse::pattern<Rows, Columns, List.size> {
  sparse::pattern<Rows, Columns, List.size> result;

  for (std::size_t i{0}; i < List.size; ++i) {
    result.positions[i] = List.items[i];
  }

  return result;
}

//! @brief Merge two patterns, keeping the positions of both, or only those in
//! common, in ascending order.
template <auto Lhs, auto Rhs, bool KeepAll>
[[nodiscard]] consteval auto merge_positions() {
  position_list<Lhs.count + Rhs.count> list;
  std::size_t i{0};
  std::size_t j{0};

  while (i < Lhs.count || j < Rhs.count) {
    if (j == Rhs.count ||
        (i < Lhs.count && Lhs.positions[i] < Rhs.positions[j])) {
      if constexpr (KeepAll) {
        list.push(Lhs.positions[i]);
      }
      ++i;
    } else if (i == Lhs.count || Rhs.positions[j] < Lhs.positions[i]) {
      if constexpr (KeepAll) {
        list.push(Rhs.positions[j]);
      }
      ++j;
    } else {
      list.push(Lhs.positions[i]);
      ++i;
      ++j;
    }
  }

  return list;
}

//! @brief The count of elementary products of a matrix product.
template <auto Lhs, auto Rhs>
inline constexpr std::size_t term_count{[] -> std::size_t {
  std::size_t result{0};

  for (std::size_t i{0}; i < Lhs.count; ++i) {
    const std::size_t row{Lhs.positions[i].column};

    result += Rhs.lower_bound({row + 1, 0}) - Rhs.lower_bound({row, 0});
  }

  return result;
}()};

//! @brief The positions reached by at least one term of a matrix product.
template <auto Lhs, auto Rhs> [[nodiscard]] consteval auto product_positions() {
  position_list<term_count<Lhs, Rhs>> list;

  for (std::size_t row{0}; row < Lhs.rows; ++row) {
    std::array<bool, Rhs.columns> reached{};
    const std::size_t end{Lhs.lower_bound({row + 1, 0})};

    for (std::size_t i{Lhs.lower_bound({row, 0})}; i < end; ++i) {
      const std::size_t inner{Lhs.positions[i].column};
      const std::size_t inner_end{Rhs.lower_bound({inner + 1, 0})};

      for (std::size_t j{Rhs.lower_bound({inner, 0})}; j < inner_end; ++j) {
        reached[Rhs.positions[j].column] = true;
      }
    }

    for (std::size_t column{0}; column < Rhs.columns; ++column) {
      if (reached[column]) {
        list.push({row, column});
      }
    }
  }

  return list;
}

//! @brief For each position of `To`, its stored element index in `From`, or
//! `absent`.
template <auto From, auto To>
inline constexpr auto projection{[] {
  std::array<std::size_t, To.count> result{};

  for (std::size_t i{0}; i < To.count; ++i) {
    result[i] = From.index_of(To.positions[i].row, To.positions[i].column);
  }

  return result;
}()};

//! @brief An elementary product of a matrix product, or of a dot product: the
//! index of the result element it accumulates into, and the indexes of its two
//! factors.
struct term {
  std::size_t result{0};
  std::size_t lhs{0};
  std::size_t rhs{0};
};

//! @brief The position along a vector of the given stored element.
template <auto Pattern>
[[nodiscard]] constexpr auto vector_index(std::size_t i) -> std::size_t {
  return (Pattern.positions[i].row * Pattern.columns) +
         Pattern.positions[i].column;
}

//! @brief The elementary products of the dot product of two vectors, of the
//! same size and of any orientation.
template <auto Lhs, auto Rhs>
inline constexpr auto dot_terms{[] {
  constexpr std::size_t count{[] -> std::size_t {
    std::size_t result{0};

    for (std::size_t i{0}, j{0}; (i < Lhs.count) && (j < Rhs.count);) {
      if (vector_index<Lhs>(i) < vector_index<Rhs>(j)) {
        ++i;
      } else if (vector_index<Rhs>(j) < vector_index<Lhs>(i)) {
        ++j;
      } else {
        ++result;
        ++i;
        ++j;
      }
    }

    return result;
  }()};
  std::array<term, count> result{};
  std::size_t index{0};

  for (std::size_t i{0}, j{0}; (i < Lhs.count) && (j < Rhs.count);) {
    if (vector_index<Lhs>(i) < vector_index<Rhs>(j)) {
      ++i;
    } else if (vector_index<Rhs>(j) < vector_index<Lhs>(i)) {
      ++j;
    } else {
      result[index] = {.result = index, .lhs = i, .rhs = j};
      ++index;
      ++i;
      ++j;
    }
  }

  return result;
}()};

//! @brief The absolute value, without a hosted implementation.
template <typename Type>
[[nodiscard]] constexpr auto absolute(const Type &value) -> Type {
  return value < Type{} ? -value : value;
}

//! @brief The square root, without a hosted implementation or at compile time
//! by Newton's method, descending from above the root.
//!
//! @note The result is not guaranteed to be identical to the last digit to
//! that of `std::sqrt`.
template <std::floating_point Type>
[[nodiscard]] constexpr auto square_root(Type value) -> Type {
#if __STDC_HOSTED__
  if (not std::is_constant_evaluated()) {
    return std::sqrt(value);
  }
#endif

  // A zero, a negative, a not-a-number, or an infinite value has no iteration.
  if (not(value > Type{0}) || value > std::numeric_limits<Type>::max()) {
    return value;
  }

  Type estimate{value < Type{1} ? Type{1} : value};
  Type next{(estimate + (value / estimate)) / Type{2}};

  // The estimates decrease strictly until the root, where the iteration stops.
  while (next < estimate) {
    estimate = next;
    next = (estimate + (value / estimate)) / Type{2};
  }

  return estimate;
}
} // namespace fcarouge::sparse_linear_algebra::internal

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP

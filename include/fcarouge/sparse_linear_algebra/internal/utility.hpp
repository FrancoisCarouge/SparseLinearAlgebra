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
#include <utility>

#if __STDC_HOSTED__

#include <cmath>

#endif

namespace fcarouge::sparse_linear_algebra::internal {
//! @brief Whether a type is a sparse matrix.
template <typename Type> inline constexpr bool is_matrix{false};

template <typename Type, sparsity_pattern auto Pattern>
inline constexpr bool is_matrix<sparse_matrix<Type, Pattern>>{true};

//! @brief The element type of a sparse matrix.
template <typename Matrix> struct element {};

template <typename Type, sparsity_pattern auto Pattern>
struct element<sparse_matrix<Type, Pattern>> {
  using type = Type;
};

//! @brief The element type of a sparse matrix.
template <typename Matrix> using element_t = element<Matrix>::type;

//! @brief Whether a pattern stores every position: a pattern without a
//! structural zero.
template <auto Pattern>
inline constexpr bool is_dense{Pattern.count == Pattern.rows * Pattern.columns};

//! @brief The type of the value of a structural zero in a construction from
//! every row, the result of the `_o` literal.
struct structural_zero {};

//! @brief Whether the characters of a literal spell a zero: zero digits, at
//! least one, and at most one decimal point.
template <char... Characters>
inline constexpr bool is_zero_literal{
    ((Characters == '0' or Characters == '.') and ...) and
    ((Characters == '0' ? 1 : 0) + ... + 0) >= 1 and
    ((Characters == '.' ? 1 : 0) + ... + 0) <= 1};

//! @brief Whether a value may be given at a linear position of a pattern: a
//! structural zero only takes the `_o` literal, a stored element only takes a
//! value of exactly the element type.
//!
//! @details Partially specialized variable templates rather than immediately
//! invoked templated lambdas, for the constraints' evaluation by MSVC.
template <typename Value, typename Type, auto Pattern, std::size_t Position,
          bool Stored = Pattern.stores(Position / Pattern.columns,
                                       Position % Pattern.columns)>
inline constexpr bool is_given_at{std::same_as<Value, Type>};

template <typename Value, typename Type, auto Pattern, std::size_t Position>
inline constexpr bool is_given_at<Value, Type, Pattern, Position, false>{
    std::same_as<Value, structural_zero>};

template <typename Type, auto Pattern, std::size_t Row, typename Columns,
          typename... Values>
inline constexpr bool are_given_in_row_at{false};

template <typename Type, auto Pattern, std::size_t Row, std::size_t... Columns,
          typename... Values>
inline constexpr bool are_given_in_row_at<
    Type, Pattern, Row, std::index_sequence<Columns...>, Values...>{
    (is_given_at<Values, Type, Pattern, Row * Pattern.columns + Columns> and
     ...)};

//! @brief Whether the values may be given at every column of a row of a
//! pattern, in order.
template <typename Type, auto Pattern, std::size_t Row, typename... Values>
inline constexpr bool are_given_in_row{
    are_given_in_row_at<Type, Pattern, Row, std::index_sequence_for<Values...>,
                        Values...>};

//! @brief The count of stored elements of a row of a pattern.
template <auto Pattern, std::size_t Row>
inline constexpr std::size_t stored_in_row{[] consteval -> std::size_t {
  constexpr std::size_t stored{Pattern.count};
  std::size_t count{0};

  for (std::size_t i{0}; i < stored; ++i) {
    if (Pattern.positions[i].row == Row) {
      ++count;
    }
  }

  return count;
}()};

//! @brief Write a given value into the next stored element, unless it is a
//! structural zero: the stored elements are the given values in row-major
//! order.
template <typename Type, std::size_t Count, typename Value>
constexpr void store_given(std::array<Type, Count> &stored, std::size_t &index,
                           const Value &value) {
  if constexpr (not std::same_as<Value, structural_zero>) {
    stored[index] = value;
    ++index;
  }
}

//! @brief Append the stored elements of a row to the stored elements of a
//! matrix.
template <typename Type, std::size_t Count, std::size_t RowCount>
constexpr void append_row(std::array<Type, Count> &stored, std::size_t &index,
                          const std::array<Type, RowCount> &elements) {
  for (const Type &element : elements) {
    stored[index] = element;
    ++index;
  }
}

//! @brief The stored elements of a row of a matrix, given at every column: a
//! value of exactly the element type at a stored element, the `0_o` literal at
//! a structural zero.
//!
//! @details Implicitly constructed from a braced row of a matrix declaration.
template <typename Type, auto Pattern, std::size_t Row> struct row {
  //! @brief The stored elements of the row, in column order.
  std::array<Type, stored_in_row<Pattern, Row>> elements{};

  //! @brief Construct from the value of every column of the row.
  template <typename... Values>
    requires(sizeof...(Values) == Pattern.columns and
             are_given_in_row<Type, Pattern, Row, Values...>)
  constexpr explicit(false) row(const Values &...values) {
    std::size_t index{0};

    (store_given(elements, index, values), ...);
  }
};

//! @brief The count of rows of a pattern.
//!
//! @details Named, for CppCheck, which misparses a member access in a template
//! argument list.
template <auto Pattern> inline constexpr std::size_t row_count{Pattern.rows};

//! @brief The storage of a matrix, and its construction from every row.
//!
//! @details A base of the matrix: the parameters of the construction from every
//! row are of a type each, the pack expansion of the rows of the pattern.
template <typename Type, auto Pattern,
          typename Rows = std::make_index_sequence<row_count<Pattern>>>
struct storage;

template <typename Type, auto Pattern, std::size_t... Rows>
struct storage<Type, Pattern, std::index_sequence<Rows...>> {
  //! @brief The stored elements, zero-initialized, in row-major order.
  std::array<Type, Pattern.count> values{};

  //! @brief Construct a zero matrix.
  constexpr storage() = default;

  //! @brief Construct from every row, in order.
  //!
  //! @param rows The rows, each of the value of every column.
  constexpr explicit storage(const row<Type, Pattern, Rows> &...rows) {
    std::size_t index{0};

    (append_row(values, index, rows.elements), ...);
  }
};

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
//!
//! @details The list is in row-major order: the construction of the pattern
//! verifies it in linear time.
template <std::size_t Rows, std::size_t Columns, auto List>
[[nodiscard]] consteval auto to_pattern()
    -> sparse::pattern<Rows, Columns, List.size> {
  if constexpr (List.size == 0) {
    return sparse::pattern{sparse::shape<Rows, Columns>};
  } else {
    sparse::position items[List.size]{};

    for (std::size_t i{0}; i < List.size; ++i) {
      items[i] = List.items[i];
    }

    return sparse::pattern{sparse::shape<Rows, Columns>, items};
  }
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

//! @brief The positions of a transposed pattern, in ascending order.
//!
//! @details The positions are bucketed by column: ascending rows within each
//! column are already what the transposed row-major order requires.
template <auto Pattern> [[nodiscard]] consteval auto transpose_positions() {
  position_list<Pattern.count> list;

  for (std::size_t column{0}; column < Pattern.columns; ++column) {
    for (std::size_t i{0}; i < Pattern.count; ++i) {
      if (Pattern.positions[i].column == column) {
        list.push({column, Pattern.positions[i].row});
      }
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

//! @brief A maximum matching of the rows and the columns of a pattern: the row
//! matched to each column, or `absent`, and the count of matched pairs.
template <std::size_t Columns> struct matching {
  std::array<std::size_t, Columns> row{};
  std::size_t size{0};
};

//! @brief Match a row to a column, by an augmenting path: a breadth-first
//! search of the columns reachable from the row through the rows of the columns
//! already matched, until a column not matched, then the reassignment of the
//! columns of the path, each to the row it was reached from.
template <auto Pattern>
[[nodiscard]] consteval auto augment(std::size_t start,
                                     matching<Pattern.columns> &matched)
    -> bool {
  std::array<bool, Pattern.columns> visited{};
  std::array<std::size_t, Pattern.columns> previous{};
  std::array<std::size_t, Pattern.columns> queue{};
  std::size_t head{0};
  std::size_t tail{0};
  std::size_t row{start};
  std::size_t from{sparse::absent};

  while (row != sparse::absent) {
    const std::size_t end{Pattern.lower_bound({row + 1, 0})};

    for (std::size_t i{Pattern.lower_bound({row, 0})}; i < end; ++i) {
      const std::size_t column{Pattern.positions[i].column};

      if (not visited[column]) {
        visited[column] = true;
        previous[column] = from;

        if (matched.row[column] == sparse::absent) {
          for (std::size_t current{column}; current != sparse::absent;
               current = previous[current]) {
            matched.row[current] = previous[current] == sparse::absent
                                       ? start
                                       : matched.row[previous[current]];
          }

          return true;
        }

        queue[tail] = column;
        ++tail;
      }
    }

    row = sparse::absent;

    if (head < tail) {
      from = queue[head];
      row = matched.row[from];
      ++head;
    }
  }

  return false;
}

//! @brief A maximum matching of the rows and the columns of a pattern.
//!
//! @details The stored diagonal positions are matched first, so that a pattern
//! storing its diagonal keeps it, then each row left is matched by an
//! augmenting path.
template <auto Pattern>
[[nodiscard]] consteval auto maximum_matching() -> matching<Pattern.columns> {
  matching<Pattern.columns> result{};
  std::array<bool, Pattern.rows> placed{};

  for (std::size_t column{0}; column < Pattern.columns; ++column) {
    result.row[column] = sparse::absent;
  }

  for (std::size_t row{0}; row < Pattern.rows; ++row) {
    if (row < Pattern.columns && Pattern.stores(row, row)) {
      result.row[row] = row;
      placed[row] = true;
      ++result.size;
    }
  }

  for (std::size_t row{0}; row < Pattern.rows; ++row) {
    if (not placed[row] && augment<Pattern>(row, result)) {
      ++result.size;
    }
  }

  return result;
}

//! @brief The block triangular form of the system `x * rhs = lhs`, of a square
//! denominator of full structural rank, read as the system `rhsᵀ * xᵀ = lhsᵀ`
//! of an equation per column of the denominator and an unknown per row.
//!
//! @details Each unknown is matched to an equation storing it, the diagonal of
//! the permuted system. An unknown depends on the other unknowns its matched
//! equation stores. The blocks are the strongly connected components of these
//! dependencies, listed in the order of their solution: a block only depends
//! on the blocks before it.
template <std::size_t Size> struct block_triangular {
  //! @brief The equation, the column of the denominator, of each unknown.
  std::array<std::size_t, Size> equation{};

  //! @brief The unknowns, block by block, in the order of their solution.
  std::array<std::size_t, Size> order{};

  //! @brief The index of each unknown in the order.
  std::array<std::size_t, Size> rank{};

  //! @brief The block of each unknown.
  std::array<std::size_t, Size> block{};

  //! @brief The index in the order of the first unknown of each block, and the
  //! end of the last block.
  std::array<std::size_t, Size + 1> begin{};

  //! @brief The count of blocks.
  std::size_t blocks{0};
};

//! @brief The state of the search for the strongly connected components of the
//! dependencies of the unknowns: Tarjan's algorithm.
template <std::size_t Size> struct component_search {
  std::array<std::size_t, Size> index{};
  std::array<std::size_t, Size> lowest{};
  std::array<bool, Size> stacked{};
  std::array<std::size_t, Size> stack{};
  std::size_t stack_size{0};
  std::size_t visits{0};
  std::size_t ordered{0};
};

//! @brief Number an unknown in the order of its visit, and stack it.
template <std::size_t Size>
consteval void enter_component(std::size_t unknown,
                               component_search<Size> &search) {
  ++search.visits;
  search.index[unknown] = search.visits;
  search.lowest[unknown] = search.visits;
  search.stack[search.stack_size] = unknown;
  ++search.stack_size;
  search.stacked[unknown] = true;
}

//! @brief Close the block of an unknown, the first visited of its block: its
//! unknowns are the stacked ones down to it, ordered ascending.
template <std::size_t Size>
consteval void close_component(std::size_t unknown,
                               component_search<Size> &search,
                               block_triangular<Size> &result) {
  const std::size_t first{search.ordered};
  std::size_t member{sparse::absent};

  result.begin[result.blocks] = first;

  while (member != unknown) {
    --search.stack_size;
    member = search.stack[search.stack_size];
    search.stacked[member] = false;
    result.order[search.ordered] = member;
    result.block[member] = result.blocks;
    ++search.ordered;
  }

  // The unknowns of a block in ascending order: the elimination of a block
  // follows the order of the system.
  for (std::size_t i{first + 1}; i < search.ordered; ++i) {
    const std::size_t moved{result.order[i]};
    std::size_t j{i};

    while (j > first && result.order[j - 1] > moved) {
      result.order[j] = result.order[j - 1];
      --j;
    }

    result.order[j] = moved;
  }

  for (std::size_t i{first}; i < search.ordered; ++i) {
    result.rank[result.order[i]] = i;
  }

  ++result.blocks;
  result.begin[result.blocks] = search.ordered;
}

//! @brief Visit an unknown and its dependencies, the columns of the transposed
//! denominator its matched equation stores, depth first, closing a block at
//! its first visited unknown once its dependencies are closed.
//!
//! @details The depth-first search keeps its own stack of the unknowns being
//! visited, and of the next dependency of each, rather than recursing.
template <auto Transposed>
consteval void visit_component(std::size_t root,
                               component_search<Transposed.rows> &search,
                               block_triangular<Transposed.rows> &result) {
  std::array<std::size_t, Transposed.rows> visiting{};
  std::array<std::size_t, Transposed.rows> next{};
  std::size_t depth{1};

  enter_component(root, search);
  visiting[0] = root;
  next[0] = Transposed.lower_bound({result.equation[root], 0});

  while (depth > 0) {
    const std::size_t unknown{visiting[depth - 1]};
    const std::size_t end{
        Transposed.lower_bound({result.equation[unknown] + 1, 0})};

    if (next[depth - 1] < end) {
      const std::size_t dependency{
          Transposed.positions[next[depth - 1]].column};

      ++next[depth - 1];

      if (search.index[dependency] == 0) {
        enter_component(dependency, search);
        visiting[depth] = dependency;
        next[depth] = Transposed.lower_bound({result.equation[dependency], 0});
        ++depth;
      } else if (search.stacked[dependency] &&
                 search.index[dependency] < search.lowest[unknown]) {
        search.lowest[unknown] = search.index[dependency];
      }
    } else {
      if (search.lowest[unknown] == search.index[unknown]) {
        close_component(unknown, search, result);
      }

      --depth;

      if (depth > 0 &&
          search.lowest[unknown] < search.lowest[visiting[depth - 1]]) {
        search.lowest[visiting[depth - 1]] = search.lowest[unknown];
      }
    }
  }
}

//! @brief The block triangular form of the division by a square matrix of the
//! pattern, of full structural rank.
template <auto Pattern>
inline constexpr auto block_triangular_form{[] consteval {
  constexpr auto transposed{to_pattern<Pattern.columns, Pattern.rows,
                                       transpose_positions<Pattern>()>()};
  constexpr auto matched{maximum_matching<Pattern>()};
  block_triangular<Pattern.rows> result{};
  component_search<Pattern.rows> search{};

  for (std::size_t column{0}; column < Pattern.columns; ++column) {
    result.equation[matched.row[column]] = column;
  }

  for (std::size_t unknown{0}; unknown < Pattern.rows; ++unknown) {
    if (search.index[unknown] == 0) {
      visit_component<transposed>(unknown, search, result);
    }
  }

  return result;
}()};

//! @brief The rows of the quotient `x` of `x * rhs = lhs` storing each of its
//! columns, the unknowns.
//!
//! @details An unknown is reached from the columns of the numerator its matched
//! equation stores, and from the unknowns it depends on. The unknowns of a
//! block depend on each other: they are reached alike.
template <auto Lhs, auto Rhs>
inline constexpr auto quotient_reach{[] {
  constexpr auto form{block_triangular_form<Rhs>};
  constexpr auto transposed{
      to_pattern<Rhs.columns, Rhs.rows, transpose_positions<Rhs>()>()};
  constexpr auto lhs_transposed{
      to_pattern<Lhs.columns, Lhs.rows, transpose_positions<Lhs>()>()};
  std::array<std::array<bool, Lhs.rows>, Rhs.rows> result{};

  for (std::size_t block{0}; block < form.blocks; ++block) {
    std::array<bool, Lhs.rows> reached{};

    for (std::size_t k{form.begin[block]}; k < form.begin[block + 1]; ++k) {
      const std::size_t equation{form.equation[form.order[k]]};
      const std::size_t lhs_end{lhs_transposed.lower_bound({equation + 1, 0})};
      const std::size_t end{transposed.lower_bound({equation + 1, 0})};

      for (std::size_t i{lhs_transposed.lower_bound({equation, 0})};
           i < lhs_end; ++i) {
        reached[lhs_transposed.positions[i].column] = true;
      }

      for (std::size_t i{transposed.lower_bound({equation, 0})}; i < end; ++i) {
        const std::size_t dependency{transposed.positions[i].column};

        for (std::size_t row{0}; row < Lhs.rows; ++row) {
          reached[row] = reached[row] || result[dependency][row];
        }
      }
    }

    for (std::size_t k{form.begin[block]}; k < form.begin[block + 1]; ++k) {
      result[form.order[k]] = reached;
    }
  }

  return result;
}()};

//! @brief The positions of the quotient `x` of `x * rhs = lhs`, in ascending
//! order.
template <auto Lhs, auto Rhs>
[[nodiscard]] consteval auto quotient_positions() {
  constexpr auto reached{quotient_reach<Lhs, Rhs>};
  position_list<Lhs.rows * Rhs.columns> list;

  for (std::size_t row{0}; row < Lhs.rows; ++row) {
    for (std::size_t column{0}; column < Rhs.columns; ++column) {
      if (reached[column][row]) {
        list.push({row, column});
      }
    }
  }

  return list;
}

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

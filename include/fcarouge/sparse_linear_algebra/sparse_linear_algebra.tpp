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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_SPARSE_LINEAR_ALGEBRA_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_SPARSE_LINEAR_ALGEBRA_TPP

//! @file
//! @brief The member definitions of the compile-time sparsity pattern matrix.

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra/pattern.hpp"
#include "fcarouge/sparse_linear_algebra/pattern_algebra.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace fcarouge {
template <typename Type, auto Pattern>
  requires sparsity_pattern<decltype(Pattern)>
template <std::same_as<Type>... Values>
  requires(sizeof...(Values) == Pattern.count and Pattern.count > 0)
constexpr sparse_matrix<Type, Pattern>::sparse_matrix(const Values &...stored) {
  values = std::array<Type, Pattern.count>{stored...};
}

template <typename Type, auto Pattern>
  requires sparsity_pattern<decltype(Pattern)>
template <auto Other>
  requires(sparsity_pattern<decltype(Other)> and
           sparse::is_subset_of<Other, Pattern>)
constexpr sparse_matrix<Type, Pattern>::sparse_matrix(
    const sparse_matrix<Type, Other> &other) {
  for (std::size_t i{0}; i < Pattern.count; ++i) {
    const std::size_t index{
        sparse_linear_algebra::internal::projection<Other, Pattern>[i]};

    if (index != sparse::absent) {
      values[i] = other.values[index];
    }
  }
}

template <typename Type, auto Pattern>
  requires sparsity_pattern<decltype(Pattern)>
constexpr auto
sparse_matrix<Type, Pattern>::operator[](std::size_t row,
                                         std::size_t column) const -> Type {
  if constexpr (sparse_linear_algebra::internal::is_dense<Pattern>) {
    return values[(row * Pattern.columns) + column];
  } else {
    const std::size_t index{Pattern.index_of(row, column)};

    return index == sparse::absent ? Type{} : values[index];
  }
}

template <typename Type, auto Pattern>
  requires sparsity_pattern<decltype(Pattern)>
constexpr auto sparse_matrix<Type, Pattern>::operator[](std::size_t row,
                                                        std::size_t column)
    -> Type &
  requires sparse_linear_algebra::internal::is_dense<Pattern>
{
  return values[(row * Pattern.columns) + column];
}

template <typename Type, auto Pattern>
  requires sparsity_pattern<decltype(Pattern)>
template <std::size_t Row, std::size_t Column>
  requires(Pattern.stores(Row, Column))
constexpr auto sparse_matrix<Type, Pattern>::operator[](
    std::integral_constant<std::size_t, Row> /*row*/,
    std::integral_constant<std::size_t, Column> /*column*/) -> Type & {
  return at<Row, Column>();
}

template <typename Type, auto Pattern>
  requires sparsity_pattern<decltype(Pattern)>
template <std::size_t Row, std::size_t Column>
  requires(Pattern.stores(Row, Column))
constexpr auto sparse_matrix<Type, Pattern>::at() -> Type & {
  static constexpr std::size_t index{Pattern.index_of(Row, Column)};

  return values[index];
}

template <typename Type, auto Pattern>
  requires sparsity_pattern<decltype(Pattern)>
template <std::size_t Row, std::size_t Column>
  requires(Row < Pattern.rows and Column < Pattern.columns)
constexpr auto sparse_matrix<Type, Pattern>::at() const -> Type {
  static constexpr std::size_t index{Pattern.index_of(Row, Column)};

  if constexpr (index == sparse::absent) {
    return Type{};
  } else {
    return values[index];
  }
}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_SPARSE_LINEAR_ALGEBRA_TPP

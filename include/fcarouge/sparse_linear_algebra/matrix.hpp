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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_MATRIX_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_MATRIX_HPP

//! @file
//! @brief The compile-time sparsity pattern matrix.

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra/pattern_algebra.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace fcarouge::sparse {
//! @name Types
//! @{

//! @brief Compile-time sparsity pattern matrix.
//!
//! @details Stores only the elements of the pattern, in row-major order, in
//! place. Every other element is a structural zero: it occupies no memory and
//! costs no operation. A read of a structural zero yields a zero value. A
//! structural zero has no element to refer to, so no write can reach it and no
//! value can be lost: a write to a structural zero does not compile. A pattern
//! without a structural zero behaves as any dense storage.
//!
//! @tparam Type The matrix element type.
//! @tparam Pattern The sparsity pattern of the matrix.
template <typename Type, sparsity_pattern auto Pattern> struct matrix {
  //! @brief The sparsity pattern of the matrix.
  static constexpr auto structure{Pattern};

  //! @brief The stored elements, zero-initialized, in row-major order.
  std::array<Type, Pattern.count> values{};

  //! @brief Construct a zero matrix.
  constexpr matrix() = default;

  //! @brief Construct from every stored element, in row-major order.
  //!
  //! @details Each element is exactly of the element type, and there is at
  //! least one: a structural zero takes no value.
  //!
  //! @param stored The stored elements, in row-major order.
  template <std::same_as<Type>... Values>
    requires(sizeof...(Values) == Pattern.count and Pattern.count > 0)
  constexpr explicit matrix(const Values &...stored) : values{stored...} {}

  //! @brief Convert from a matrix whose stored positions are all stored by
  //! this pattern, which cannot lose an element.
  //!
  //! @details The positions this pattern stores and the other does not are
  //! zero.
  //!
  //! @param other The matrix of a narrower, or the same, pattern.
  template <sparsity_pattern auto Other>
    requires is_subset_of<Other, Pattern>
  constexpr explicit(false) matrix(const matrix<Type, Other> &other) {
    for (std::size_t i{0}; i < Pattern.count; ++i) {
      const std::size_t index{sla::projection<Other, Pattern>[i]};

      if (index != absent) {
        values[i] = other.values[index];
      }
    }
  }

  //! @brief Read the element at the given position, by value.
  //!
  //! @details A structural zero reads as a zero of the element type.
  //!
  //! @param row The row of the element.
  //! @param column The column of the element.
  //!
  //! @complexity Logarithmic in the count of stored elements.
  [[nodiscard]] constexpr auto operator[](std::size_t row,
                                          std::size_t column) const -> Type {
    const std::size_t index{Pattern.index_of(row, column)};

    return index == absent ? Type{} : values[index];
  }

  //! @brief The element at the given position, of a pattern without a
  //! structural zero: every position is addressable, as in any dense storage.
  //!
  //! @param row The row of the element.
  //! @param column The column of the element.
  //!
  //! @complexity Constant.
  [[nodiscard]] constexpr auto operator[](std::size_t row, std::size_t column)
      -> Type &
    requires sla::is_dense<Pattern>
  {
    return values[(row * Pattern.columns) + column];
  }

  //! @brief The stored element at the given compile-time position.
  //!
  //! @details A structural zero has nothing to refer to, and no such access.
  //! The position's type is the standard type of the compile-time index of a
  //! linear algebra library, for a generic access by position.
  //!
  //! @tparam Row The row of the element.
  //! @tparam Column The column of the element.
  //!
  //! @complexity Constant.
  template <std::size_t Row, std::size_t Column>
    requires(Pattern.stores(Row, Column))
  [[nodiscard]] constexpr auto
  operator[](std::integral_constant<std::size_t, Row> /*row*/,
             std::integral_constant<std::size_t, Column> /*column*/) -> Type & {
    return at<Row, Column>();
  }

  //! @brief The stored element at the given compile-time position.
  //!
  //! @details A structural zero has nothing to refer to, and no such access:
  //! the read-only overload reads it, by value.
  //!
  //! @tparam Row The row of the element.
  //! @tparam Column The column of the element.
  //!
  //! @complexity Constant.
  template <std::size_t Row, std::size_t Column>
    requires(Pattern.stores(Row, Column))
  [[nodiscard]] constexpr auto at() -> Type & {
    static constexpr std::size_t index{Pattern.index_of(Row, Column)};

    return values[index];
  }

  //! @brief Read the element at the given compile-time position, by value.
  //!
  //! @details A structural zero reads as a zero of the element type. Selected
  //! for a structural zero of a modifiable matrix too: the read is a value, not
  //! a reference, so that an assignment to the element of an arithmetic type, a
  //! write that would be lost, does not compile.
  //!
  //! @tparam Row The row of the element.
  //! @tparam Column The column of the element.
  //!
  //! @complexity Constant.
  template <std::size_t Row, std::size_t Column>
    requires(Row < Pattern.rows and Column < Pattern.columns)
  [[nodiscard]] constexpr auto at() const -> Type {
    static constexpr std::size_t index{Pattern.index_of(Row, Column)};

    if constexpr (index == absent) {
      return Type{};
    } else {
      return values[index];
    }
  }
};

//! @}
} // namespace fcarouge::sparse

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_MATRIX_HPP

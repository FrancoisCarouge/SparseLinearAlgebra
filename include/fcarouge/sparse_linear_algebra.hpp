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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_HPP

//! @file
//! @brief Compile-time sparsity pattern linear algebra.
//!
//! @details The positions of the stored elements are part of the matrix type.
//! Every other element is a structural zero: it occupies no memory and costs no
//! operation. Each operation computes the pattern of its result, and the exact
//! arithmetic terms it needs, at compile time. There is no memory allocation,
//! no exception, no input-output, and every operation is `constexpr`. A pattern
//! without a structural zero behaves as any dense storage.
//!
//! @note Freestanding compatible. The formatting support requires a hosted
//! implementation.
//!
//! @note A structural zero stays an exact zero: no term exists for it. Unlike
//! a dense computation, a non-finite operand does not propagate through it, for
//! example `0 * infinity` is never computed, and is not a not-a-number value.

#include "sparse_linear_algebra/internal/utility.hpp"
#include "sparse_linear_algebra/pattern.hpp"
#include "sparse_linear_algebra/pattern_algebra.hpp"
#include "sparse_linear_algebra_forward.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace fcarouge::sparse::literals {
//! @name Literals
//! @{

//! @brief The structural zero literal, the value of a structural zero in a
//! construction from every row: `0_o`, or `0._o`.
//!
//! @details Does not participate in overload resolution unless the literal is
//! a zero: `1_o` does not compile.
template <char... Characters>
  requires sparse_linear_algebra::internal::is_zero_literal<Characters...>
[[nodiscard]] consteval auto operator""_o() noexcept
    -> sparse_linear_algebra::internal::structural_zero {
  return {};
}

//! @}
} // namespace fcarouge::sparse::literals

namespace fcarouge {
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
template <typename Type, auto Pattern>
  requires sparse::sparsity_pattern<decltype(Pattern)>
class sparse_matrix
    : public sparse_linear_algebra::internal::storage<Type, Pattern> {
public:
  //! @name Public Member Variables
  //! @{

  //! @brief The sparsity pattern of the matrix.
  static constexpr auto structure{Pattern};

  //! @brief The stored elements, zero-initialized, in row-major order.
  using sparse_linear_algebra::internal::storage<Type, Pattern>::values;

  //! @}

  //! @name Public Member Functions
  //! @{

  //! @brief Construct a zero matrix.
  constexpr sparse_matrix() = default;

  //! @brief Construct from every row, in order, each a braced list of the value
  //! of every column: the matrix as it reads.
  //!
  //! @details Each stored element is exactly of the element type, and each
  //! structural zero is the `0_o` literal: a value at a structural zero, even
  //! a zero, does not compile, nor does `0_o` at a stored element.
  //!
  //! @complexity Linear in the count of positions.
  using sparse_linear_algebra::internal::storage<Type, Pattern>::storage;

  //! @brief Construct from every stored element, in row-major order.
  //!
  //! @details Each element is exactly of the element type, and there is at
  //! least one: a structural zero takes no value.
  //!
  //! @param stored The stored elements, in row-major order.
  template <std::same_as<Type>... Values>
    requires(sizeof...(Values) == Pattern.count and Pattern.count > 0)
  constexpr explicit sparse_matrix(const Values &...stored);

  //! @brief Convert from a matrix whose stored positions are all stored by
  //! this pattern, which cannot lose an element.
  //!
  //! @details The positions this pattern stores and the other does not are
  //! zero.
  //!
  //! @param other The matrix of a narrower, or the same, pattern.
  template <auto Other>
    requires(sparse::sparsity_pattern<decltype(Other)> and
             sparse::is_subset_of<Other, Pattern>)
  constexpr explicit(false)
      sparse_matrix(const sparse_matrix<Type, Other> &other);

  //! @brief Read the element at the given position, by value.
  //!
  //! @details A structural zero reads as a zero of the element type.
  //!
  //! @param row The row of the element.
  //! @param column The column of the element.
  //!
  //! @complexity Logarithmic in the count of stored elements.
  [[nodiscard]] constexpr auto operator[](std::size_t row,
                                          std::size_t column) const -> Type;

  //! @brief The element at the given position, of a pattern without a
  //! structural zero: every position is addressable, as in any dense storage.
  //!
  //! @param row The row of the element.
  //! @param column The column of the element.
  //!
  //! @complexity Constant.
  [[nodiscard]] constexpr auto operator[](std::size_t row, std::size_t column)
      -> Type &
    requires sparse_linear_algebra::internal::is_dense<Pattern>;

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
             std::integral_constant<std::size_t, Column> /*column*/) -> Type &;

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
  [[nodiscard]] constexpr auto at() -> Type &;

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
  [[nodiscard]] constexpr auto at() const -> Type;

  //! @}
};

//! @}
} // namespace fcarouge

#include "sparse_linear_algebra/algorithm/add.tpp"
#include "sparse_linear_algebra/algorithm/divide.tpp"
#include "sparse_linear_algebra/algorithm/dot.tpp"
#include "sparse_linear_algebra/algorithm/equal_to.tpp"
#include "sparse_linear_algebra/algorithm/matrix_frob_norm.tpp"
#include "sparse_linear_algebra/algorithm/matrix_inf_norm.tpp"
#include "sparse_linear_algebra/algorithm/matrix_one_norm.tpp"
#include "sparse_linear_algebra/algorithm/minus.tpp"
#include "sparse_linear_algebra/algorithm/product.tpp"
#include "sparse_linear_algebra/algorithm/subtract.tpp"
#include "sparse_linear_algebra/algorithm/transposed.tpp"
#include "sparse_linear_algebra/format.tpp"
#include "sparse_linear_algebra/sparse_linear_algebra.tpp"

namespace fcarouge {
//! @name Algorithms
//! @{

template <typename Type, sparse::sparsity_pattern auto Lhs,
          sparse::sparsity_pattern auto Rhs>
  requires((Lhs.rows == 1 or Lhs.columns == 1) and
           (Rhs.rows == 1 or Rhs.columns == 1) and
           Lhs.rows * Lhs.columns == Rhs.rows * Rhs.columns)
[[nodiscard]] constexpr auto dot(const sparse_matrix<Type, Lhs> &lhs,
                                 const sparse_matrix<Type, Rhs> &rhs) -> Type;

template <std::floating_point Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
matrix_frob_norm(const sparse_matrix<Type, Pattern> &value) -> Type;

template <typename Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
matrix_inf_norm(const sparse_matrix<Type, Pattern> &value) -> Type;

template <typename Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
matrix_one_norm(const sparse_matrix<Type, Pattern> &value) -> Type;

template <typename Type, sparse::sparsity_pattern auto Lhs,
          sparse::sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator==(const sparse_matrix<Type, Lhs> &lhs,
                                        const sparse_matrix<Type, Rhs> &rhs)
    -> bool;

template <typename Type, sparse::sparsity_pattern auto Lhs,
          sparse::sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator+(const sparse_matrix<Type, Lhs> &lhs,
                                       const sparse_matrix<Type, Rhs> &rhs)
    -> sparse_matrix<Type, sparse::pattern_union<Lhs, Rhs>>;

template <typename Type, sparse::sparsity_pattern auto Lhs,
          sparse::sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator-(const sparse_matrix<Type, Lhs> &lhs,
                                       const sparse_matrix<Type, Rhs> &rhs)
    -> sparse_matrix<Type, sparse::pattern_union<Lhs, Rhs>>;

template <typename Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
operator-(const sparse_matrix<Type, Pattern> &value)
    -> sparse_matrix<Type, Pattern>;

template <typename Type, sparse::sparsity_pattern auto Lhs,
          sparse::sparsity_pattern auto Rhs>
  requires(Lhs.columns == Rhs.rows)
[[nodiscard]] constexpr auto operator*(const sparse_matrix<Type, Lhs> &lhs,
                                       const sparse_matrix<Type, Rhs> &rhs)
    -> sparse_matrix<Type, sparse::pattern_product<Lhs, Rhs>>;

template <typename Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator*(const sparse_matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> sparse_matrix<Type, Pattern>;

template <typename Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator*(const std::type_identity_t<Type> &lhs,
                                       const sparse_matrix<Type, Pattern> &rhs)
    -> sparse_matrix<Type, Pattern>;

template <typename Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator/(const sparse_matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> sparse_matrix<Type, Pattern>;

template <typename Type, sparse::sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
transposed(const sparse_matrix<Type, Pattern> &value)
    -> sparse_matrix<Type, sparse::pattern_transpose<Pattern>>;

//! @}
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_HPP

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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_MAKE_MATRIX_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_MAKE_MATRIX_HPP

//! @file
//! @brief The construction of a matrix from a value at every position.

#include "fcarouge/sparse_linear_algebra/internal/utility.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace fcarouge {
namespace sparse_linear_algebra::internal {
//! @brief The type of the value of a structural zero in `make_matrix`, the
//! result of the `_o` literal.
struct structural_zero {};

//! @brief Whether the characters of a literal spell a zero: zero digits, at
//! least one, and at most one decimal point.
template <char... Characters>
inline constexpr bool is_zero_literal{
    ((Characters == '0' or Characters == '.') and ...) and
    ((Characters == '0' ? 1 : 0) + ... + 0) >= 1 and
    ((Characters == '.' ? 1 : 0) + ... + 0) <= 1};
} // namespace sparse_linear_algebra::internal

namespace sparse::literals {
//! @name Element Construction
//! @{

//! @brief The structural zero literal, the value of a structural zero in
//! `make_matrix`: `0_o`, or `0._o`.
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
} // namespace sparse::literals

namespace sparse_linear_algebra::internal {
//! @brief Whether a value may be given at a linear position of a matrix: a
//! structural zero only takes the `_o` literal, a stored element only takes an
//! implicitly convertible value.
//!
//! @details Partially specialized variable templates rather than immediately
//! invoked templated lambdas, for the constraints' evaluation by MSVC.
template <typename Value, typename Matrix, std::size_t Position,
          bool Stored =
              Matrix::structure.stores(Position / Matrix::structure.columns,
                                       Position % Matrix::structure.columns)>
inline constexpr bool is_given_at{
    std::is_convertible_v<Value, element_t<Matrix>>};

template <typename Value, typename Matrix, std::size_t Position>
inline constexpr bool is_given_at<Value, Matrix, Position, false>{
    std::same_as<Value, structural_zero>};

template <typename Matrix, typename Positions, typename... Values>
inline constexpr bool are_given_at{false};

template <typename Matrix, std::size_t... Positions, typename... Values>
inline constexpr bool
    are_given_at<Matrix, std::index_sequence<Positions...>, Values...>{
        (is_given_at<Values, Matrix, Positions> and ...)};

//! @brief Whether the values may be given at every position of a matrix, in
//! row-major order.
template <typename Matrix, typename... Values>
inline constexpr bool are_given{
    are_given_at<Matrix, std::index_sequence_for<Values...>, Values...>};

//! @brief Write the given value of a linear position into its stored element,
//! if the position is stored.
template <std::size_t Position, typename Matrix, typename Value>
constexpr void write_given(Matrix &result, const Value &value) {
  if constexpr (not std::same_as<Value, structural_zero>) {
    static constexpr std::size_t index{
        Matrix::structure.index_of(Position / Matrix::structure.columns,
                                   Position % Matrix::structure.columns)};

    // The implicit conversion the constraints verified, stated: a conversion
    // to the element type may narrow, or promote, which braces reject or an
    // implicit conversion warns of.
    result.values[index] = static_cast<element_t<Matrix>>(value);
  }
}
} // namespace sparse_linear_algebra::internal

namespace sparse {
//! @name Element Construction
//! @{

//! @brief A matrix from a value at every position, in row-major order.
//!
//! @details Each position is known at compile time, and the type of each value
//! says what it is: a value for a structural zero, even a zero, or
//! `0_o` for a stored element, does not compile. A run-time value
//! needs no run-time check.
//!
//! @tparam Matrix The sparse matrix type to construct.
//! @tparam Values The types of the values.
//!
//! @param values The value of each position, in row-major order: a value
//! implicitly convertible to the element type at a stored element,
//! the `0_o` literal at a structural zero.
//!
//! @complexity Linear in the count of positions.
template <typename Matrix, typename... Values>
  requires sparse_linear_algebra::internal::is_matrix<Matrix> and
           (sizeof...(Values) ==
            Matrix::structure.rows * Matrix::structure.columns) and
           sparse_linear_algebra::internal::are_given<Matrix, Values...>
[[nodiscard]] constexpr auto make_matrix(const Values &...values) -> Matrix {
  Matrix result;

  [&result, &values...]<std::size_t... Positions>(
      std::index_sequence<Positions...>) -> void {
    (sparse_linear_algebra::internal::write_given<Positions>(result, values),
     ...);
  }(std::index_sequence_for<Values...>{});

  return result;
}

//! @}
} // namespace sparse
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_MAKE_MATRIX_HPP

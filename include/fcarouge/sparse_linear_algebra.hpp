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
#include "sparse_linear_algebra/make_matrix.hpp"
#include "sparse_linear_algebra/matrix.hpp"
#include "sparse_linear_algebra/pattern.hpp"
#include "sparse_linear_algebra/pattern_algebra.hpp"
#include "sparse_linear_algebra_forward.hpp"

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

#include <concepts>
#include <type_traits>

namespace fcarouge::sparse {
//! @name Algorithms
//! @{

template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires((Lhs.rows == 1 or Lhs.columns == 1) and
           (Rhs.rows == 1 or Rhs.columns == 1) and
           Lhs.rows * Lhs.columns == Rhs.rows * Rhs.columns)
[[nodiscard]] constexpr auto dot(const matrix<Type, Lhs> &lhs,
                                 const matrix<Type, Rhs> &rhs) -> Type;

template <std::floating_point Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto
matrix_frob_norm(const matrix<Type, Pattern> &value) -> Type;

template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto matrix_inf_norm(const matrix<Type, Pattern> &value)
    -> Type;

template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto matrix_one_norm(const matrix<Type, Pattern> &value)
    -> Type;

template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator==(const matrix<Type, Lhs> &lhs,
                                        const matrix<Type, Rhs> &rhs) -> bool;

template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator+(const matrix<Type, Lhs> &lhs,
                                       const matrix<Type, Rhs> &rhs)
    -> matrix<Type, pattern_union<Lhs, Rhs>>;

template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
[[nodiscard]] constexpr auto operator-(const matrix<Type, Lhs> &lhs,
                                       const matrix<Type, Rhs> &rhs)
    -> matrix<Type, pattern_union<Lhs, Rhs>>;

template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator-(const matrix<Type, Pattern> &value)
    -> matrix<Type, Pattern>;

template <typename Type, sparsity_pattern auto Lhs, sparsity_pattern auto Rhs>
  requires(Lhs.columns == Rhs.rows)
[[nodiscard]] constexpr auto operator*(const matrix<Type, Lhs> &lhs,
                                       const matrix<Type, Rhs> &rhs)
    -> matrix<Type, pattern_product<Lhs, Rhs>>;

template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator*(const matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> matrix<Type, Pattern>;

template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator*(const std::type_identity_t<Type> &lhs,
                                       const matrix<Type, Pattern> &rhs)
    -> matrix<Type, Pattern>;

template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto operator/(const matrix<Type, Pattern> &lhs,
                                       const std::type_identity_t<Type> &rhs)
    -> matrix<Type, Pattern>;

template <typename Type, sparsity_pattern auto Pattern>
[[nodiscard]] constexpr auto transposed(const matrix<Type, Pattern> &value)
    -> matrix<Type, pattern_transpose<Pattern>>;

//! @}
} // namespace fcarouge::sparse

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_HPP

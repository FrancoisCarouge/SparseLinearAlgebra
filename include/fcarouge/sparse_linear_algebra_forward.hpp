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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_FORWARD_HPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_FORWARD_HPP

//! @file
//! @brief Authoritative forward declaration header.
//!
//! @details Use this authoritative header to forward declare the types of this
//! project and avoid inconsistent declarations. The sparsity pattern concept
//! is defined here, since the matrix declaration is constrained by it.

#include <cstddef>
#include <type_traits>

namespace fcarouge {
namespace sparse {
struct position;
template <std::size_t Rows, std::size_t Columns, std::size_t Count>
struct pattern;
struct structural_zero_t;
} // namespace sparse

namespace sparse_linear_algebra::internal {
template <typename Type> inline constexpr bool is_pattern{false};

template <std::size_t Rows, std::size_t Columns, std::size_t Count>
inline constexpr bool is_pattern<sparse::pattern<Rows, Columns, Count>>{true};
} // namespace sparse_linear_algebra::internal

namespace sla = sparse_linear_algebra::internal;

namespace sparse {
//! @brief Concept of a compile-time sparsity pattern.
template <typename Type>
concept sparsity_pattern = sla::is_pattern<std::remove_cv_t<Type>>;

template <typename Type, sparsity_pattern auto Pattern> class sparse_matrix;
} // namespace sparse
} // namespace fcarouge

#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_FORWARD_HPP

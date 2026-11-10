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

#include "fcarouge/sparse_linear_algebra.hpp"

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace fcarouge::test {
namespace {
using sparse::pattern;
using sparse::shape;

using diagonal_matrix = sparse_matrix<double, sparse::diagonal<2>>;
inline constexpr pattern column_pattern{shape<2, 1>, {{1, 0}}};
using column = sparse_matrix<double, column_pattern>;
inline constexpr pattern row_pattern{shape<1, 2>, {{0, 1}}};
using row = sparse_matrix<double, row_pattern>;
using singleton = sparse_matrix<double, sparse::dense<1, 1>>;
using empty_singleton = sparse_matrix<double, sparse::zero<1, 1>>;

template <typename Matrix, std::size_t Row, std::size_t Column>
concept writable_at =
    std::same_as<decltype(std::declval<Matrix &>()
                              [std::integral_constant<std::size_t, Row>{},
                               std::integral_constant<std::size_t, Column>{}]),
                 double &> and
    std::same_as<decltype(std::declval<Matrix &>().template at<Row, Column>()),
                 double &>;

template <typename Matrix, std::size_t Row, std::size_t Column>
concept readable_at = requires(const Matrix value) {
  { value.template at<Row, Column>() } -> std::same_as<double>;
};

template <typename Matrix, std::size_t Row, std::size_t Column>
concept constant_readable_at = requires {
  static_cast<double (Matrix::*)(std::integral_constant<std::size_t, Row>,
                                 std::integral_constant<std::size_t, Column>)
                  const>(&Matrix::template operator[]<Row, Column>);
};

template <typename Matrix, std::size_t Row, std::size_t Column>
concept assignable_at =
    requires(Matrix value) { value.template at<Row, Column>() = 1.; };

//! @test Verifies every position can be read, a structural zero included, but
//! only a stored position has a reference for a write, at compile time.
static_assert(readable_at<diagonal_matrix, 0, 0>);
static_assert(readable_at<diagonal_matrix, 0, 1>);
static_assert(writable_at<diagonal_matrix, 0, 0>);
static_assert(writable_at<diagonal_matrix, 1, 1>);
static_assert(not writable_at<diagonal_matrix, 0, 1>);
static_assert(not writable_at<diagonal_matrix, 1, 0>);

//! @test Verifies every position has its own read by value at a compile-time
//! position, a structural zero included, rather than a conversion of the
//! position to a run-time one.
static_assert(constant_readable_at<diagonal_matrix, 0, 0>);
static_assert(constant_readable_at<diagonal_matrix, 0, 1>);
static_assert(constant_readable_at<column, 0, 0>);
static_assert(constant_readable_at<column, 1, 0>);
static_assert(constant_readable_at<singleton, 0, 0>);
static_assert(constant_readable_at<empty_singleton, 0, 0>);

//! @test Verifies a position out of the matrix is not accessible.
static_assert(not readable_at<diagonal_matrix, 2, 0>);
static_assert(not readable_at<diagonal_matrix, 0, 2>);
static_assert(not constant_readable_at<diagonal_matrix, 2, 0>);
static_assert(not constant_readable_at<diagonal_matrix, 0, 2>);

//! @test Verifies the rule for a vector.
static_assert(readable_at<column, 0, 0>);
static_assert(not writable_at<column, 0, 0>);
static_assert(writable_at<column, 1, 0>);
static_assert(readable_at<row, 0, 0>);
static_assert(not writable_at<row, 0, 0>);
static_assert(writable_at<row, 0, 1>);

//! @test Verifies the rule for the singleton.
static_assert(writable_at<singleton, 0, 0>);
static_assert(readable_at<empty_singleton, 0, 0>);
static_assert(not writable_at<empty_singleton, 0, 0>);

//! @test Verifies a read of a structural zero of a modifiable matrix is not
//! assignable: the write would be lost.
static_assert(assignable_at<diagonal_matrix, 0, 0>);
static_assert(not assignable_at<diagonal_matrix, 0, 1>);

//! @test Verifies a run-time position only has a reference without a
//! structural zero.
static_assert(
    std::same_as<decltype(std::declval<diagonal_matrix &>()[0, 1]), double>);
static_assert(
    std::same_as<decltype(std::declval<sparse_matrix<
                              double, sparse::dense<2, 2>> &>()[0, 1]),
                 double &>);
} // namespace
} // namespace fcarouge::test

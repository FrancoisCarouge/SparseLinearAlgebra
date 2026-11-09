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

namespace fcarouge::test {
namespace {
using sparse::banded;
using sparse::dense;
using sparse::diagonal;
using sparse::lower_triangular;
using sparse::pattern;
using sparse::shape;
using sparse::tridiagonal;
using sparse::upper_triangular;
using sparse::zero;

template <auto Pattern> using matrix = sparse_matrix<double, Pattern>;

inline constexpr pattern rectangular{shape<2, 3>, {{0, 2}, {1, 0}}};
inline constexpr pattern rectangular_transposed{shape<3, 2>, {{0, 1}, {2, 0}}};
inline constexpr pattern first_element{shape<3, 1>, {{0, 0}}};

//! @test Verifies the sum and the difference store the union of their terms.
static_assert(
    std::same_as<decltype(matrix<diagonal<3>>{} + matrix<tridiagonal<3>>{}),
                 matrix<tridiagonal<3>>>);
static_assert(
    std::same_as<decltype(matrix<tridiagonal<3>>{} - matrix<diagonal<3>>{}),
                 matrix<tridiagonal<3>>>);
static_assert(std::same_as<decltype(matrix<lower_triangular<3>>{} +
                                    matrix<upper_triangular<3>>{}),
                           matrix<dense<3, 3>>>);
static_assert(
    std::same_as<decltype(matrix<zero<3, 3>>{} + matrix<diagonal<3>>{}),
                 matrix<diagonal<3>>>);

//! @test Verifies the product stores the structural product of its factors.
static_assert(
    std::same_as<decltype(matrix<diagonal<3>>{} * matrix<diagonal<3>>{}),
                 matrix<diagonal<3>>>);
static_assert(
    std::same_as<decltype(matrix<diagonal<3>>{} * matrix<tridiagonal<3>>{}),
                 matrix<tridiagonal<3>>>);
static_assert(
    std::same_as<decltype(matrix<tridiagonal<4>>{} * matrix<tridiagonal<4>>{}),
                 matrix<banded<4, 4, 2, 2>>>);
static_assert(std::same_as<decltype(matrix<lower_triangular<3>>{} *
                                    matrix<lower_triangular<3>>{}),
                           matrix<lower_triangular<3>>>);
static_assert(
    std::same_as<decltype(matrix<zero<3, 3>>{} * matrix<dense<3, 3>>{}),
                 matrix<zero<3, 3>>>);

//! @test Verifies the transposition transposes the pattern, and the shape.
static_assert(std::same_as<decltype(transposed(matrix<lower_triangular<3>>{})),
                           matrix<upper_triangular<3>>>);
static_assert(std::same_as<decltype(transposed(matrix<rectangular>{})),
                           matrix<rectangular_transposed>>);

//! @test Verifies the negation and the scaling keep the pattern.
static_assert(
    std::same_as<decltype(-matrix<tridiagonal<3>>{}), matrix<tridiagonal<3>>>);
static_assert(std::same_as<decltype(matrix<tridiagonal<3>>{} * 2.),
                           matrix<tridiagonal<3>>>);
static_assert(std::same_as<decltype(2. * matrix<tridiagonal<3>>{}),
                           matrix<tridiagonal<3>>>);
static_assert(std::same_as<decltype(matrix<tridiagonal<3>>{} / 2.),
                           matrix<tridiagonal<3>>>);

//! @test Verifies the division by a matrix stores the structural quotient.
static_assert(
    std::same_as<decltype(matrix<tridiagonal<3>>{} / matrix<diagonal<3>>{}),
                 matrix<tridiagonal<3>>>);
static_assert(std::same_as<decltype(matrix<diagonal<3>>{} /
                                    matrix<lower_triangular<3>>{}),
                           matrix<lower_triangular<3>>>);
static_assert(
    std::same_as<decltype(matrix<diagonal<3>>{} / matrix<tridiagonal<3>>{}),
                 matrix<dense<3, 3>>>);

//! @test Verifies a matrix-vector product is sparse when the factors are.
static_assert(
    std::same_as<decltype(matrix<diagonal<3>>{} * matrix<first_element>{}),
                 matrix<first_element>>);
inline constexpr pattern first_rows{shape<3, 1>, {{0, 0}, {1, 0}}};
static_assert(
    std::same_as<decltype(matrix<tridiagonal<3>>{} * matrix<first_element>{}),
                 matrix<first_rows>>);

//! @test Verifies the element type is kept whatever the pattern.
static_assert(std::same_as<decltype(sparse_matrix<float, diagonal<3>>{} *
                                    sparse_matrix<float, dense<3, 3>>{}),
                           sparse_matrix<float, dense<3, 3>>>);
} // namespace
} // namespace fcarouge::test

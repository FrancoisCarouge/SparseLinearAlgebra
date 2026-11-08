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

namespace fcarouge::test {
namespace {
using sparse::diagonal;
using sparse::make_pattern;
using sparse::literals::operator""_o;
using sparse::tridiagonal;

using diagonal_matrix = sparse_matrix<double, diagonal<3>>;
using tridiagonal_matrix = sparse_matrix<double, tridiagonal<3>>;

inline constexpr diagonal_matrix lhs{2., 3., 4.};
inline constexpr tridiagonal_matrix rhs{1., 2., 3., 4., 5., 6., 7.};

//! @test Verifies each operation is available in a constant expression:
//! everything is known at compile time, down to the values.
inline constexpr auto sum{lhs + rhs};

static_assert(sum.at<0, 0>() == 3.);
static_assert(sum.at<0, 1>() == 2.);
static_assert(sum.at<0, 2>() == 0.);
static_assert(sum.at<2, 2>() == 11.);

inline constexpr auto difference{lhs - rhs};

static_assert(difference.at<0, 0>() == 1.);
static_assert(difference.at<1, 0>() == -3.);

inline constexpr auto negated{-rhs};

static_assert(negated.at<1, 2>() == -5.);
static_assert(negated.at<0, 2>() == 0.);

inline constexpr auto product{lhs * rhs};

static_assert(product.at<0, 1>() == 4.);
static_assert(product.at<2, 1>() == 24.);
static_assert(product.at<2, 0>() == 0.);

inline constexpr auto scaled{lhs * 2.};

static_assert(scaled.at<1, 1>() == 6.);
static_assert((2. * lhs).at<2, 2>() == 8.);
static_assert((lhs / 2.).at<2, 2>() == 2.);
static_assert(sum == rhs + lhs);

inline constexpr sparse_matrix<double, make_pattern<2, 3>({{0, 2}, {1, 0}})>
    rectangular{{0_o, 0_o, 5.},  //
                {6., 0_o, 0_o}}; //
inline constexpr auto flipped{transposed(rectangular)};

static_assert(flipped.at<2, 0>() == 5.);
static_assert(flipped.at<0, 1>() == 6.);
static_assert(flipped.at<1, 1>() == 0.);

//! @test Verifies the norms and the dot product in a constant expression.
static_assert(matrix_one_norm(rhs) == 12.);
static_assert(matrix_inf_norm(rhs) == 13.);
static_assert(matrix_frob_norm(
                  sparse_matrix<double, make_pattern<2, 2>({{0, 1}, {1, 0}})>{
                      {0_o, 3.},          //
                      {4., 0_o}}) == 5.); //

inline constexpr sparse_matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})>
    first_third{2., 3.};

static_assert(dot(first_third, first_third) == 13.);

//! @test Verifies writes in a constant expression.
inline constexpr auto written{[] -> diagonal_matrix {
  diagonal_matrix value;

  value.at<1, 1>() = 7.;

  return value;
}()};

static_assert(written.at<1, 1>() == 7.);
static_assert(written.at<0, 0>() == 0.);
} // namespace
} // namespace fcarouge::test

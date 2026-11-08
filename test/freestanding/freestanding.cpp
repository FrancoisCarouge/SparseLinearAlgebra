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
//! @test Verifies the library compiles for a freestanding implementation,
//! where the compiler supports one.
#ifndef _MSC_VER
static_assert(__STDC_HOSTED__ == 0);
#endif

using literals::operator""_o;

using band_matrix = sparse_matrix<double, sparse::tridiagonal<3>>;

inline constexpr band_matrix band{1., 2., 3., 4., 5., 6., 7.};

//! @test Verifies the operations without a hosted implementation, the
//! square root's included.
static_assert((band + band).at<0, 1>() == 4.);
static_assert((band - (band * 2.)).at<0, 1>() == -2.);
static_assert((-band).at<0, 1>() == -2.);
static_assert((band * band).at<1, 1>() == 52.);
static_assert((2. * band).at<2, 2>() == 14.);
static_assert((band / 2.).at<0, 1>() == 1.);
static_assert(transposed(band).at<1, 0>() == 2.);
static_assert(band == band_matrix{1., 2., 3., 4., 5., 6., 7.});
static_assert(matrix_one_norm(band) == 12.);
static_assert(matrix_inf_norm(band) == 13.);
static_assert(
    matrix_frob_norm(
        sparse_matrix<double, sparse::make_pattern<2, 2>({{0, 1}, {1, 0}})>{
            3., 4.}) == 5.);
static_assert(dot(sparse_matrix<double, sparse::dense<2, 1>>{3., 4.},
                  sparse_matrix<double, sparse::dense<1, 2>>{3., 4.}) == 25.);
static_assert(band_matrix{{1., 2., 0_o},           //
                          {3., 4., 5.},            //
                          {0_o, 6., 7.}} == band); //
} // namespace
} // namespace fcarouge::test

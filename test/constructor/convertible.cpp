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
using sparse::dense;
using sparse::diagonal;
using sparse::sparse_matrix;
using sparse::tridiagonal;

//! @test Verifies a matrix converts implicitly when it cannot lose an element:
//! to a pattern storing every one of its positions.
static_assert(std::convertible_to<sparse_matrix<double, diagonal<2>>,
                                  sparse_matrix<double, dense<2, 2>>>);
static_assert(std::convertible_to<sparse_matrix<double, diagonal<3>>,
                                  sparse_matrix<double, tridiagonal<3>>>);
static_assert(std::convertible_to<sparse_matrix<double, sparse::zero<3, 3>>,
                                  sparse_matrix<double, diagonal<3>>>);

//! @test Verifies a matrix does not convert to a narrower pattern, nor to
//! another shape, nor to another element type.
static_assert(not std::constructible_from<sparse_matrix<double, diagonal<2>>,
                                          sparse_matrix<double, dense<2, 2>>>);
static_assert(not std::constructible_from<sparse_matrix<double, diagonal<3>>,
                                          sparse_matrix<double, diagonal<2>>>);
static_assert(not std::constructible_from<sparse_matrix<double, dense<2, 2>>,
                                          sparse_matrix<float, diagonal<2>>>);

//! @test Verifies the construction from the stored elements is explicit, of
//! exactly the element type, and of every stored element.
static_assert(std::constructible_from<sparse_matrix<double, diagonal<2>>,
                                      double, double>);
static_assert(not std::convertible_to<
              double, sparse_matrix<double, sparse::dense<1, 1>>>);
static_assert(not std::constructible_from<sparse_matrix<double, diagonal<2>>,
                                          float, float>);
static_assert(not std::constructible_from<sparse_matrix<double, diagonal<2>>,
                                          double, double, double>);
static_assert(
    not std::constructible_from<sparse_matrix<double, diagonal<2>>, double>);
} // namespace
} // namespace fcarouge::test

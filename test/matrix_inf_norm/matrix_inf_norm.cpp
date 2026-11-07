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

#include <cassert>

namespace fcarouge::test {
namespace {
//! @test Verifies the infinity norm, the largest absolute row sum. The operands
//! discriminate a wrong implementation: negative elements for an absolute
//! value, a non-square matrix for a row and column mix-up, a row vector for a
//! vector norm special case, a position shared between rows and columns.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_pattern;
  using sparse::matrix;

  // The largest absolute row sum is not the largest absolute column sum.
  const matrix<double, make_pattern<2, 3>({{0, 0}, {0, 1}, {1, 0}, {1, 2}})>
      value{1., -2., -3., 4.};

  assert(matrix_inf_norm(value) == 7.);

  // A row vector: the infinity norm is the sum of its absolute elements.
  const matrix<double, make_pattern<1, 3>({{0, 0}, {0, 2}})> row{-3., 4.};

  assert(matrix_inf_norm(row) == 7.);

  // A column vector: the infinity norm is its largest absolute element.
  const matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})> column{-3., 4.};

  assert(matrix_inf_norm(column) == 4.);

  // A row without any stored element.
  const matrix<double, make_pattern<2, 2>({{0, 0}})> corner{-2.};

  assert(matrix_inf_norm(corner) == 2.);

  // Every negative row: the maximum is not of the signed sums.
  const matrix<double, make_pattern<2, 2>({{1, 0}, {1, 1}})> last_row{-1., -2.};

  assert(matrix_inf_norm(last_row) == 3.);

  // The zero matrix.
  const matrix<double, sparse::zero<2, 2>> none;

  assert(matrix_inf_norm(none) == 0.);

  // An integer element type, needing no square root.
  const matrix<int, make_pattern<2, 2>({{0, 0}, {1, 0}})> integers{-2, 3};

  assert(matrix_inf_norm(integers) == 3);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

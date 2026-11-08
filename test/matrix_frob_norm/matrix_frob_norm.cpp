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
#include <cmath>

namespace fcarouge::test {
namespace {
//! @test Verifies the Frobenius norm, of the stored elements only, of every
//! shape.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::pattern;
  using sparse::shape;

  constexpr pattern value_pattern{shape<2, 3>,
                                  {{0, 0}, {0, 1}, {1, 0}, {1, 2}}};
  const sparse_matrix<double, value_pattern> value{1., -2., -3., 4.};
  constexpr pattern row_pattern{shape<1, 3>, {{0, 0}, {0, 2}}};
  const sparse_matrix<double, row_pattern> row{-3., 4.};
  constexpr pattern column_pattern{shape<3, 1>, {{0, 0}, {2, 0}}};
  const sparse_matrix<double, column_pattern> column{-3., 4.};
  constexpr pattern corner_pattern{shape<2, 2>, {{0, 0}}};
  const sparse_matrix<double, corner_pattern> corner{-2.};
  const sparse_matrix<double, sparse::zero<2, 2>> none;
  const sparse_matrix<float, sparse::diagonal<2>> single_precision{3.F, 4.F};

  assert(matrix_frob_norm(value) == std::sqrt(30.));
  assert(matrix_frob_norm(row) == 5.);
  assert(matrix_frob_norm(column) == 5.);
  assert(matrix_frob_norm(corner) == 2.);
  assert(matrix_frob_norm(none) == 0.);
  assert(matrix_frob_norm(single_precision) == 5.F);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

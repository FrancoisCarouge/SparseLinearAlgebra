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
//! @test Verifies the assignments of a matrix, from the same pattern and from a
//! narrower pattern.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::diagonal;
  using sparse::matrix;

  using diagonal_matrix = matrix<double, diagonal<2>>;
  using dense_matrix = matrix<double, sparse::dense<2, 2>>;

  const diagonal_matrix listed{1., 2.};

  diagonal_matrix copy{3., 4.};
  dense_matrix assigned{5., 6., 7., 8.};

  // An assignment of the same pattern.
  copy = listed;

  assert((copy.at<0, 0>() == 1.));
  assert((copy.at<1, 1>() == 2.));

  // A wider pattern is assigned from a narrower one, with its other elements
  // zero.
  assigned = listed;

  assert((assigned.at<0, 0>() == 1.));
  assert((assigned.at<0, 1>() == 0.));
  assert((assigned.at<1, 0>() == 0.));
  assert((assigned.at<1, 1>() == 2.));

  // The assignment of a result, of the pattern of the result.
  assigned = listed + listed;

  assert((assigned.at<1, 1>() == 4.));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

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
//! @test Verifies the construction of a matrix from a value at every position,
//! a structural zero explicitly.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_pattern;
  using sparse::literals::operator""_o;

  using band_matrix = sparse_matrix<double, sparse::tridiagonal<3>>;
  using column = sparse_matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})>;
  using row = sparse_matrix<double, make_pattern<1, 3>({{0, 1}})>;

  const band_matrix band{1., 2., 0_o, 3., 4., 5., 0_o, 6., 7.};

  assert((band == band_matrix{1., 2., 3., 4., 5., 6., 7.}));

  // A vector, of either orientation.
  const column first_last{1., 0_o, 3.};
  const row middle{0_o, 2., 0_o};

  assert((first_last.at<0, 0>() == 1.));
  assert((first_last.at<1, 0>() == 0.));
  assert((first_last.at<2, 0>() == 3.));
  assert((middle.at<0, 0>() == 0.));
  assert((middle.at<0, 1>() == 2.));

  // A pattern without a stored element.
  const sparse_matrix<double, sparse::zero<1, 2>> none{0_o, 0_o};

  assert((none.at<0, 1>() == 0.));

  // Every spelling of a zero literal.
  const column spelled{1., 0._o, 3.};
  const column decimal{1., 0.0_o, 3.};
  const column fraction{1., .0_o, 3.};

  assert((spelled == first_last));
  assert((decimal == first_last));
  assert((fraction == first_last));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

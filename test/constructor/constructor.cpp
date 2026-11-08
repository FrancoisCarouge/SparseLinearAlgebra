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
//! @test Verifies the constructions of a sparse matrix.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::diagonal;
  using sparse::literals::operator""_o;
  using sparse::tridiagonal;

  using diagonal_matrix = sparse_matrix<double, diagonal<2>>;
  using dense_matrix = sparse_matrix<double, sparse::dense<2, 2>>;

  // A default construction is zero-initialized.
  const diagonal_matrix zero_value;

  assert((zero_value.at<0, 0>() == 0.));
  assert((zero_value.at<1, 1>() == 0.));

  // Every position is given, a structural zero explicitly.
  const diagonal_matrix listed{1., 0_o, 0_o, 2.};

  assert((listed.at<0, 0>() == 1.));
  assert((listed.at<0, 1>() == 0.));
  assert((listed.at<1, 0>() == 0.));
  assert((listed.at<1, 1>() == 2.));

  // A construction from the stored elements, in row-major order.
  const sparse_matrix<double, tridiagonal<3>> stored{1., 2., 3., 4.,
                                                     5., 6., 7.};

  assert((stored.at<0, 0>() == 1.));
  assert((stored.at<0, 1>() == 2.));
  assert((stored.at<1, 0>() == 3.));
  assert((stored.at<1, 1>() == 4.));
  assert((stored.at<1, 2>() == 5.));
  assert((stored.at<2, 1>() == 6.));
  assert((stored.at<2, 2>() == 7.));
  assert((stored.at<0, 2>() == 0.));
  assert((stored.at<2, 0>() == 0.));

  // A copy owns its storage independently of the source.
  diagonal_matrix copy{listed};

  copy.at<0, 0>() = 9.;

  assert((copy.at<0, 0>() == 9.));
  assert((listed.at<0, 0>() == 1.));

  // A wider pattern is constructed from a narrower one, with its other
  // elements zero, explicitly or implicitly.
  const dense_matrix widened{listed};
  const dense_matrix converted = listed;

  assert((widened.at<0, 0>() == 1.));
  assert((widened.at<0, 1>() == 0.));
  assert((widened.at<1, 1>() == 2.));
  assert((converted == widened));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

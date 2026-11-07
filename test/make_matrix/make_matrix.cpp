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
#include <concepts>

namespace fcarouge::test {
namespace {
//! @test Verifies the construction of a matrix from a value at every position,
//! a structural zero explicitly.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::make_matrix;
  using sparse::make_pattern;
  using sparse::matrix;
  using sparse::structural_zero;

  using band_matrix = matrix<double, sparse::tridiagonal<3>>;
  using column = matrix<double, make_pattern<3, 1>({{0, 0}, {2, 0}})>;
  using row = matrix<double, make_pattern<1, 3>({{0, 1}})>;

  const auto band{make_matrix<band_matrix>(1., 2., structural_zero, 3., 4., 5.,
                                           structural_zero, 6., 7.)};

  static_assert(std::same_as<decltype(make_matrix<band_matrix>(
                                 1., 2., structural_zero, 3., 4., 5.,
                                 structural_zero, 6., 7.)),
                             band_matrix>);
  assert((band == band_matrix{1., 2., 3., 4., 5., 6., 7.}));

  // A vector, of either orientation.
  const auto first_last{make_matrix<column>(1., structural_zero, 3.)};
  const auto middle{make_matrix<row>(structural_zero, 2., structural_zero)};

  assert((first_last.at<0, 0>() == 1.));
  assert((first_last.at<1, 0>() == 0.));
  assert((first_last.at<2, 0>() == 3.));
  assert((middle.at<0, 0>() == 0.));
  assert((middle.at<0, 1>() == 2.));

  // A value implicitly convertible to the element type.
  const auto converted{make_matrix<column>(1, structural_zero, 3.F)};

  assert((converted == first_last));

  // A pattern without a structural zero, and one without a stored element.
  const auto dense_value{
      make_matrix<matrix<double, sparse::dense<1, 2>>>(1., 2.)};
  const auto none{make_matrix<matrix<double, sparse::zero<1, 2>>>(
      structural_zero, structural_zero)};

  assert((dense_value.at<0, 1>() == 2.));
  assert((none.at<0, 1>() == 0.));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

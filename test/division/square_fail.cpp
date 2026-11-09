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
//! @test Verifies a division by a rectangular denominator does not compile:
//! the quotient of `x * rhs = lhs` is unique only for a square denominator.
[[maybe_unused]] const auto test{[] -> int {
  const sparse_matrix<double, sparse::dense<2, 3>> lhs{{1., 2., 3.},
                                                       {4., 5., 6.}};
  const sparse_matrix<double, sparse::dense<3, 3>> square{
      {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  const sparse_matrix<double, sparse::dense<2, 3>> rectangular{{1., 0., 0.},
                                                               {0., 1., 0.}};

  // Intended:
  // const auto divided{lhs / square};
  const auto divided{lhs / rectangular};

  assert((divided.at<0, 0>() == 1.));
  assert((square.at<0, 0>() == 1.));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

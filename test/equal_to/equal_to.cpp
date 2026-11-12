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
//! @test Verifies the equality compares values, not patterns.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using sparse::dense;
  using sparse::diagonal;

  const sparse_matrix<double, diagonal<2>> narrow{1., 2.};
  const sparse_matrix<double, dense<2, 2>> wide{1., 0., 0., 2.};
  const sparse_matrix<double, dense<2, 2>> different{1., 3., 0., 2.};
  const sparse_matrix<double, sparse::tridiagonal<2>> band{1., 0., 0., 2.};

  // A structural zero equals a stored zero, in both directions.
  assert((narrow == wide));
  assert((wide == narrow));
  assert((narrow == band));

  // A stored non-zero differs from a structural zero, in both directions.
  assert((not(narrow == different)));
  assert((not(different == narrow)));
  assert((narrow != different));

  // A difference at a stored position, of the same pattern, and of different
  // patterns, in both directions.
  assert((narrow != sparse_matrix<double, diagonal<2>>{1., 3.}));
  assert((not(narrow == sparse_matrix<double, dense<2, 2>>{1., 0., 0., 3.})));
  assert((not(sparse_matrix<double, dense<2, 2>>{1., 0., 0., 3.} == narrow)));

  // The zero matrices of every pattern.
  static_assert(sparse_matrix<double, sparse::zero<2, 2>>{} ==
                sparse_matrix<double, dense<2, 2>>{});

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

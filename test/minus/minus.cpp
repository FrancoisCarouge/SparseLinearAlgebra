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
#include <type_traits>

namespace fcarouge::test {
namespace {
//! @test Verifies the negation keeps the pattern and negates the stored
//! elements.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using first_row =
      sparse_matrix<double, sparse::make_pattern<2, 2>({{0, 0}, {0, 1}})>;

  const first_row value{1., -2.};

  const auto negated{-value};

  static_assert(
      std::same_as<std::remove_cvref_t<decltype(negated)>, first_row>);
  assert((negated.at<0, 0>() == -1.));
  assert((negated.at<0, 1>() == 2.));
  assert((negated.at<1, 1>() == 0.));
  assert((-negated == value));

  // An integer element type.
  const sparse_matrix<int, sparse::diagonal<2>> integers{1, -2};

  assert(((-integers).at<1, 1>() == 2));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

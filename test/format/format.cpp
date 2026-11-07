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

#include "fcarouge/realtime.hpp"
#include "fcarouge/sparse_linear_algebra.hpp"

#include <cassert>
#include <format>

namespace fcarouge::test {
namespace {
//! @test Verifies a structural zero is formatted as an element, so the format
//! does not depend on the pattern.
[[maybe_unused]] const auto test{[] -> int {
  const not_realtime opt_out;

  using sparse::make_pattern;
  using sparse::matrix;

  const matrix<double, sparse::diagonal<2>> value{1., 2.5};
  const matrix<double, sparse::dense<2, 2>> dense_value{1., 0., 0., 2.5};
  const matrix<double, make_pattern<3, 1>({{1, 0}})> column{4.};
  const matrix<double, make_pattern<1, 3>({{0, 2}})> row{5.};
  const matrix<double, sparse::zero<1, 1>> singleton;
  const matrix<int, sparse::dense<1, 1>> stored_singleton{7};

  assert(std::format("{}", value) == "[[1, 0], [0, 2.5]]");
  assert(std::format("{}", value) == std::format("{}", dense_value));
  assert(std::format("{}", column) == "[[0], [4], [0]]");
  assert(std::format("{}", row) == "[0, 0, 5]");
  assert(std::format("{}", singleton) == "0");
  assert(std::format("{}", stored_singleton) == "7");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

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

#include <array>
#include <cassert>
#include <limits>

namespace fcarouge::test {
namespace {
//! @test Verifies the run-time square root of a freestanding implementation,
//! computed by Newton's method rather than by `std::sqrt`: the exact root of a
//! perfect square, below and above one, and the identity on zero and on an
//! infinite value.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: the square root is `constexpr`, and this test
  // would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  constexpr double infinity{std::numeric_limits<double>::infinity()};
  constexpr std::array<std::array<double, 2>, 6> roots{{{0., 0.},
                                                        {1., 1.},
                                                        {0.25, 0.5},
                                                        {144., 12.},
                                                        {1e10, 1e5},
                                                        {infinity, infinity}}};
  [[maybe_unused]] const volatile float quarter{0.25F};

  for ([[maybe_unused]] const auto &[value, root] : roots) {
    assert(sparse_linear_algebra::internal::square_root(value) == root);
  }

  assert(sparse_linear_algebra::internal::square_root(quarter) == 0.5F);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

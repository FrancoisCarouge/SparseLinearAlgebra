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

//! @file
//! @brief A constant velocity state prediction, with a triangular transition.
//!
//! @details The prediction of a one-dimension constant velocity state estimate
//! and of its uncertainty. The state transition is upper triangular, the
//! uncertainty and the process noise are diagonal: the sparsity patterns
//! follow the computation, so the product only computes the terms the
//! structure admits, and the pattern of each result is known at compile time.

#include "fcarouge/realtime.hpp"
#include "fcarouge/sparse_linear_algebra.hpp"

#include <cassert>
#include <concepts>
#include <print>
#include <type_traits>

namespace fcarouge::sample {
namespace {
//! @example constant_velocity.cpp
[[maybe_unused]] const auto sample{[] -> int {
  const not_realtime opt_out;

  using sparse::dense;
  using sparse::diagonal;
  using sparse::upper_triangular;
  using sparse::literals::operator""_o;

  using transition = sparse_matrix<double, upper_triangular<2>>;
  using noise = sparse_matrix<double, diagonal<2>>;
  using state = sparse_matrix<double, dense<2, 1>>;
  using covariance = sparse_matrix<double, dense<2, 2>>;

  constexpr double period{0.5};

  // The bottom left element is a structural zero, without storage.
  static_assert(sizeof(transition) == 3 * sizeof(double));

  const transition model{1., period, 0_o, 1.};
  const state initial{1., 3.};
  const noise uncertainty{4., 0_o, 0_o, 1.};
  const noise process{0.25, 0_o, 0_o, 0.5};

  // The prediction of the state: the position moves by the velocity.
  const auto predicted{model * initial};

  // The prediction of the uncertainty: upper triangular by diagonal by lower
  // triangular, plus diagonal, is dense.
  const auto propagated{model * uncertainty * transposed(model) + process};

  static_assert(std::same_as<std::remove_cvref_t<decltype(predicted)>, state>);
  static_assert(
      std::same_as<std::remove_cvref_t<decltype(propagated)>, covariance>);

  std::println("State: {}", predicted);
  std::println("Uncertainty: {}", propagated);

  assert((predicted == state{2.5, 3.}));
  assert((propagated == covariance{4.5, 0.5, 0.5, 1.5}));

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample

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
//! @brief The heat diffusion along a rod, by a tridiagonal stencil.
//!
//! @details An explicit time integration of the one-dimension heat equation
//! over eight cells, of fixed zero temperature boundaries. Each step multiplies
//! the temperatures by a tridiagonal stencil: a cell keeps half of its heat
//! and receives a quarter of each neighbor's. The stencil stores 22 of its 64
//! elements, and a step costs 22 elementary products instead of 64.

#include "fcarouge/realtime.hpp"
#include "fcarouge/sparse_linear_algebra.hpp"

#include <cassert>
#include <cstddef>
#include <print>

namespace fcarouge::sample {
namespace {
//! @example heat_diffusion.cpp
[[maybe_unused]] const auto sample{[] -> int {
  const not_realtime opt_out;

  using sparse::dense;
  using sparse::sparse_matrix;
  using sparse::tridiagonal;

  constexpr std::size_t cells{8};

  using stencil = sparse_matrix<double, tridiagonal<cells>>;
  using temperatures = sparse_matrix<double, dense<cells, 1>>;

  constexpr std::size_t stored{stencil::structure.count};

  static_assert(stored == 22);
  static_assert(
      sparse::product_term_count<tridiagonal<cells>, dense<cells, 1>> == 22);

  stencil diffusion;
  temperatures rod{0., 0., 0., 16., 16., 0., 0., 0.};

  // Each stored element by its position: a cell or one of its neighbors.
  for (std::size_t i{0}; i < stored; ++i) {
    const auto [row, column]{stencil::structure.positions.at(i)};

    diffusion.values.at(i) = row == column ? 0.5 : 0.25;
  }

  std::println("Step 0: {}", rod);

  for (int step{1}; step <= 3; ++step) {
    rod = diffusion * rod;

    std::println("Step {}: {}", step, rod);
  }

  assert((rod == temperatures{0.25, 1.75, 5.25, 8.75, 8.75, 5.25, 1.75, 0.25}));

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample

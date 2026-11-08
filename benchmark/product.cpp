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

#include <nanobench.h>

#include <cstddef>
#include <format>
#include <fstream>
#include <ios>
#include <random>
#include <string>

namespace fcarouge::benchmark {
namespace {
//! @brief The benchmark measures, prints, and tears down its counters past its
//! entry point: it opts out of the real-time verification for its lifetime.
const not_realtime opt_out;

//! @brief Nanoseconds in a second, nanobench measuring elapsed time in seconds.
constexpr double nanoseconds_per_second{1e9};

//! @brief Percent in a ratio, nanobench reporting its median absolute
//! percentage error as a ratio.
constexpr double percent_per_ratio{100.};

//! @brief Square matrix-matrix product of the given pattern, of random stored
//! elements, appended as a Markdown table row of its median elapsed time in
//! nanoseconds and its median absolute percentage error in percent.
template <auto Size, auto Pattern> void bench(const char *title) {
  sparse::matrix<double, Pattern> lhs;
  sparse::matrix<double, Pattern> rhs;
  std::random_device device;
  std::mt19937 generator{device()};
  std::uniform_real_distribution<> distribution{0., 1.};

  for (std::size_t i{0}; i < Pattern.count; ++i) {
    lhs.values[i] = distribution(generator);
    rhs.values[i] = distribution(generator);
  }

  ankerl::nanobench::Bench runner;
  runner.output(nullptr).title(title).run([&lhs, &rhs]() -> void {
    const auto product{lhs * rhs};
    ankerl::nanobench::doNotOptimizeAway(product);
  });

  const auto &result{runner.results().front()};
  std::ofstream results{"results.txt", std::ios::app};
  results << std::format(
      "| {} | {}x{} | {:.2f} | {:.2f} |\n", title, Size, Size,
      result.median(ankerl::nanobench::Result::Measure::elapsed) *
          nanoseconds_per_second,
      result.medianAbsolutePercentError(
          ankerl::nanobench::Result::Measure::elapsed) *
          percent_per_ratio);
}
} // namespace
} // namespace fcarouge::benchmark

//! @benchmark The square matrix-matrix product of a dense, a tridiagonal, and a
//! diagonal pattern, of one size.
auto main() -> int {
  constexpr std::size_t size{${SIZE}};

  fcarouge::benchmark::bench<size, fcarouge::sparse::dense<size, size>>(
      "dense");
  fcarouge::benchmark::bench<size, fcarouge::sparse::tridiagonal<size>>(
      "tridiagonal");
  fcarouge::benchmark::bench<size, fcarouge::sparse::diagonal<size>>(
      "diagonal");

  return 0;
}

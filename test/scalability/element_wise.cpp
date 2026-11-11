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

#include <concepts>
#include <cstddef>

namespace fcarouge::test {
namespace {
//! @test Verifies the element-wise operations of a banded and a dense matrix of
//! the supported scale in constant evaluation, within the compilers' default
//! constant evaluation limits.
inline constexpr std::size_t size{100};

using band = sparse_matrix<double, sparse::tridiagonal<size>>;
using dense = sparse_matrix<double, sparse::dense<size, size>>;

inline constexpr auto lhs{[] -> band {
  band value;

  for (auto &element : value.values) {
    element = 2.;
  }

  return value;
}()};

inline constexpr auto rhs{[] -> dense {
  dense value;

  for (auto &element : value.values) {
    element = 1.;
  }

  return value;
}()};

inline constexpr auto sum{lhs + rhs};
inline constexpr auto difference{lhs - rhs};
inline constexpr dense converted{lhs};

static_assert(std::same_as<decltype(sum), const dense>);
static_assert(sum[0, 0] == 3.);
static_assert(sum[0, 2] == 1.);
static_assert(sum[99, 98] == 3.);
static_assert(difference[0, 1] == 1.);
static_assert(difference[99, 0] == -1.);
static_assert(converted[99, 99] == 2.);
static_assert(converted[0, 99] == 0.);
static_assert(converted == lhs);
static_assert(lhs == converted);
static_assert(not(lhs == rhs));
static_assert(sum - rhs == lhs);
} // namespace
} // namespace fcarouge::test

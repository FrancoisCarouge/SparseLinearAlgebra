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

#include <bit>
#include <cstdint>
#include <limits>

namespace fcarouge::test {
namespace {
//! @test Verifies the compile-time square root reaches the exact root of a
//! perfect square, and the identity on zero and one.
static_assert(sparse_linear_algebra::internal::square_root(0.) == 0.);
static_assert(sparse_linear_algebra::internal::square_root(1.) == 1.);
static_assert(sparse_linear_algebra::internal::square_root(0.25) == 0.5);
static_assert(sparse_linear_algebra::internal::square_root(144.) == 12.);
static_assert(sparse_linear_algebra::internal::square_root(1e10) == 1e5);
static_assert(sparse_linear_algebra::internal::square_root(0.25F) == 0.5F);

//! @test Verifies the compile-time square root of an infinite or a
//! not-a-number value, as `std::sqrt`'s.
static_assert(sparse_linear_algebra::internal::square_root(
                  std::numeric_limits<double>::infinity()) ==
              std::numeric_limits<double>::infinity());

// A not-a-number value is recognized by its encoding, all exponent bits set
// and a non-zero significand, rather than by comparisons: MSVC does not
// evaluate a comparison with a not-a-number value correctly in a constant
// expression.
inline constexpr std::uint64_t exponent{0x7FF0'0000'0000'0000};
inline constexpr std::uint64_t significand{0x000F'FFFF'FFFF'FFFF};
inline constexpr std::uint64_t not_a_number_root{
    std::bit_cast<std::uint64_t>(sparse_linear_algebra::internal::square_root(
        std::numeric_limits<double>::quiet_NaN()))};

static_assert((not_a_number_root & exponent) == exponent);
static_assert((not_a_number_root & significand) != 0);
} // namespace
} // namespace fcarouge::test

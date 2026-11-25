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

#include <cstddef>

namespace fcarouge::benchmark {
//! @brief The order of the square matrices.
inline constexpr std::size_t size{${SIZE}};

//! @brief The number of sub-diagonals and super-diagonals stored.
inline constexpr std::size_t band{${BAND}};

//! @brief The pattern of the square matrices.
inline constexpr auto structure{sparse::banded<size, size, band, band>};

//! @brief The pattern of the row vectors, the first row of the matrices.
inline constexpr auto row_structure{sparse::banded<1, size, band, band>};

//! @brief The pattern of the column vectors, the first column of the matrices.
inline constexpr auto column_structure{sparse::banded<size, 1, band, band>};

using matrix = sparse_matrix<float, structure>;
using row = sparse_matrix<float, row_structure>;
using column = sparse_matrix<float, column_structure>;
using product_matrix =
    sparse_matrix<float, sparse::pattern_product<structure, structure>>;
using quotient_matrix =
    sparse_matrix<float, sparse::pattern_quotient<structure, structure>>;
using transposed_matrix =
    sparse_matrix<float, sparse::pattern_transpose<structure>>;
} // namespace fcarouge::benchmark

//! @brief The entry points of the operations, each linked alone, named after
//! their operation.
//!
//! @details The operands are the caller's. A matrix result is statically
//! allocated, a scalar result returned. The library algorithms are qualified,
//! not to call the entry point of the same name.
extern "C" {
using fcarouge::benchmark::column;
using fcarouge::benchmark::matrix;
using fcarouge::benchmark::product_matrix;
using fcarouge::benchmark::quotient_matrix;
using fcarouge::benchmark::row;
using fcarouge::benchmark::transposed_matrix;

[[gnu::flatten]] auto addition(const matrix &lhs, const matrix &rhs)
    -> const matrix * {
  static matrix result;
  result = lhs + rhs;
  return &result;
}

[[gnu::flatten]] auto subtraction(const matrix &lhs, const matrix &rhs)
    -> const matrix * {
  static matrix result;
  result = lhs - rhs;
  return &result;
}

[[gnu::flatten]] auto negation(const matrix &value) -> const matrix * {
  static matrix result;
  result = -value;
  return &result;
}

[[gnu::flatten]] auto product(const matrix &lhs, const matrix &rhs)
    -> const product_matrix * {
  static product_matrix result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto scale(const matrix &lhs, float rhs) -> const matrix * {
  static matrix result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto division(const matrix &lhs, const matrix &rhs)
    -> const quotient_matrix * {
  static quotient_matrix result;
  result = lhs / rhs;
  return &result;
}

[[gnu::flatten]] auto scalar_division(const matrix &lhs, float rhs)
    -> const matrix * {
  static matrix result;
  result = lhs / rhs;
  return &result;
}

[[gnu::flatten]] auto transposed(const matrix &value)
    -> const transposed_matrix * {
  static transposed_matrix result;
  result = fcarouge::transposed(value);
  return &result;
}

[[gnu::flatten]] auto dot(const row &lhs, const column &rhs) -> float {
  return fcarouge::dot(lhs, rhs);
}

[[gnu::flatten]] auto equal_to(const matrix &lhs, const matrix &rhs) -> bool {
  return lhs == rhs;
}

[[gnu::flatten]] auto matrix_one_norm(const matrix &value) -> float {
  return fcarouge::matrix_one_norm(value);
}

[[gnu::flatten]] auto matrix_inf_norm(const matrix &value) -> float {
  return fcarouge::matrix_inf_norm(value);
}

[[gnu::flatten]] auto matrix_frob_norm(const matrix &value) -> float {
  return fcarouge::matrix_frob_norm(value);
}
}

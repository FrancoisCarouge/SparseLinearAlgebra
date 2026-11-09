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
using literals::operator""_o;

//! @test Verifies the division by a matrix, the quotient `x` of
//! `x * rhs = lhs`, by a singleton, a diagonal, triangular, block triangular,
//! and irreducible denominators, of square and rectangular numerators.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  // A singleton denominator divides each stored element.
  const sparse_matrix<double, sparse::dense<2, 1>> column{{1.}, {-3.}};
  const sparse_matrix<double, sparse::dense<1, 1>> singleton{{4.}};
  const auto scaled{column / singleton};

  assert((scaled.at<0, 0>() == 0.25));
  assert((scaled.at<1, 0>() == -0.75));

  // A diagonal denominator divides each column, and keeps the pattern.
  using tridiagonal_matrix = sparse_matrix<double, sparse::tridiagonal<3>>;

  const tridiagonal_matrix band{1., 2., 3., 4., 5., 6., 7.};
  const sparse_matrix<double, sparse::diagonal<3>> diagonal{2., 4., 0.5};
  const auto columns{band / diagonal};

  static_assert(
      std::same_as<std::remove_cvref_t<decltype(columns)>, tridiagonal_matrix>);
  assert((columns.at<0, 0>() == 0.5));
  assert((columns.at<0, 1>() == 0.5));
  assert((columns.at<1, 2>() == 10.));
  assert((columns * diagonal == band));

  // An upper triangular denominator substitutes forward along the row.
  constexpr sparse::pattern first_of_row{sparse::shape<1, 3>, {{0, 0}}};

  const sparse_matrix<double, first_of_row> first{{2., 0_o, 0_o}};
  const sparse_matrix<double, sparse::upper_triangular<3>> upper{
      {2., 1., 0.5}, {0_o, 4., 1.}, {0_o, 0_o, 0.5}};
  const auto substituted{first / upper};

  assert((substituted.at<0, 0>() == 1.));
  assert((substituted.at<0, 1>() == -0.25));
  assert((substituted.at<0, 2>() == -0.5));
  assert((substituted * upper == first));

  // A lower triangular denominator divides the identity into its inverse.
  const sparse_matrix<double, sparse::diagonal<2>> identity{1., 1.};
  const sparse_matrix<double, sparse::lower_triangular<2>> lower{{2., 0_o},
                                                                 {1., 4.}};
  const auto inverse{identity / lower};

  static_assert(
      std::same_as<std::remove_cvref_t<decltype(inverse)>,
                   sparse_matrix<double, sparse::lower_triangular<2>>>);
  assert((inverse.at<0, 0>() == 0.5));
  assert((inverse.at<1, 0>() == -0.125));
  assert((inverse.at<1, 1>() == 0.25));
  assert((inverse * lower == identity));

  // A dense denominator whose elimination exchanges its rows to pivot on the
  // largest element, of a rectangular numerator.
  const sparse_matrix<double, sparse::dense<2, 2>> dense{{1., 2.}, {1., 4.}};
  const sparse_matrix<double, sparse::dense<3, 2>> rectangular{
      {1., 0.}, {0., 1.}, {1., 1.}};
  const auto solved{rectangular / dense};

  assert((solved.at<0, 0>() == 2.));
  assert((solved.at<0, 1>() == -1.));
  assert((solved.at<1, 0>() == -0.5));
  assert((solved.at<1, 1>() == 0.5));
  assert((solved.at<2, 0>() == 1.5));
  assert((solved.at<2, 1>() == -0.5));
  assert((solved * dense == rectangular));

  // A block triangular denominator solves its last unknown first, then the
  // irreducible block of the others, updated by it.
  constexpr sparse::pattern coupled{
      sparse::shape<3, 3>, {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {2, 1}, {2, 2}}};

  const sparse_matrix<double, sparse::dense<1, 3>> row{{1., 1., 2.}};
  const sparse_matrix<double, coupled> blocks{
      {1., 2., 0_o}, {1., 4., 0_o}, {0_o, 1., 2.}};
  const auto unknowns{row / blocks};

  assert((unknowns.at<0, 0>() == 2.));
  assert((unknowns.at<0, 1>() == -1.));
  assert((unknowns.at<0, 2>() == 1.));
  assert((unknowns * blocks == row));

  // An irreducible tridiagonal denominator has a dense inverse.
  const sparse_matrix<double, sparse::diagonal<3>> ones{1., 1., 1.};
  const tridiagonal_matrix stencil{1., 1., 1., 2., 1., 1., 2.};
  const auto filled{ones / stencil};

  static_assert(std::same_as<std::remove_cvref_t<decltype(filled)>,
                             sparse_matrix<double, sparse::dense<3, 3>>>);
  assert((filled.at<0, 0>() == 3.));
  assert((filled.at<0, 2>() == 1.));
  assert((filled.at<1, 1>() == 2.));
  assert((filled.at<2, 1>() == -1.));
  assert((filled * stencil == ones));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

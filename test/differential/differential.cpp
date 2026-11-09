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
#include <cmath>
#include <cstddef>
#include <utility>

namespace fcarouge::test {
namespace {
//! @brief The verification of the operations against the textbook dense
//! computations, of explicit loops over arrays.
struct differential {
  template <std::size_t Rows, std::size_t Columns>
  using reference = std::array<std::array<double, Columns>, Rows>;

  //! @brief A matrix of the pattern filled with distinct, signed, integral
  //! values of exact products and sums.
  template <auto Pattern>
  [[nodiscard]] static constexpr auto fill(double seed)
      -> sparse_matrix<double, Pattern> {
    sparse_matrix<double, Pattern> result;

    for (std::size_t i{0}; i < Pattern.count; ++i) {
      const double value{seed + static_cast<double>(i)};

      result.values.at(i) = i % 2 == 0 ? value : -value;
    }

    return result;
  }

  //! @brief The matrix as an explicit dense array, read at each position.
  template <auto Pattern>
  [[nodiscard]] static constexpr auto
  to_reference(const sparse_matrix<double, Pattern> &value)
      -> reference<Pattern.rows, Pattern.columns> {
    reference<Pattern.rows, Pattern.columns> result{};

    for (std::size_t row{0}; row < Pattern.rows; ++row) {
      for (std::size_t column{0}; column < Pattern.columns; ++column) {
        result.at(row).at(column) = value[row, column];
      }
    }

    return result;
  }

  //! @brief The sum, the difference, the transposition, and the negation.
  template <auto Lhs, auto Rhs>
    requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
  static constexpr void verify_additive() {
    const auto lhs{fill<Lhs>(3.)};
    const auto rhs{fill<Rhs>(11.)};
    const auto lhs_values{to_reference(lhs)};
    const auto rhs_values{to_reference(rhs)};

    const auto sum{to_reference(lhs + rhs)};
    const auto difference{to_reference(lhs - rhs)};
    const auto negated{to_reference(-lhs)};
    const auto mirrored{to_reference(transposed(lhs))};

    for (std::size_t row{0}; row < Lhs.rows; ++row) {
      for (std::size_t column{0}; column < Lhs.columns; ++column) {
        [[maybe_unused]] const double left{lhs_values.at(row).at(column)};
        [[maybe_unused]] const double right{rhs_values.at(row).at(column)};

        assert(sum.at(row).at(column) == left + right);
        assert(difference.at(row).at(column) == left - right);
        assert(negated.at(row).at(column) == -left);
        assert(mirrored.at(column).at(row) == left);
      }
    }
  }

  //! @brief The product, and its pattern exactly storing the reached positions.
  template <auto Lhs, auto Rhs>
    requires(Lhs.columns == Rhs.rows)
  static constexpr void verify_product() {
    const auto lhs{fill<Lhs>(3.)};
    const auto rhs{fill<Rhs>(11.)};
    const auto lhs_values{to_reference(lhs)};
    const auto rhs_values{to_reference(rhs)};
    const auto product{lhs * rhs};
    const auto values{to_reference(product)};

    for (std::size_t row{0}; row < Lhs.rows; ++row) {
      for (std::size_t column{0}; column < Rhs.columns; ++column) {
        double expected{0.};
        bool reached{false};

        for (std::size_t inner{0}; inner < Lhs.columns; ++inner) {
          expected +=
              lhs_values.at(row).at(inner) * rhs_values.at(inner).at(column);
          reached =
              reached || (Lhs.stores(row, inner) && Rhs.stores(inner, column));
        }

        assert(values.at(row).at(column) == expected);
        assert(decltype(product)::structure.stores(row, column) == reached);
      }
    }
  }

  //! @brief An invertible denominator of the pattern: its stored diagonal
  //! elements dominate.
  template <auto Pattern>
  [[nodiscard]] static constexpr auto fill_invertible()
      -> sparse_matrix<double, Pattern> {
    sparse_matrix<double, Pattern> result{fill<Pattern>(11.)};

    for (std::size_t i{0}; i < Pattern.count; ++i) {
      if (Pattern.positions.at(i).row == Pattern.positions.at(i).column) {
        result.values.at(i) += 1000.;
      }
    }

    return result;
  }

  //! @brief The quotient `x` of `x * rhs = lhs` by a textbook Gauss-Jordan
  //! elimination with partial pivoting of the dense system `rhsᵀ * xᵀ = lhsᵀ`.
  template <std::size_t Rows, std::size_t Size>
  [[nodiscard]] static auto solve(const reference<Rows, Size> &lhs,
                                  const reference<Size, Size> &rhs)
      -> reference<Rows, Size> {
    reference<Size, Size> system{};
    reference<Size, Rows> solution{};

    for (std::size_t row{0}; row < Size; ++row) {
      for (std::size_t column{0}; column < Size; ++column) {
        system.at(row).at(column) = rhs.at(column).at(row);
      }

      for (std::size_t column{0}; column < Rows; ++column) {
        solution.at(row).at(column) = lhs.at(column).at(row);
      }
    }

    for (std::size_t pivot{0}; pivot < Size; ++pivot) {
      std::size_t largest{pivot};

      for (std::size_t row{pivot + 1}; row < Size; ++row) {
        if (std::abs(system.at(row).at(pivot)) >
            std::abs(system.at(largest).at(pivot))) {
          largest = row;
        }
      }

      std::swap(system.at(pivot), system.at(largest));
      std::swap(solution.at(pivot), solution.at(largest));

      for (std::size_t row{0}; row < Size; ++row) {
        if (row != pivot) {
          const double factor{system.at(row).at(pivot) /
                              system.at(pivot).at(pivot)};

          for (std::size_t column{0}; column < Size; ++column) {
            system.at(row).at(column) -= factor * system.at(pivot).at(column);
          }

          for (std::size_t column{0}; column < Rows; ++column) {
            solution.at(row).at(column) -=
                factor * solution.at(pivot).at(column);
          }
        }
      }
    }

    reference<Rows, Size> result{};

    for (std::size_t row{0}; row < Size; ++row) {
      for (std::size_t column{0}; column < Rows; ++column) {
        result.at(column).at(row) =
            solution.at(row).at(column) / system.at(row).at(row);
      }
    }

    return result;
  }

  //! @brief The division by a matrix, against the dense solution, and its
  //! pattern storing every position the dense solution does not cancel.
  template <auto Lhs, auto Rhs>
    requires(Lhs.columns == Rhs.rows)
  static void verify_quotient() {
    const auto lhs{fill<Lhs>(3.)};
    const auto rhs{fill_invertible<Rhs>()};
    const auto quotient{lhs / rhs};
    const auto values{to_reference(quotient)};
    const auto expected{solve(to_reference(lhs), to_reference(rhs))};

    for (std::size_t row{0}; row < Lhs.rows; ++row) {
      for (std::size_t column{0}; column < Rhs.columns; ++column) {
        [[maybe_unused]] const double value{values.at(row).at(column)};
        [[maybe_unused]] const double reference_value{
            expected.at(row).at(column)};

        assert(std::abs(value - reference_value) <=
               1e-12 * (1. + std::abs(reference_value)));
        assert(decltype(quotient)::structure.stores(row, column) ||
               std::abs(reference_value) <= 1e-12);
      }
    }
  }

  template <auto Lhs, auto... Denominators> static void verify_quotients_of() {
    (verify_quotient<Lhs, Denominators>(), ...);
  }

  //! @brief The quotients of each numerator by every denominator.
  template <auto... Denominators> struct quotients_by {
    template <auto... Numerators> static void verify() {
      (verify_quotients_of<Numerators, Denominators...>(), ...);
    }
  };

  //! @brief A matrix computes the same as a matrix without structural zeros.
  template <auto Lhs, auto Rhs> static constexpr void verify_dense() {
    using dense_matrix =
        sparse_matrix<double, sparse::dense<Lhs.rows, Lhs.columns>>;

    const auto lhs{fill<Lhs>(3.)};
    const auto rhs{fill<Rhs>(11.)};

    assert((lhs + rhs) == (dense_matrix{lhs} + rhs));
    assert((lhs - rhs) == (dense_matrix{lhs} - rhs));
    assert((lhs * rhs) == (dense_matrix{lhs} * rhs));
    assert(transposed(lhs) == transposed(dense_matrix{lhs}));
    assert((lhs * 2.) == (dense_matrix{lhs} * 2.));
  }

  template <auto Lhs, auto... Patterns> static constexpr void verify_with() {
    (verify_additive<Lhs, Patterns>(), ...);
    (verify_product<Lhs, Patterns>(), ...);
    (verify_dense<Lhs, Patterns>(), ...);
  }

  template <auto... Patterns> static constexpr void verify_all_pairs() {
    (verify_with<Patterns, Patterns...>(), ...);
  }
};

inline constexpr auto is_arrow{[](std::size_t row, std::size_t column) -> bool {
  return row == 0 || column == 0 || row == column;
}};
inline constexpr sparse::pattern arrow{sparse::shape<3, 3>, is_arrow};
inline constexpr sparse::pattern anti_diagonal{sparse::shape<3, 3>,
                                               {{0, 2}, {1, 1}, {2, 0}}};
inline constexpr sparse::pattern corner{sparse::shape<3, 3>, {{2, 2}}};
inline constexpr sparse::pattern chained{
    sparse::shape<3, 3>, {{0, 0}, {0, 1}, {1, 1}, {1, 2}, {2, 0}}};

//! @test Verifies every pair of a subset of the patterns at compile time,
//! within the compilers' constant evaluation limits.
static_assert([] -> bool {
  differential::verify_all_pairs<sparse::tridiagonal<3>, sparse::zero<3, 3>,
                                 arrow, corner>();

  return true;
}());

//! @test Verifies every pair of a family of patterns, and rectangular factors,
//! at run time.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  constexpr sparse::pattern wide{sparse::shape<2, 4>,
                                 {{0, 0}, {0, 3}, {1, 1}, {1, 2}}};
  constexpr sparse::pattern tall{sparse::shape<4, 3>,
                                 {{0, 2}, {1, 0}, {2, 1}, {3, 0}, {3, 2}}};

  differential::verify_all_pairs<
      sparse::diagonal<3>, sparse::tridiagonal<3>, sparse::lower_triangular<3>,
      sparse::upper_triangular<3>, sparse::dense<3, 3>, sparse::zero<3, 3>,
      arrow, anti_diagonal, corner>();

  differential::verify_product<wide, tall>();

  differential::quotients_by<
      sparse::diagonal<3>, sparse::tridiagonal<3>, sparse::lower_triangular<3>,
      sparse::upper_triangular<3>, sparse::dense<3, 3>, arrow, anti_diagonal,
      chained>::verify<sparse::diagonal<3>, sparse::tridiagonal<3>,
                       sparse::lower_triangular<3>, sparse::upper_triangular<3>,
                       sparse::dense<3, 3>, sparse::zero<3, 3>, arrow,
                       anti_diagonal, corner, tall>();
  differential::verify_product<wide, sparse::dense<4, 3>>();
  differential::verify_product<sparse::dense<2, 4>, tall>();
  differential::verify_product<wide, sparse::zero<4, 3>>();
  differential::verify_additive<wide, sparse::dense<2, 4>>();
  differential::verify_additive<wide, sparse::zero<2, 4>>();

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

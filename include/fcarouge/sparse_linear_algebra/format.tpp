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

#ifndef FCAROUGE_SPARSE_LINEAR_ALGEBRA_FORMAT_TPP
#define FCAROUGE_SPARSE_LINEAR_ALGEBRA_FORMAT_TPP

//! @file
//! @brief Formatting support for the sparse matrix.
//!
//! @note The standard `<format>` library is not part of the freestanding
//! subset: a freestanding implementation omits this formatting support.

#if __STDC_HOSTED__

#include "fcarouge/sparse_linear_algebra/matrix.hpp"
#include "fcarouge/sparse_linear_algebra_forward.hpp"

#include <cstddef>
#include <format>

//! @brief Specialization of the standard formatter for the sparse matrix.
//!
//! @details Formats every element, structural zeros included, so that the
//! format does not depend on the pattern: a singleton as its element, a row
//! vector as a list of elements, a column vector or a matrix as a list of rows.
template <typename Type, fcarouge::sparse::sparsity_pattern auto Pattern,
          typename Char>
struct std::formatter<fcarouge::sparse::matrix<Type, Pattern>, Char> {
  //! @todo Add a pretty, multi-line, aligned format option?
  //! @todo Add a debug, showing the stored elements, format option?
  static constexpr auto
  parse(std::basic_format_parse_context<Char> &parse_context) {
    return parse_context.begin();
  }

  template <typename FormatContext>
  constexpr auto format(const fcarouge::sparse::matrix<Type, Pattern> &value,
                        FormatContext &format_context) const
      -> FormatContext::iterator {
    if constexpr (Pattern.rows == 1 and Pattern.columns == 1) {
      format_context.advance_to(
          std::format_to(format_context.out(), "{}", value[0, 0]));
    } else if constexpr (Pattern.rows == 1) {
      format_context.advance_to(std::format_to(format_context.out(), "["));

      for (std::size_t column{0}; column < Pattern.columns; ++column) {
        if (column > 0) {
          format_context.advance_to(std::format_to(format_context.out(), ", "));
        }

        format_context.advance_to(
            std::format_to(format_context.out(), "{}", value[0, column]));
      }

      format_context.advance_to(std::format_to(format_context.out(), "]"));
    } else {
      format_context.advance_to(std::format_to(format_context.out(), "["));

      for (std::size_t row{0}; row < Pattern.rows; ++row) {
        if (row > 0) {
          format_context.advance_to(std::format_to(format_context.out(), ", "));
        }

        format_context.advance_to(std::format_to(format_context.out(), "["));

        for (std::size_t column{0}; column < Pattern.columns; ++column) {
          if (column > 0) {
            format_context.advance_to(
                std::format_to(format_context.out(), ", "));
          }

          format_context.advance_to(
              std::format_to(format_context.out(), "{}", value[row, column]));
        }

        format_context.advance_to(std::format_to(format_context.out(), "]"));
      }

      format_context.advance_to(std::format_to(format_context.out(), "]"));
    }

    return format_context.out();
  }
};

#endif // __STDC_HOSTED__
#endif // FCAROUGE_SPARSE_LINEAR_ALGEBRA_FORMAT_TPP

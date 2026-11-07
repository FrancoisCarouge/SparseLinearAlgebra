#[[ Sparse Linear Algebra
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

For more information, please refer to <https://unlicense.org> ]]

if(NOT BUILD_TESTING)
  return()
endif()

separate_arguments(TEST_COMMAND UNIX_COMMAND $ENV{COMMAND})

add_executable(sparse_linear_algebra_amalgamate_sparse_linear_algebra_driver
               "sparse_linear_algebra.cpp")
target_include_directories(
  sparse_linear_algebra_amalgamate_sparse_linear_algebra_driver
  PRIVATE "${CMAKE_BINARY_DIR}/amalgamate/fcarouge")
target_link_libraries(
  sparse_linear_algebra_amalgamate_sparse_linear_algebra_driver
  PRIVATE sparse_linear_algebra_main sparse_linear_algebra_options)
add_dependencies(sparse_linear_algebra_amalgamate_sparse_linear_algebra_driver
                 amalgamate_sparse_linear_algebra)
add_test(
  NAME sparse_linear_algebra_amalgamate_sparse_linear_algebra
  COMMAND
    ${TEST_COMMAND}
    $<TARGET_FILE:sparse_linear_algebra_amalgamate_sparse_linear_algebra_driver>
)

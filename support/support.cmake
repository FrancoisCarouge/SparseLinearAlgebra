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

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # Workaround: CMake does not record MSVC as supporting the C++26 language
  # dialect yet, so requesting it through `CMAKE_CXX_STANDARD` hard-fails
  # configuration, including inside fetched dependencies' own checks. Leave the
  # standard unset on MSVC and rely on `/std:c++latest`, set in
  # `support/CMakeLists.txt`, to opt in to C++26 features directly.
  unset(CMAKE_CXX_STANDARD)
endif()

check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT _)

# Interprocedural optimization defaults on for the benchmark drivers, unless
# explicitly disabled: static analysis captures the compilation of each
# translation unit and has no use for link-time optimized objects.
if(DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION
   AND NOT CMAKE_INTERPROCEDURAL_OPTIMIZATION)
  set(IPO_SUPPORTED OFF)
endif()

# Add a given run-time test.
#
# The test passes when its file compiles, links, and runs to completion: its
# verifications are assertions. The test is named after the backend, the
# directory it is declared in, and its file name:
# `sparse_linear_algebra_[<backend>_]<directory>_<name>_pass`.
#
# * NAME The name of the test file without extension.
# * BACKENDS Optional list of backends to use against the test.
function(pass NAME)
  cmake_parse_arguments(PARSE_ARGV 1 TEST "" "" "BACKENDS")

  get_filename_component(CALLER "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  if(NOT TEST_BACKENDS)
    set(TEST_BACKENDS "none")
  endif()

  foreach(BACKEND IN LISTS TEST_BACKENDS)
    if(BACKEND STREQUAL "none")
      set(TEST_TARGET "sparse_linear_algebra_${CALLER}_${NAME}")
      set(TEST_BACKEND "")
    else()
      set(TEST_TARGET "sparse_linear_algebra_${BACKEND}_${CALLER}_${NAME}")
      set(TEST_BACKEND "sparse_linear_algebra_${BACKEND}")
    endif()

    add_executable(${TEST_TARGET}_pass_driver "${NAME}.cpp")
    target_link_libraries(
      ${TEST_TARGET}_pass_driver
      PRIVATE
        slinalg sparse_linear_algebra_main sparse_linear_algebra_options
        ${TEST_BACKEND})
    separate_arguments(TEST_COMMAND UNIX_COMMAND $ENV{COMMAND})
    add_test(
      NAME ${TEST_TARGET}_pass
      COMMAND ${TEST_COMMAND} $<TARGET_FILE:${TEST_TARGET}_pass_driver>)
  endforeach()
endfunction(pass)

# Add a given compile-time-only test.
#
# The test passes when its file compiles: its verifications are static
# assertions. No executable is linked nor run. The test is named as a run-time
# test is: `sparse_linear_algebra_[<backend>_]<directory>_<name>_pass`.
#
# * NAME The name of the test file without extension.
# * BACKENDS Optional list of backends to use against the test.
function(build NAME)
  cmake_parse_arguments(PARSE_ARGV 1 TEST "" "" "BACKENDS")

  get_filename_component(CALLER "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  if(NOT TEST_BACKENDS)
    set(TEST_BACKENDS "none")
  endif()

  foreach(BACKEND IN LISTS TEST_BACKENDS)
    if(BACKEND STREQUAL "none")
      set(TEST_TARGET "sparse_linear_algebra_${CALLER}_${NAME}")
      set(TEST_BACKEND "")
    else()
      set(TEST_TARGET "sparse_linear_algebra_${BACKEND}_${CALLER}_${NAME}")
      set(TEST_BACKEND "sparse_linear_algebra_${BACKEND}")
    endif()

    add_library(${TEST_TARGET}_build_driver OBJECT "${NAME}.cpp")
    target_link_libraries(
      ${TEST_TARGET}_build_driver
      PRIVATE slinalg sparse_linear_algebra_options ${TEST_BACKEND})
    set_target_properties(
      ${TEST_TARGET}_build_driver
      PROPERTIES EXCLUDE_FROM_ALL TRUE)
    add_test(
      NAME ${TEST_TARGET}_pass
      COMMAND
        ${CMAKE_COMMAND}
        --build
        ${CMAKE_BINARY_DIR}
        --config
        $<CONFIG>
        --target
        ${TEST_TARGET}_build_driver)
    set_tests_properties(
      ${TEST_TARGET}_pass
      PROPERTIES RESOURCE_LOCK "sparse_linear_algebra_build_tree")
  endforeach()
endfunction(build)

# Add a given compile-fail test.
#
# The test passes when its file fails to compile. By convention, the file name
# ends in `_fail`, and the valid form of the rejected line is kept in an
# `// Intended:` comment above it. The test is named after the backend, the
# directory it is declared in, and its file name:
# `sparse_linear_algebra_[<backend>_]<directory>_<name>`.
#
# * NAME The name of the test file without extension.
# * BACKENDS Optional list of backends to use against the test.
function(fail NAME)
  cmake_parse_arguments(PARSE_ARGV 1 TEST "" "" "BACKENDS")

  get_filename_component(CALLER "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  if(NOT TEST_BACKENDS)
    set(TEST_BACKENDS "none")
  endif()

  foreach(BACKEND IN LISTS TEST_BACKENDS)
    if(BACKEND STREQUAL "none")
      set(TEST_TARGET "sparse_linear_algebra_${CALLER}_${NAME}")
      set(TEST_BACKEND "")
    else()
      set(TEST_TARGET "sparse_linear_algebra_${BACKEND}_${CALLER}_${NAME}")
      set(TEST_BACKEND "sparse_linear_algebra_${BACKEND}")
    endif()

    add_library(${TEST_TARGET}_driver OBJECT "${NAME}.cpp")
    target_link_libraries(
      ${TEST_TARGET}_driver
      PRIVATE slinalg sparse_linear_algebra_options ${TEST_BACKEND})
    set_target_properties(
      ${TEST_TARGET}_driver
      PROPERTIES EXCLUDE_FROM_ALL TRUE)
    add_test(
      NAME ${TEST_TARGET}
      COMMAND
        ${CMAKE_COMMAND}
        --build
        ${CMAKE_BINARY_DIR}
        --config
        $<CONFIG>
        --target
        ${TEST_TARGET}_driver)
    set_tests_properties(
      ${TEST_TARGET}
      PROPERTIES RESOURCE_LOCK "sparse_linear_algebra_build_tree" WILL_FAIL TRUE
    )
  endforeach()
endfunction(fail)

# Add a given sample.
#
# The sample is a run-time test named `sparse_linear_algebra_sample_<name>`.
#
# * NAME The name of the sample file without extension.
function(sample NAME)
  add_executable(sparse_linear_algebra_sample_${NAME}_driver "${NAME}.cpp")
  target_link_libraries(
    sparse_linear_algebra_sample_${NAME}_driver
    PRIVATE slinalg sparse_linear_algebra_main sparse_linear_algebra_options)
  separate_arguments(SAMPLE_COMMAND UNIX_COMMAND $ENV{COMMAND})
  add_test(
    NAME sparse_linear_algebra_sample_${NAME}
    COMMAND
      ${SAMPLE_COMMAND}
      $<TARGET_FILE:sparse_linear_algebra_sample_${NAME}_driver>)
endfunction(sample)

# Add a given benchmark.
#
# The benchmark's `${SIZE}` placeholder is replaced by its size, and it is a
# run-time test named `sparse_linear_algebra_<directory>_<name>_<size>_bench`.
#
# * NAME The name of the benchmark file without extension.
# * SIZE The size factor under measurement.
function(bench NAME SIZE)
  get_filename_component(CALLER "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  set(BENCH_TARGET "sparse_linear_algebra_${CALLER}_${NAME}_${SIZE}_bench")

  configure_file("${NAME}.cpp" "${NAME}_${SIZE}.cpp")
  add_executable(${BENCH_TARGET}_driver "${NAME}_${SIZE}.cpp")
  target_link_libraries(
    ${BENCH_TARGET}_driver
    PRIVATE
      slinalg sparse_linear_algebra_options sparse_linear_algebra_realtime
      nanobench::nanobench)
  if(IPO_SUPPORTED)
    set_target_properties(
      ${BENCH_TARGET}_driver
      PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
  separate_arguments(BENCH_COMMAND UNIX_COMMAND $ENV{COMMAND})
  add_test(
    NAME ${BENCH_TARGET}
    COMMAND ${BENCH_COMMAND} $<TARGET_FILE:${BENCH_TARGET}_driver>)
endfunction(bench)

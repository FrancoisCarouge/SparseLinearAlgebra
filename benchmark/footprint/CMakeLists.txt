#[[ Sparse Linear Algebra
Version 0.3.0
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

# Add a given embedded footprint measurement.
#
# Cross-compiled, the test passes when the operation's entry point, linked
# alone, measures the expected sizes in bytes. Otherwise, the configuration is
# only compiled. The configuration `<pattern>_<size>` is configured from the
# `configuration.cpp` template: a square matrix of the given size of the
# `dense`, `tridiagonal`, or `diagonal` band pattern. The test is named after
# the directory it is declared in, its configuration, operation, and
# optimization.
#
# * NAME The configuration: `<pattern>_<size>`.
# * OPERATION The operation of the entry point.
# * OPTIMIZATION The optimization level: `Os` or `O2`.
# * FLASH The expected read-only memory size: `.text`, `.rodata`, `.data`.
# * RAM The expected read-write memory size: `.data`, `.bss`.
# * STACK The expected worst-case stack size.
function(footprint NAME OPERATION OPTIMIZATION)
  cmake_parse_arguments(PARSE_ARGV 3 EXPECTED "" "FLASH;RAM;STACK" "")

  if(NOT NAME MATCHES "^(dense|tridiagonal|diagonal)_([0-9]+)$")
    message(FATAL_ERROR "Unknown footprint configuration: ${NAME}")
  endif()
  set(SIZE ${CMAKE_MATCH_2})
  if(CMAKE_MATCH_1 STREQUAL "dense")
    set(BAND ${SIZE})
  elseif(CMAKE_MATCH_1 STREQUAL "tridiagonal")
    set(BAND 1)
  else()
    set(BAND 0)
  endif()
  configure_file("configuration.cpp" "${NAME}.cpp")

  if(NOT CMAKE_CROSSCOMPILING)
    if(NOT TARGET
       sparse_linear_algebra_freestanding_footprint_${NAME}_build_driver)
      build("${NAME}" BACKENDS "freestanding")
      # The entry points' `[[gnu::flatten]]` shapes the cross-compiled image
      # only: MSVC does not recognize it (C5030).
      target_compile_options(
        sparse_linear_algebra_freestanding_footprint_${NAME}_build_driver
        PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/wd5030>)
    endif()
    return()
  endif()

  set(MEMORY sparse_linear_algebra_footprint_memory_${OPTIMIZATION})
  if(NOT TARGET ${MEMORY})
    add_library(${MEMORY} OBJECT "memory.cpp")
    target_compile_options(
      ${MEMORY}
      PRIVATE "-${OPTIMIZATION}" "-fno-tree-loop-distribute-patterns")
    target_link_libraries(
      ${MEMORY}
      PRIVATE sparse_linear_algebra_footprint_options)
  endif()

  set(CONFIGURATION sparse_linear_algebra_footprint_${NAME}_${OPTIMIZATION})
  if(NOT TARGET ${CONFIGURATION})
    add_library(${CONFIGURATION} OBJECT "${NAME}.cpp")
    target_compile_options(${CONFIGURATION} PRIVATE "-${OPTIMIZATION}")
    target_link_libraries(
      ${CONFIGURATION}
      PRIVATE
        slinalg sparse_linear_algebra_freestanding
        sparse_linear_algebra_footprint_options)
  endif()

  string(TOLOWER "${NAME}_${OPERATION}_${OPTIMIZATION}" MEASUREMENT)
  set(IMAGE sparse_linear_algebra_footprint_${MEASUREMENT})
  set(ROW "${CMAKE_CURRENT_BINARY_DIR}/row/${MEASUREMENT}.html")
  add_executable(${IMAGE} $<TARGET_OBJECTS:${CONFIGURATION}>
                 $<TARGET_OBJECTS:${MEMORY}>)
  set_target_properties(${IMAGE} PROPERTIES LINKER_LANGUAGE "CXX")
  target_link_options(
    ${IMAGE}
    PRIVATE
      "-nostdlib" "LINKER:--gc-sections" "LINKER:--entry=${OPERATION}"
      "LINKER:--defsym=memcpy=footprint_copy"
      "LINKER:--defsym=memset=footprint_fill")
  add_test(
    NAME sparse_linear_algebra_freestanding_footprint_${MEASUREMENT}_pass
    COMMAND
      ${CMAKE_COMMAND}
      "-DNM=${CMAKE_NM}"
      "-DREADELF=${CMAKE_READELF}"
      "-DIMAGE=$<TARGET_FILE:${IMAGE}>"
      "-DENTRY_USAGE=$<TARGET_OBJECTS:${CONFIGURATION}>"
      "-DLEAF_USAGE=$<TARGET_OBJECTS:${MEMORY}>"
      "-DNAME=${NAME}"
      "-DENTRY=${OPERATION}"
      "-DLABEL=${OPERATION}"
      "-DOPTIMIZATION=${OPTIMIZATION}"
      "-DFLASH=${EXPECTED_FLASH}"
      "-DRAM=${EXPECTED_RAM}"
      "-DSTACK=${EXPECTED_STACK}"
      "-DROW=${ROW}"
      -P
      "${CMAKE_CURRENT_SOURCE_DIR}/footprint.cmake")
  set_tests_properties(
    sparse_linear_algebra_freestanding_footprint_${MEASUREMENT}_pass
    PROPERTIES FIXTURES_REQUIRED "sparse_linear_algebra_footprint_page")
  set_property(DIRECTORY APPEND PROPERTY FOOTPRINT_ROWS "${ROW}")
endfunction(footprint)

# Add the documentation page test, assembling the rows of every declared
# measurement after them, passed or failed.
function(footprint_page)
  get_property(ROWS DIRECTORY PROPERTY FOOTPRINT_ROWS)
  add_test(
    NAME sparse_linear_algebra_freestanding_footprint_page_pass
    COMMAND
      ${CMAKE_COMMAND}
      "-DROWS=${ROWS}"
      "-DTEMPLATE=${PROJECT_SOURCE_DIR}/documentation/html/footprint/index.html.in"
      "-DPAGE=${CMAKE_BINARY_DIR}/footprint/index.html"
      "-DCOMPILER=${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
      "-DFLAGS=${CMAKE_CXX_FLAGS}"
      "-DVERSION=${PROJECT_VERSION}"
      -P
      "${CMAKE_CURRENT_SOURCE_DIR}/footprint_page.cmake")
  set_tests_properties(
    sparse_linear_algebra_freestanding_footprint_page_pass
    PROPERTIES FIXTURES_CLEANUP "sparse_linear_algebra_footprint_page")
endfunction(footprint_page)

if(CMAKE_CROSSCOMPILING)
  add_library(sparse_linear_algebra_footprint_options INTERFACE)
  target_compile_options(
    sparse_linear_algebra_footprint_options
    INTERFACE "-ffunction-sections" "-fdata-sections" "-fstack-usage")
  cmake_language(DEFER CALL footprint_page)
else()
  add_library(sparse_linear_algebra_footprint_memory OBJECT "memory.cpp")
  set_target_properties(
    sparse_linear_algebra_footprint_memory
    PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

# Changing a criterion requires its reason in the commit message.
footprint("dense_4" "addition" "O2" FLASH 104 RAM 64 STACK 80)
footprint("dense_4" "addition" "Os" FLASH 100 RAM 64 STACK 80)
footprint("dense_4" "division" "O2" FLASH 924 RAM 64 STACK 240)
footprint("dense_4" "division" "Os" FLASH 960 RAM 64 STACK 232)
footprint("dense_4" "dot" "O2" FLASH 60 RAM 0 STACK 0)
footprint("dense_4" "dot" "Os" FLASH 60 RAM 0 STACK 0)
footprint("dense_4" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("dense_4" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("dense_4" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("dense_4" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("dense_4" "matrix_inf_norm" "O2" FLASH 112 RAM 0 STACK 4)
footprint("dense_4" "matrix_inf_norm" "Os" FLASH 108 RAM 0 STACK 12)
footprint("dense_4" "matrix_one_norm" "O2" FLASH 104 RAM 0 STACK 16)
footprint("dense_4" "matrix_one_norm" "Os" FLASH 128 RAM 0 STACK 24)
footprint("dense_4" "negation" "O2" FLASH 100 RAM 64 STACK 72)
footprint("dense_4" "negation" "Os" FLASH 88 RAM 64 STACK 80)
footprint("dense_4" "product" "O2" FLASH 572 RAM 64 STACK 48)
footprint("dense_4" "product" "Os" FLASH 572 RAM 64 STACK 72)
footprint("dense_4" "scalar_division" "O2" FLASH 100 RAM 64 STACK 72)
footprint("dense_4" "scalar_division" "Os" FLASH 88 RAM 64 STACK 80)
footprint("dense_4" "scale" "O2" FLASH 100 RAM 64 STACK 72)
footprint("dense_4" "scale" "Os" FLASH 88 RAM 64 STACK 80)
footprint("dense_4" "subtraction" "O2" FLASH 104 RAM 64 STACK 80)
footprint("dense_4" "subtraction" "Os" FLASH 100 RAM 64 STACK 80)
footprint("dense_4" "transposed" "O2" FLASH 120 RAM 64 STACK 80)
footprint("dense_4" "transposed" "Os" FLASH 108 RAM 64 STACK 80)
footprint("dense_8" "addition" "O2" FLASH 120 RAM 256 STACK 272)
footprint("dense_8" "addition" "Os" FLASH 106 RAM 256 STACK 280)
footprint("dense_8" "division" "O2" FLASH 2140 RAM 256 STACK 824)
footprint("dense_8" "division" "Os" FLASH 2118 RAM 256 STACK 824)
footprint("dense_8" "dot" "O2" FLASH 108 RAM 0 STACK 0)
footprint("dense_8" "dot" "Os" FLASH 108 RAM 0 STACK 0)
footprint("dense_8" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("dense_8" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("dense_8" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("dense_8" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("dense_8" "matrix_inf_norm" "O2" FLASH 128 RAM 0 STACK 4)
footprint("dense_8" "matrix_inf_norm" "Os" FLASH 124 RAM 0 STACK 12)
footprint("dense_8" "matrix_one_norm" "O2" FLASH 128 RAM 0 STACK 40)
footprint("dense_8" "matrix_one_norm" "Os" FLASH 128 RAM 0 STACK 40)
footprint("dense_8" "negation" "O2" FLASH 88 RAM 256 STACK 264)
footprint("dense_8" "negation" "Os" FLASH 82 RAM 256 STACK 272)
footprint("dense_8" "product" "O2" FLASH 5984 RAM 256 STACK 440)
footprint("dense_8" "product" "Os" FLASH 5044 RAM 256 STACK 440)
footprint("dense_8" "scalar_division" "O2" FLASH 100 RAM 256 STACK 272)
footprint("dense_8" "scalar_division" "Os" FLASH 94 RAM 256 STACK 280)
footprint("dense_8" "scale" "O2" FLASH 100 RAM 256 STACK 272)
footprint("dense_8" "scale" "Os" FLASH 94 RAM 256 STACK 280)
footprint("dense_8" "subtraction" "O2" FLASH 120 RAM 256 STACK 272)
footprint("dense_8" "subtraction" "Os" FLASH 106 RAM 256 STACK 280)
footprint("dense_8" "transposed" "O2" FLASH 136 RAM 256 STACK 264)
footprint("dense_8" "transposed" "Os" FLASH 114 RAM 256 STACK 280)
footprint("dense_16" "addition" "O2" FLASH 124 RAM 1024 STACK 1040)
footprint("dense_16" "addition" "Os" FLASH 110 RAM 1024 STACK 1048)
footprint("dense_16" "division" "O2" FLASH 6740 RAM 1024 STACK 3112)
footprint("dense_16" "division" "Os" FLASH 6722 RAM 1024 STACK 3120)
footprint("dense_16" "dot" "O2" FLASH 204 RAM 0 STACK 0)
footprint("dense_16" "dot" "Os" FLASH 204 RAM 0 STACK 0)
footprint("dense_16" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("dense_16" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("dense_16" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("dense_16" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("dense_16" "matrix_inf_norm" "O2" FLASH 164 RAM 0 STACK 4)
footprint("dense_16" "matrix_inf_norm" "Os" FLASH 156 RAM 0 STACK 12)
footprint("dense_16" "matrix_one_norm" "O2" FLASH 132 RAM 0 STACK 72)
footprint("dense_16" "matrix_one_norm" "Os" FLASH 132 RAM 0 STACK 72)
footprint("dense_16" "negation" "O2" FLASH 96 RAM 1024 STACK 1032)
footprint("dense_16" "negation" "Os" FLASH 86 RAM 1024 STACK 1040)
footprint("dense_16" "product" "O2" FLASH 53396 RAM 1024 STACK 1936)
footprint("dense_16" "product" "Os" FLASH 40324 RAM 1024 STACK 2016)
footprint("dense_16" "scalar_division" "O2" FLASH 108 RAM 1024 STACK 1040)
footprint("dense_16" "scalar_division" "Os" FLASH 98 RAM 1024 STACK 1048)
footprint("dense_16" "scale" "O2" FLASH 108 RAM 1024 STACK 1040)
footprint("dense_16" "scale" "Os" FLASH 98 RAM 1024 STACK 1048)
footprint("dense_16" "subtraction" "O2" FLASH 124 RAM 1024 STACK 1040)
footprint("dense_16" "subtraction" "Os" FLASH 110 RAM 1024 STACK 1048)
footprint("dense_16" "transposed" "O2" FLASH 140 RAM 1024 STACK 1032)
footprint("dense_16" "transposed" "Os" FLASH 118 RAM 1024 STACK 1048)
footprint("tridiagonal_4" "addition" "O2" FLASH 100 RAM 40 STACK 56)
footprint("tridiagonal_4" "addition" "Os" FLASH 96 RAM 40 STACK 56)
footprint("tridiagonal_4" "division" "O2" FLASH 836 RAM 64 STACK 232)
footprint("tridiagonal_4" "division" "Os" FLASH 864 RAM 64 STACK 232)
footprint("tridiagonal_4" "dot" "O2" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_4" "dot" "Os" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_4" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_4" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("tridiagonal_4" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("tridiagonal_4" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("tridiagonal_4" "matrix_inf_norm" "O2" FLASH 112 RAM 0 STACK 4)
footprint("tridiagonal_4" "matrix_inf_norm" "Os" FLASH 108 RAM 0 STACK 12)
footprint("tridiagonal_4" "matrix_one_norm" "O2" FLASH 208 RAM 0 STACK 16)
footprint("tridiagonal_4" "matrix_one_norm" "Os" FLASH 220 RAM 0 STACK 24)
footprint("tridiagonal_4" "negation" "O2" FLASH 88 RAM 40 STACK 48)
footprint("tridiagonal_4" "negation" "Os" FLASH 80 RAM 40 STACK 56)
footprint("tridiagonal_4" "product" "O2" FLASH 320 RAM 56 STACK 16)
footprint("tridiagonal_4" "product" "Os" FLASH 316 RAM 56 STACK 16)
footprint("tridiagonal_4" "scalar_division" "O2" FLASH 88 RAM 40 STACK 48)
footprint("tridiagonal_4" "scalar_division" "Os" FLASH 80 RAM 40 STACK 56)
footprint("tridiagonal_4" "scale" "O2" FLASH 88 RAM 40 STACK 48)
footprint("tridiagonal_4" "scale" "Os" FLASH 80 RAM 40 STACK 56)
footprint("tridiagonal_4" "subtraction" "O2" FLASH 100 RAM 40 STACK 56)
footprint("tridiagonal_4" "subtraction" "Os" FLASH 96 RAM 40 STACK 56)
footprint("tridiagonal_4" "transposed" "O2" FLASH 144 RAM 40 STACK 48)
footprint("tridiagonal_4" "transposed" "Os" FLASH 132 RAM 40 STACK 56)
footprint("tridiagonal_8" "addition" "O2" FLASH 116 RAM 88 STACK 104)
footprint("tridiagonal_8" "addition" "Os" FLASH 102 RAM 88 STACK 112)
footprint("tridiagonal_8" "division" "O2" FLASH 1460 RAM 256 STACK 816)
footprint("tridiagonal_8" "division" "Os" FLASH 1446 RAM 256 STACK 824)
footprint("tridiagonal_8" "dot" "O2" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_8" "dot" "Os" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_8" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_8" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("tridiagonal_8" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("tridiagonal_8" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("tridiagonal_8" "matrix_inf_norm" "O2" FLASH 128 RAM 0 STACK 4)
footprint("tridiagonal_8" "matrix_inf_norm" "Os" FLASH 124 RAM 0 STACK 12)
footprint("tridiagonal_8" "matrix_one_norm" "O2" FLASH 332 RAM 0 STACK 40)
footprint("tridiagonal_8" "matrix_one_norm" "Os" FLASH 320 RAM 0 STACK 40)
footprint("tridiagonal_8" "negation" "O2" FLASH 84 RAM 88 STACK 96)
footprint("tridiagonal_8" "negation" "Os" FLASH 78 RAM 88 STACK 104)
footprint("tridiagonal_8" "product" "O2" FLASH 812 RAM 136 STACK 72)
footprint("tridiagonal_8" "product" "Os" FLASH 784 RAM 136 STACK 96)
footprint("tridiagonal_8" "scalar_division" "O2" FLASH 96 RAM 88 STACK 104)
footprint("tridiagonal_8" "scalar_division" "Os" FLASH 90 RAM 88 STACK 112)
footprint("tridiagonal_8" "scale" "O2" FLASH 96 RAM 88 STACK 104)
footprint("tridiagonal_8" "scale" "Os" FLASH 90 RAM 88 STACK 112)
footprint("tridiagonal_8" "subtraction" "O2" FLASH 116 RAM 88 STACK 104)
footprint("tridiagonal_8" "subtraction" "Os" FLASH 102 RAM 88 STACK 112)
footprint("tridiagonal_8" "transposed" "O2" FLASH 204 RAM 88 STACK 96)
footprint("tridiagonal_8" "transposed" "Os" FLASH 186 RAM 88 STACK 104)
footprint("tridiagonal_16" "addition" "O2" FLASH 116 RAM 184 STACK 200)
footprint("tridiagonal_16" "addition" "Os" FLASH 102 RAM 184 STACK 208)
footprint("tridiagonal_16" "division" "O2" FLASH 3392 RAM 1024 STACK 3120)
footprint("tridiagonal_16" "division" "Os" FLASH 3358 RAM 1024 STACK 3120)
footprint("tridiagonal_16" "dot" "O2" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_16" "dot" "Os" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_16" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("tridiagonal_16" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("tridiagonal_16" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("tridiagonal_16" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("tridiagonal_16" "matrix_inf_norm" "O2" FLASH 164 RAM 0 STACK 4)
footprint("tridiagonal_16" "matrix_inf_norm" "Os" FLASH 156 RAM 0 STACK 12)
footprint("tridiagonal_16" "matrix_one_norm" "O2" FLASH 524 RAM 0 STACK 72)
footprint("tridiagonal_16" "matrix_one_norm" "Os" FLASH 512 RAM 0 STACK 72)
footprint("tridiagonal_16" "negation" "O2" FLASH 84 RAM 184 STACK 192)
footprint("tridiagonal_16" "negation" "Os" FLASH 78 RAM 184 STACK 200)
footprint("tridiagonal_16" "product" "O2" FLASH 2280 RAM 296 STACK 240)
footprint("tridiagonal_16" "product" "Os" FLASH 1900 RAM 296 STACK 256)
footprint("tridiagonal_16" "scalar_division" "O2" FLASH 96 RAM 184 STACK 200)
footprint("tridiagonal_16" "scalar_division" "Os" FLASH 90 RAM 184 STACK 208)
footprint("tridiagonal_16" "scale" "O2" FLASH 96 RAM 184 STACK 200)
footprint("tridiagonal_16" "scale" "Os" FLASH 90 RAM 184 STACK 208)
footprint("tridiagonal_16" "subtraction" "O2" FLASH 116 RAM 184 STACK 200)
footprint("tridiagonal_16" "subtraction" "Os" FLASH 102 RAM 184 STACK 208)
footprint("tridiagonal_16" "transposed" "O2" FLASH 300 RAM 184 STACK 192)
footprint("tridiagonal_16" "transposed" "Os" FLASH 282 RAM 184 STACK 200)
footprint("diagonal_4" "addition" "O2" FLASH 80 RAM 16 STACK 56)
footprint("diagonal_4" "addition" "Os" FLASH 96 RAM 16 STACK 64)
footprint("diagonal_4" "division" "O2" FLASH 72 RAM 16 STACK 0)
footprint("diagonal_4" "division" "Os" FLASH 76 RAM 16 STACK 0)
footprint("diagonal_4" "dot" "O2" FLASH 24 RAM 0 STACK 0)
footprint("diagonal_4" "dot" "Os" FLASH 24 RAM 0 STACK 0)
footprint("diagonal_4" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("diagonal_4" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("diagonal_4" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("diagonal_4" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("diagonal_4" "matrix_inf_norm" "O2" FLASH 112 RAM 0 STACK 4)
footprint("diagonal_4" "matrix_inf_norm" "Os" FLASH 108 RAM 0 STACK 12)
footprint("diagonal_4" "matrix_one_norm" "O2" FLASH 160 RAM 0 STACK 16)
footprint("diagonal_4" "matrix_one_norm" "Os" FLASH 172 RAM 0 STACK 24)
footprint("diagonal_4" "negation" "O2" FLASH 56 RAM 16 STACK 40)
footprint("diagonal_4" "negation" "Os" FLASH 60 RAM 16 STACK 48)
footprint("diagonal_4" "product" "O2" FLASH 92 RAM 16 STACK 0)
footprint("diagonal_4" "product" "Os" FLASH 92 RAM 16 STACK 0)
footprint("diagonal_4" "scalar_division" "O2" FLASH 56 RAM 16 STACK 40)
footprint("diagonal_4" "scalar_division" "Os" FLASH 60 RAM 16 STACK 48)
footprint("diagonal_4" "scale" "O2" FLASH 56 RAM 16 STACK 40)
footprint("diagonal_4" "scale" "Os" FLASH 60 RAM 16 STACK 48)
footprint("diagonal_4" "subtraction" "O2" FLASH 80 RAM 16 STACK 56)
footprint("diagonal_4" "subtraction" "Os" FLASH 96 RAM 16 STACK 64)
footprint("diagonal_4" "transposed" "O2" FLASH 92 RAM 16 STACK 40)
footprint("diagonal_4" "transposed" "Os" FLASH 104 RAM 16 STACK 48)
footprint("diagonal_8" "addition" "O2" FLASH 96 RAM 32 STACK 48)
footprint("diagonal_8" "addition" "Os" FLASH 92 RAM 32 STACK 48)
footprint("diagonal_8" "division" "O2" FLASH 144 RAM 32 STACK 8)
footprint("diagonal_8" "division" "Os" FLASH 140 RAM 32 STACK 0)
footprint("diagonal_8" "dot" "O2" FLASH 24 RAM 0 STACK 0)
footprint("diagonal_8" "dot" "Os" FLASH 24 RAM 0 STACK 0)
footprint("diagonal_8" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("diagonal_8" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("diagonal_8" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("diagonal_8" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("diagonal_8" "matrix_inf_norm" "O2" FLASH 128 RAM 0 STACK 4)
footprint("diagonal_8" "matrix_inf_norm" "Os" FLASH 124 RAM 0 STACK 12)
footprint("diagonal_8" "matrix_one_norm" "O2" FLASH 220 RAM 0 STACK 40)
footprint("diagonal_8" "matrix_one_norm" "Os" FLASH 208 RAM 0 STACK 40)
footprint("diagonal_8" "negation" "O2" FLASH 76 RAM 32 STACK 40)
footprint("diagonal_8" "negation" "Os" FLASH 72 RAM 32 STACK 48)
footprint("diagonal_8" "product" "O2" FLASH 180 RAM 32 STACK 8)
footprint("diagonal_8" "product" "Os" FLASH 172 RAM 32 STACK 0)
footprint("diagonal_8" "scalar_division" "O2" FLASH 76 RAM 32 STACK 40)
footprint("diagonal_8" "scalar_division" "Os" FLASH 72 RAM 32 STACK 48)
footprint("diagonal_8" "scale" "O2" FLASH 76 RAM 32 STACK 40)
footprint("diagonal_8" "scale" "Os" FLASH 72 RAM 32 STACK 48)
footprint("diagonal_8" "subtraction" "O2" FLASH 96 RAM 32 STACK 48)
footprint("diagonal_8" "subtraction" "Os" FLASH 92 RAM 32 STACK 48)
footprint("diagonal_8" "transposed" "O2" FLASH 128 RAM 32 STACK 40)
footprint("diagonal_8" "transposed" "Os" FLASH 120 RAM 32 STACK 48)
footprint("diagonal_16" "addition" "O2" FLASH 104 RAM 64 STACK 80)
footprint("diagonal_16" "addition" "Os" FLASH 100 RAM 64 STACK 80)
footprint("diagonal_16" "division" "O2" FLASH 276 RAM 64 STACK 8)
footprint("diagonal_16" "division" "Os" FLASH 276 RAM 64 STACK 8)
footprint("diagonal_16" "dot" "O2" FLASH 24 RAM 0 STACK 0)
footprint("diagonal_16" "dot" "Os" FLASH 24 RAM 0 STACK 0)
footprint("diagonal_16" "equal_to" "O2" FLASH 36 RAM 0 STACK 0)
footprint("diagonal_16" "equal_to" "Os" FLASH 34 RAM 0 STACK 0)
footprint("diagonal_16" "matrix_frob_norm" "O2" FLASH 140 RAM 0 STACK 0)
footprint("diagonal_16" "matrix_frob_norm" "Os" FLASH 136 RAM 0 STACK 0)
footprint("diagonal_16" "matrix_inf_norm" "O2" FLASH 164 RAM 0 STACK 4)
footprint("diagonal_16" "matrix_inf_norm" "Os" FLASH 156 RAM 0 STACK 12)
footprint("diagonal_16" "matrix_one_norm" "O2" FLASH 284 RAM 0 STACK 72)
footprint("diagonal_16" "matrix_one_norm" "Os" FLASH 272 RAM 0 STACK 72)
footprint("diagonal_16" "negation" "O2" FLASH 100 RAM 64 STACK 72)
footprint("diagonal_16" "negation" "Os" FLASH 88 RAM 64 STACK 80)
footprint("diagonal_16" "product" "O2" FLASH 344 RAM 64 STACK 8)
footprint("diagonal_16" "product" "Os" FLASH 340 RAM 64 STACK 8)
footprint("diagonal_16" "scalar_division" "O2" FLASH 100 RAM 64 STACK 72)
footprint("diagonal_16" "scalar_division" "Os" FLASH 88 RAM 64 STACK 80)
footprint("diagonal_16" "scale" "O2" FLASH 100 RAM 64 STACK 72)
footprint("diagonal_16" "scale" "Os" FLASH 88 RAM 64 STACK 80)
footprint("diagonal_16" "subtraction" "O2" FLASH 104 RAM 64 STACK 80)
footprint("diagonal_16" "subtraction" "Os" FLASH 100 RAM 64 STACK 80)
footprint("diagonal_16" "transposed" "O2" FLASH 172 RAM 64 STACK 72)
footprint("diagonal_16" "transposed" "Os" FLASH 160 RAM 64 STACK 80)

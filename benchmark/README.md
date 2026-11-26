# Benchmarks

## Embedded Footprint

The [embedded footprint](https://francoiscarouge.github.io/SparseLinearAlgebra/footprint.xhtml) measures the flash, RAM, and worst-case stack sizes of each operation, for the [configurations](https://github.com/FrancoisCarouge/SparseLinearAlgebra/tree/master/benchmark/footprint) of a dense, a tridiagonal, and a diagonal pattern of single precision floating point elements, per optimization level, on a bare-metal Arm Cortex-M4 target.

Measure on Linux with the Arm GNU Toolchain:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/SparseLinearAlgebra" "sparse"
cmake -S "sparse" -B "build" --toolchain "cmake/toolchain_arm_cortex_m4.cmake"
cmake --build "build" --parallel
ctest --test-dir "build" --output-on-failure --parallel
```

The tests write the measurements to `build/footprint/index.html`, and fail on any size change, larger or smaller, from its `footprint()` line in [footprint/CMakeLists.txt](https://github.com/FrancoisCarouge/SparseLinearAlgebra/blob/master/benchmark/footprint/CMakeLists.txt). The criteria hold for the toolchain of the continuous integration runner, `gcc-arm-none-eabi` 14.2 on Ubuntu 26.04. To accept a change, update the line, with its reason in the commit message.

## Speed

Build and run the benchmarks on Linux platform, building the benchmark drivers only:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/SparseLinearAlgebra" "sparse"
rm -rf "build"
CXX=clang++-21 cmake -S "sparse" -B "build" -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKING=ON
cmake --build "build" --parallel --target $(ninja -C "build" -t targets all | grep -o '^sparse_linear_algebra_benchmark_[a-z_0-9]*_bench_driver' | sort -u)
ctest --test-dir "build" --tests-regex "bench" --parallel 1
```

The results are appended to `build/benchmark/results.txt`.

# Results

Disclaimer: naive benchmark results for illustration purposes only.

The square matrix-matrix product of a dense, a tridiagonal, and a diagonal pattern, of random double precision floating point elements. The cost follows the count of elementary products of each pattern, not the shape: cubic for the dense pattern, linear for the tridiagonal and the diagonal patterns.

Clang 21.1.8, libc++, x86-64.

| Title | Size | Median Elapsed Time (ns) | Median Absolute Percentage Error (%) |
| --- | --- | --- | --- |
| dense | 2x2 | 0.88 | 7.96 |
| tridiagonal | 2x2 | 0.84 | 4.63 |
| diagonal | 2x2 | 0.20 | 1.38 |
| dense | 4x4 | 39.00 | 5.70 |
| tridiagonal | 4x4 | 16.17 | 2.17 |
| diagonal | 4x4 | 0.34 | 6.18 |
| dense | 8x8 | 237.69 | 2.74 |
| tridiagonal | 8x8 | 39.56 | 2.91 |
| diagonal | 8x8 | 0.60 | 2.66 |
| dense | 16x16 | 2061.06 | 4.22 |
| tridiagonal | 16x16 | 72.84 | 2.66 |
| diagonal | 16x16 | 1.21 | 2.44 |

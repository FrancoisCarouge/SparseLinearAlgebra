# Benchmarks

Build and run the benchmarks on Linux platform, building the benchmark drivers only: a release build defines `NDEBUG`, under which the tests' assertions vanish.

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

| Title | Size | Median Elapsed Time (s) | Median Absolute Error (%) |
| --- | --- | --- | --- |
| dense |     2x2     | 8.84068119194387e-10 | 0.0796378221605947 |
| tridiagonal |     2x2     | 8.37212e-10 | 0.0463061723057106 |
| diagonal |     2x2     | 1.96428904515864e-10 | 0.0137714697195507 |
| dense |     4x4     | 3.90028651014723e-08 | 0.0569639195063688 |
| tridiagonal |     4x4     | 1.61721394451299e-08 | 0.021749930594574 |
| diagonal |     4x4     | 3.38985522296087e-10 | 0.0618330477142346 |
| dense |     8x8     | 2.37689876425856e-07 | 0.0273598610942169 |
| tridiagonal |     8x8     | 3.95618884991693e-08 | 0.0291214483268116 |
| diagonal |     8x8     | 6.02481303955703e-10 | 0.0265591817204245 |
| dense |    16x16    | 2.06106387665198e-06 | 0.0421850606649183 |
| tridiagonal |    16x16    | 7.28395095367847e-08 | 0.0265965628835399 |
| diagonal |    16x16    | 1.20535200077102e-09 | 0.0243913759322464 |

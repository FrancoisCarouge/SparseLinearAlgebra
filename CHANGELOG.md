# Changelog

All notable changes to SparseLinearAlgebra are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-11-07

### Added

- Compile-time sparsity patterns: `position`, `pattern`, the `sparsity_pattern` concept, the construction of a pattern from its `shape` and its listed positions or a predicate, and the `dense`, `zero`, `banded`, `diagonal`, `tridiagonal`, `lower_triangular`, and `upper_triangular` named patterns.
- Pattern algebra: `pattern_union`, `pattern_intersection`, `pattern_transpose`, `pattern_product`, `product_term_count`, and `is_subset_of`.
- The `sparse_matrix` storing only the elements of its pattern, with its compile-time `at<Row, Column>()` and run-time `operator[]` element accesses, and its construction from every row with the `0_o` structural zero literal.
- Operations: `+`, `-`, unary `-`, `*` of matrices and by a scalar, `/` by a scalar, `==`, `transposed`, `dot`, `matrix_frob_norm`, `matrix_inf_norm`, and `matrix_one_norm`, each of a result pattern and arithmetic terms computed at compile time.
- `std::formatter` support, on a hosted implementation.
- Freestanding, real-time, and constant expression support: no memory allocation, no exception, no input-output.
- Samples, benchmarks, single-header amalgamation, CMake and pkg-config packages, and AI agent guidance (`AGENTS.md`).

## François Carouge / Sparse Linear Algebra

Most matrices of estimation, control, and simulation problems are mostly zeros, and their structure is known before the program runs: a diagonal noise, a triangular transition, a banded stencil. Ordinary dense linear algebra stores and multiplies those zeros anyway, and ordinary sparse linear algebra discovers them at run time, with memory allocations and indirections. This library is a C++ compile-time sparsity pattern linear algebra: the positions of the stored elements are part of the matrix type. Every other element is a structural zero, occupying no memory and costing no operation. Each operation computes the pattern of its result, and the exact arithmetic terms it needs, at compile time. There is no memory allocation, no exception, no input-output, and every operation is `constexpr`, freestanding, and real-time.

```cpp
using band = matrix<double, tridiagonal<3>>;

band value{1., 2., 3., 4., 5., 6., 7.}; // Stores 7 of the 9 elements.

value.at<1, 2>() = 8.;     // Writes a stored element.
// value.at<0, 2>() = 8.;  // Does not compile: a structural zero.

std::println("{}", value * value); // Computes 17 of the 27 products.

// [[7, 10, 16], [15, 70, 88], [18, 66, 97]]
```

[![Pipeline](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/pipeline.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/pipeline.yml)
[![Freestanding](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/freestanding.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/freestanding.yml)
[![Sanitizer](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/sanitizer.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/sanitizer.yml)
[![Format](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/format.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/format.yml)
[![ClangTidy](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/clang_tidy.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/clang_tidy.yml)
[![CppCheck](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/cppcheck.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/cppcheck.yml)
[![Coverity Scan](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/coverity.yml/badge.svg)](https://scan.coverity.com/projects/francoiscarouge-sparselinearalgebra)
[![Doxygen](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/doxygen.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/doxygen.yml)
[![Valgrind](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/valgrind_memory.yml/badge.svg)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/valgrind_memory.yml)
[![Public Domain](https://img.shields.io/badge/License-Public%20Domain%20%F0%9F%94%97-brightgreen)](https://raw.githubusercontent.com/francoiscarouge/SparseLinearAlgebra/master/LICENSE.txt)
[![Deploy Unit Test Code Coverage](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/deploy_coverage.yml/badge.svg)](https://francoiscarouge.github.io/SparseLinearAlgebra/unit_test_coverage.xhtml)
[![Deploy Doxygen](https://github.com/FrancoisCarouge/SparseLinearAlgebra/actions/workflows/deploy_doxygen.yml/badge.svg)](https://francoiscarouge.github.io/SparseLinearAlgebra/index.xhtml)
[![Sponsor](https://img.shields.io/badge/Support-Sponsor%20%F0%9F%94%97-brightgreen)](http://paypal.me/francoiscarouge)
[![Discussions](https://img.shields.io/github/discussions/FrancoisCarouge/SparseLinearAlgebra)](https://github.com/FrancoisCarouge/SparseLinearAlgebra/discussions)

## Installation & Usage

See [INSTALL.md](https://github.com/FrancoisCarouge/SparseLinearAlgebra/blob/master/INSTALL.md) for installation instructions.

## Reference

### Layers

![Where Sparse Linear Algebra sits in the technical stack: below the application and its optional typed or estimation libraries, building on the C++ language alone and storing arithmetic element types.](https://francoiscarouge.github.io/SparseLinearAlgebra/sparse_linear_algebra/layers.svg)

The library is the compile-time sparsity layer between the application and the C++ language, without dependency.

### CMake

Link against the CMake target `fcarouge-sparse-linear-algebra::slinalg`.

### Include

Include the library header in your sources:

```cpp
#include "fcarouge/sparse_linear_algebra.hpp"
```

For forward declarations only, include instead:

```cpp
#include "fcarouge/sparse_linear_algebra_forward.hpp"
```

The library is header-only and has no dependency. A single-header amalgamation, `sparse_linear_algebra.h`, is also [published](https://francoiscarouge.github.io/SparseLinearAlgebra/amalgamate/fcarouge/sparse_linear_algebra.h).

### Namespace

The public API lives in the `fcarouge::sparse` namespace. Import it with a using-declaration or namespace alias as your project's style prefers.

### Samples

| Use Case | Sample |
| --- | --- |
| Constant velocity prediction | [Usage example: `sample/constant_velocity.cpp`](https://github.com/FrancoisCarouge/SparseLinearAlgebra/blob/master/sample/constant_velocity.cpp). |
| Heat diffusion stencil | [Usage example: `sample/heat_diffusion.cpp`](https://github.com/FrancoisCarouge/SparseLinearAlgebra/blob/master/sample/heat_diffusion.cpp). |

### Class Matrix

Compile-time sparsity pattern matrix. Stores the elements of the pattern in place, in row-major order. Every other element is a structural zero: it reads as a zero, and a write to it does not compile.

Also documented in the [fcarouge/sparse_linear_algebra/matrix.hpp](https://github.com/FrancoisCarouge/SparseLinearAlgebra/blob/master/include/fcarouge/sparse_linear_algebra/matrix.hpp) header.

#### Declaration

```cpp
namespace fcarouge::sparse {
template <typename Type, sparsity_pattern auto Pattern>
struct matrix;
}
```

#### Template Parameters

| Template Parameter | Definition |
| --- | --- |
| `Type` | The type of the stored elements. |
| `Pattern` | The sparsity pattern of the matrix. |

#### Member Variables

| Member Variable | Definition |
| --- | --- |
| `structure` | Static. The sparsity pattern of the matrix. |
| `values` | The `std::array` of the stored elements, zero-initialized, in row-major order. |

#### Member Functions

| Member Function | Definition |
| --- | --- |
| `(default constructor)` | Construct a zero matrix. |
| `(explicit constructor)` | Construct from every stored element, in row-major order. Each one is exactly of `Type`, and there is at least one. |
| `(conversion constructor)` | Implicitly convert from a matrix of the same element type and shape, whose stored positions this pattern all stores. |
| `at<Row, Column>()` | Access the stored element at the compile-time position, by reference. Only for a stored element. |
| `at<Row, Column>() const` | Read the element at the compile-time position, by value. A structural zero reads as `Type{}`. |
| `operator[row, column] const` | Read the element at the run-time position, by value. A structural zero reads as `Type{}`. |
| `operator[row, column]` | Access the element at the run-time position, by reference. Only without a structural zero. |
| `operator[std::integral_constant<std::size_t, Row>, std::integral_constant<std::size_t, Column>]` | Access the stored element at the compile-time position, by reference. Only for a stored element. The generic access of a linear algebra library by compile-time position. |

### Class Pattern

The positions of the stored elements of a matrix. A structural type, usable as a template argument, of strictly ascending row-major positions: two equal patterns are one and the same type of matrix. Prefer the [named patterns and factories](#patterns), which establish that order.

#### Declaration

```cpp
namespace fcarouge::sparse {
struct position;

inline constexpr std::size_t absent;

template <std::size_t Rows, std::size_t Columns, std::size_t Count>
struct pattern;
}
```

#### Template Parameters

| Template Parameter | Definition |
| --- | --- |
| `Rows` | The count of rows of the matrix. |
| `Columns` | The count of columns of the matrix. |
| `Count` | The count of stored elements of the matrix. |

#### Member Variables

| Member Variable | Definition |
| --- | --- |
| `rows` | The count of rows. |
| `columns` | The count of columns. |
| `count` | The count of stored elements. |
| `positions` | The `std::array` of the stored positions, in row-major order. A `position` has a `row` and a `column`. |

#### Member Functions

| Member Function | Definition |
| --- | --- |
| `lower_bound(position)` | The index of the first stored position not before the given one. |
| `index_of(row, column)` | The index of the stored element at the given position, or `absent` for a structural zero. |
| `stores(row, column)` | Whether the given position is a stored element. |
| `operator==` | Whether two patterns, of any sizes, store the same positions of the same shape. |

### Patterns

```cpp
namespace fcarouge::sparse {
template <std::size_t Rows, std::size_t Columns>
consteval auto make_pattern();

template <std::size_t Rows, std::size_t Columns, std::size_t Count>
consteval auto make_pattern(const position (&positions)[Count]);

template <std::size_t Rows, std::size_t Columns, auto Predicate>
consteval auto generate();
}
```

| Pattern | Definition |
| --- | --- |
| `make_pattern<R, C>()` | No position. |
| `make_pattern<R, C>({{i, j}, ...})` | The listed positions, in any order. A duplicate or out-of-range position does not compile. |
| `generate<R, C, predicate>()` | The positions satisfying a stateless predicate of a row and a column. |
| `dense<R, C>` | Every position: a matrix without structural zeros. |
| `zero<R, C>` | No position: the zero matrix. |
| `banded<R, C, L, U>` | The `L` sub-diagonals, the main diagonal, and the `U` super-diagonals. |
| `diagonal<N>` | The main diagonal of a square matrix. |
| `tridiagonal<N>` | The three central diagonals of a square matrix. |
| `lower_triangular<N>` | The main diagonal and below of a square matrix. |
| `upper_triangular<N>` | The main diagonal and above of a square matrix. |

### Pattern Algebra

The compile-time computations the operations derive their results with, each a variable template of the `fcarouge::sparse` namespace.

| Variable | Definition |
| --- | --- |
| `pattern_union<Lhs, Rhs>` | The positions stored by either pattern of the same shape: the pattern of a sum or a difference. |
| `pattern_intersection<Lhs, Rhs>` | The positions stored by both patterns of the same shape. |
| `pattern_transpose<Pattern>` | The transposed pattern. |
| `pattern_product<Lhs, Rhs>` | The structural pattern of a matrix product: a position is stored whenever a term reaches it. |
| `product_term_count<Lhs, Rhs>` | The count of elementary products of a matrix product. |
| `is_subset_of<Inner, Outer>` | Whether a pattern stores every position of another, of the same shape. |

### Element Construction

#### `make_matrix`

```cpp
namespace fcarouge::sparse {
struct structural_zero_t;

inline constexpr structural_zero_t structural_zero;

template <typename Matrix, typename... Values>
[[nodiscard]] constexpr auto make_matrix(const Values &...values) -> Matrix;
}
```

Construct a matrix from a value at every position, in row-major order. Does not participate in overload resolution unless `Matrix` is a sparse matrix, every position has a value, every stored element's value is implicitly convertible to the element type, and every structural zero's value is `structural_zero`: a `0.` at a structural zero does not compile, nor does `structural_zero` at a stored element.

```cpp
using lower = matrix<double, lower_triangular<2>>;

const lower value{make_matrix<lower>(1., structural_zero,
                                     2., 3.)};
```

| Parameter | Definition |
| --- | --- |
| `Matrix` | The sparse matrix type to construct. |
| `values` | The value of each position, in row-major order. |

### Element Access

A structural zero has nothing to refer to. A read is by value, and a write that would be lost does not compile, whatever the value. A pattern without a structural zero is an ordinary dense storage.

| Access | Stored Element | Structural Zero |
| --- | --- | --- |
| Read: `at<i, j>()`, `operator[i, j]` | Yes. | Yes, by value. A zero of the element type. |
| Write: `at<i, j>() = value` | Yes. | Does not compile. |
| Write: `operator[i, j] = value` | Without a structural zero only: its position is only known at run time. | Does not compile. |
| Construction from the stored elements | Yes, in row-major order. | Takes no value. |
| `make_matrix<Matrix>(values...)` | Yes, a value. | Yes, `structural_zero`. |
| Conversion from a matrix of a narrower pattern | Yes. | Yes, zero. |

The rejections of the writes and of the constructions are compile errors, which a `requires` expression detects. To ask whether a position can be written, use `decltype(value)::structure.stores(i, j)`.

### Operations

The operations below are supported, each as a function or an operator of the `fcarouge::sparse` namespace, found by argument-dependent lookup. `K` is the count of stored elements, `R` and `C` the counts of rows and columns, and `T` the [`product_term_count`](#pattern-algebra). Each result pattern is computed at compile time. The algorithms mirror the names of the C++26 standard linear algebra library, `std::linalg`, and return their results by value.

| Operation | Result Pattern | Run-Time Complexity |
| --- | --- | --- |
| `+`, binary `-` | The union of the patterns, of the same shape. | `O(K)` of the result. |
| Unary `-` | The same pattern. | `O(K)`. |
| `*`, of matrices and vectors | The structural product of the patterns. | `O(K + T)` of the result. |
| `*`, `/`, by a scalar | The same pattern. | `O(K)`. |
| `==` | Not applicable. Compares every position of the same shape: a structural zero equals a stored zero. | `O(K log K)` of the union. |
| `transposed` | The transposed pattern. Returns a matrix, not a view. | `O(K)`. |
| `dot` | Not applicable. Returns an element. Vectors of the same size, of any orientation. | `O(K)`, of the elements stored by both vectors. |
| `matrix_frob_norm` | Not applicable. Returns an element. Floating-point elements only. | `O(K)`. |
| `matrix_one_norm` | Not applicable. Returns an element, the largest absolute column sum. | `O(K + C)`. |
| `matrix_inf_norm` | Not applicable. Returns an element, the largest absolute row sum. | `O(K + R)`. |
| `make_matrix` | The pattern of the matrix to construct. | `O(R * C)`. |
| Read of an element at a run-time position | Not applicable. Returns an element. | `O(log K)`. |
| Access of an element at a compile-time position | Not applicable. Returns an element. | `O(1)`. |

The following are not supported yet:

- Matrix division, which needs a solver, a factorization.
- The `std::linalg` out-parameter algorithms: `add`, `matrix_product`, `matrix_rank_1_update`, `matrix_vector_product`, and `scale`.

### Concepts

| Concept | Definition |
| --- | --- |
| `sparsity_pattern` | Concept of a compile-time sparsity pattern. |

### Format

A matrix formats every element, structural zeros included, so the format does not depend on the pattern: a singleton as its element, a row vector as a list of elements, a column vector or a matrix as a list of rows. Hosted implementation only.

```cpp
const matrix<double, diagonal<2>> value{1., 2.5};

std::println("{}", value); // [[1, 0], [0, 2.5]]
```

## More

### Use Cases

The library serves fixed-size computations whose structure is known at compile time and whose budget is tight: embedded, real-time, safety-critical, or high-rate software, where a dense computation wastes most of its operations on zeros and a run-time sparse computation is not an option.

- **Estimation and filtering**: diagonal noises, triangular or block transitions, and selection observation models of Kalman filters.
- **Guidance, navigation, and control**: decoupled axes and cascaded loops whose gain and dynamics matrices are block or band structured.
- **Simulation**: finite difference stencils, banded Jacobians, and explicit time integration.
- **Robotics**: kinematic chains whose Jacobians only couple neighboring joints.
- **Linear algebra backend**: a storage of fixed patterns beneath a typed linear algebra facade, such as [TypedLinearAlgebra](https://github.com/FrancoisCarouge/TypedLinearAlgebra).

### Notes

- A pattern is structural: elements cancelling out at run time stay stored.
- A structural zero stays an exact zero, for which no term exists. Unlike a dense computation, a non-finite operand does not propagate through it: `0 * infinity` is never computed, and is not a not-a-number value.
- A matrix converts implicitly to a pattern storing all of its positions, and never to a narrower one, which could lose an element. Matrices of different element types do not mix: a conversion must be decided by the user.
- `sizeof` of a matrix with stored elements is that of its `std::array` of stored elements. A matrix of an empty pattern still occupies the storage of an empty `std::array`, which is not guaranteed to be nothing.
- The library is freestanding. On GCC and Clang, the `sparse_linear_algebra_freestanding` target verifies it with `-ffreestanding -fno-exceptions -fno-rtti`, and the continuous integration cross-compiles it for a bare-metal Arm Cortex-M4.
- At compile time, and on a freestanding implementation, the Frobenius norm uses its own square root, whose result is not guaranteed to be identical to the last digit to that of `std::sqrt`.
- The compile-time cost suits tens to low hundreds of stored elements. The pattern operations are near-linear in the stored elements and verified at 100 x 100, but a product precomputes each of its elementary products: a 16 x 16 dense product is verified, while a 32 x 32 dense product exceeds Clang's default constant evaluation limit.

### Roadmap & Vision

The library's patterns, pattern algebra, and matrix are stable in practice but not yet frozen. The longer-term ambition is to be the compile-time sparse storage beneath the [TypedLinearAlgebra](https://github.com/FrancoisCarouge/TypedLinearAlgebra) facade, where dimensional and semantic type safety meets structural sparsity, and beneath the [Kalman](https://github.com/FrancoisCarouge/Kalman) filters.

Toward version 1.0:

- Integrate as a TypedLinearAlgebra backend, with its strong element types.
- Solve: triangular substitutions and factorizations, for the matrix division.
- Interoperate with `std::mdspan` and the `std::linalg` out-parameter algorithms.
- Compose block patterns from smaller patterns.
- Lower the compile-time cost of large products.
- Stabilize the public API.

### Lessons Learned

**The pattern is in the type**, not in the object. A run-time sparse format stores its indexes beside its values and pays an indirection for each element. A compile-time pattern stores nothing but its values: the indexes, the result patterns, and the exact list of elementary products of each operation are constants of the program, computed once, by the compiler.

**A structural zero has nothing to refer to.** A reference to a structural zero would either write to nowhere, silently losing a value, or need a run-time check. Element writes are therefore compile-time positions only, and a run-time position is read by value. A pattern without structural zeros is an ordinary dense storage, with every element addressable.

**A structural zero is not a stored zero.** It is a property of the type: the product of a structural zero and an infinity is never computed. A stored zero is a value, which may become non-zero, and multiplies as any other value.

**Construction states the structure.** A value given for a structural zero, even a zero, is a defect: it says the author expected an element there. The `structural_zero` tag makes the intent explicit and checks it at compile time, at no run-time cost.

### Performance

The cost of an operation follows its count of stored elements and of elementary products, not its shape. See the [benchmarks](https://github.com/FrancoisCarouge/SparseLinearAlgebra/blob/master/benchmark/README.md).

### Projects

The library is used in projects:

*Your project link here!*

### Resources

1. D. Withopf, "Compile-time Is the New Constexpr: Leveraging Compile-time Sparsity for Vectors and Matrices," presented at C++Now, Aspen, CO, USA, May 2023. [Online]. Available: <https://schedule.cppnow.org/session/2023/compile-time-is-the-new-constexpr-leveraging-compile-time-sparsity-for-vectors-and-matrices/>
2. M. Hoemmen, "std::linalg: Linear Algebra Coming to Standard C++" (WG21 P1673), presented at CppCon, Aurora, CO, USA, Oct. 2023. [Online]. Available: <https://www.youtube.com/watch?v=-UXHMlAMXNk>
3. T. A. Davis, *Direct Methods for Sparse Linear Systems*. Philadelphia, PA, USA: SIAM, 2006.

### Third Party Acknowledgement

The library is designed, developed, and tested with the help of third-party tools and services acknowledged and thanked here:

- [actions-gh-pages](https://github.com/peaceiris/actions-gh-pages) to upload the documentation to GitHub pages.
- [Clang](https://clang.llvm.org) for compilation and code sanitizers.
- [clang-format](https://clang.llvm.org/docs/ClangFormat.html) for code formatting.
- [clang-tidy](https://clang.llvm.org/extra/clang-tidy/) for static analysis.
- [CMake](https://cmake.org) for build automation.
- [cmakefmt](https://cmakefmt.dev) for pretty CMake list files.
- [Coverity Scan](https://scan.coverity.com/projects/francoiscarouge-sparselinearalgebra) for static analysis.
- [cppcheck](https://cppcheck.sourceforge.io) for static analysis.
- [Doxygen](https://doxygen.nl) for documentation generation.
- [Doxygen Awesome](https://github.com/jothepro/doxygen-awesome-css) for pretty documentation.
- [GCC](https://gcc.gnu.org) for compilation and code sanitizers.
- [lcov](http://ltp.sourceforge.net/coverage/lcov.php) to process coverage information.
- [MSVC](https://docs.microsoft.com/en-US/cpp/windows/latest-supported-vc-redist) for compilation and code sanitizers.
- [nanobench](https://github.com/martinus/nanobench) for benchmarking.
- [Valgrind](https://valgrind.org) to check for correct memory management.

### Contributors

Thank you to everyone who has contributed code, tests, or documentation!

- [François Carouge](https://github.com/FrancoisCarouge) — creator and maintainer.

See [CONTRIBUTORS.md](https://github.com/FrancoisCarouge/SparseLinearAlgebra/blob/master/CONTRIBUTORS.md) for the full, credited list.

### Sponsors

Become a sponsor today! Support this project with coffee and infrastructure!

[![Sponsor](https://img.shields.io/badge/Support-Sponsor-brightgreen)](http://paypal.me/francoiscarouge)

#### Corporations & Institutions

*Your group logo and link here!*

#### Individuals

*Your name and link here!*

Thanks everyone!

### Commercial Support & Consulting

Need help integrating, tuning, or extending sparse linear algebra in your product? Commercial support, consulting, training, and custom development are available. [Get in touch](https://www.linkedin.com/in/francoiscarouge/).

## License

<img align="right" src="http://opensource.org/trademarks/opensource/OSI-Approved-License-100x137.png">

SparseLinearAlgebra is public domain:

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

For more information, please refer to <https://unlicense.org>

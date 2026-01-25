# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

### Primary Build Targets
- `make clang-san`: Build and test with Clang using AddressSanitizer and UBSan using libstdc++
- `make clang-no-san`: Build and test with Clang without sanitizers using libstdc++
- `make clang-san-libcpp`: Build and test with Clang using AddressSanitizer and UBSan using libc++
- `make clang-no-san-libcpp`: Build and test with Clang without sanitizers using libc++

### Architecture-Specific Builds
- `make clang-avx2`: Build targeting AVX2 instructions (x86-64-v3)
- `make clang-avx512`: Build targeting AVX512 instructions (x86-64-v4) with libc++
- `make clang-base-arch`: Build for baseline x86-64 architecture with libc++

### Other Builds
- `make clang-release`: Release build with Clang
- `make clang-type`: Build with TypeSanitizer (Clang only)
- `make clang-bench`: Build benchmarks with Clang

### Testing and Linting
- Tests run automatically with build targets
- `clang-tidy path/to/file`: Run the linter on specific files
- Individual test executables are in `build-*/` directories
- Run specific test patterns: `./build-clang/no-san/test_name`
- Benchmark executables are in `build-*/` directories

### Coverage Analysis
- `make coverage-local`: Generate local HTML coverage reports with lcov
- `make coverage-clean`: Clean coverage data files
- `make coverage-check`: Verify coverage meets 85% threshold
- Coverage reports available at: `build-clang/coverage/html/index.html`
- Coverage data excludes test files, external dependencies, and system headers

### Clean
- `make clean`: Remove all build directories

## Code Architecture

### Build System
The library uses C++23 modules:
- Module interface files (`.cxxm`) in `mod/` define module interfaces
- Module implementation files (`.cxx`) in `mod/` provide implementations for complex modules
- Headers in `include/` provide utility macros and allocator overrides
- Requires CMake 3.31+ for C++23 module support

### Core Components

**Math Module (`mod/Math/`)**
- `Array.cxxm`: Core array types, expression templates, and operations
- `ArrayConcepts.cxxm`: C++20 concepts for type safety and constraints
- `StaticArrays.cxxm`: Compile-time sized arrays for stack allocation
- `ManagedArray.cxxm`: Runtime-sized arrays with dynamic memory management
- `LinearAlgebra.cxxm`: Matrix operations, LU decomposition, etc.
- `Simplex.cxxm`: Rational simplex solver for linear programming
- `NormalForm.cxxm`: Matrix normal form algorithms
- `Constraints.cxxm`: Constraint handling for polyhedral analysis
- `Reductions.cxxm`, `Comparisons.cxxm`: Array reduction and comparison operations

**SIMD Module (`mod/SIMD/`)**
- `SIMD.cxxm`: SIMD vector types using Clang's `ext_vector_type`
- Platform-specific optimizations for x86-64 (SSE, AVX, AVX512)
- Automatic SIMD width detection and vectorization

**Container Module (`mod/Containers/`)**
- `TinyVector.cxxm`: Small vector optimizations
- `Storage.cxxm`: Memory layout and allocation strategies
- `Tuple.cxxm`: Multi-dimensional indexing utilities
- `BitSets.cxxm`: Efficient bit manipulation containers
- `UnrolledList.cxxm`: Unrolled linked list implementation
- `Permutation.cxxm`: Permutation utilities
- `Buffer.cxxm`: Buffer container

**Numbers Module (`mod/Numbers/`)**
- `Int8.cxxm`: Integer type utilities

**Allocation Module (`mod/Alloc/`)**
- `Arena.cxxm`: Arena-based memory allocator
- `Mallocator.cxxm`: Custom memory allocation strategies

**Utilities Module (`mod/Utilities/`)**
- `ArrayPrint.cxxm`: Array printing and formatting utilities
- `MatrixStringParse.cxxm`: Parsing matrices from string representations
- `Optional.cxxm`: Optional type implementations
- `Valid.cxxm`: Validation utilities
- `Sort.cxxm`: Sorting utilities

### Key Design Patterns

**Expression Templates**
- Operations build expression trees evaluated on assignment
- Use `<<` operator for data copying: `A << B + C * D`
- Regular `=` operator reserved for view/reference assignment
- Enables operation fusion and SIMD vectorization

**Strong Typing System**
- `Length<N>`: 1D vectors with compile/runtime size
- `DenseDims<M,N>`: Dense matrices
- `StridedDims<M,N,X>`: Strided matrices  
- `SquareDims<N>`: Square matrices
- Prevents dimension mismatches at compile time

**Array Types**
- `Array<T,S>`: Immutable array view
- `MutArray<T,S>`: Mutable array view
- `StaticArray<T,M,N>`: Stack-allocated arrays
- `ManagedArray<T,S>`: Heap-allocated with RAII

**Slicing and Indexing**
- Julia-style slicing: `_` (entire slice), `_(i,j)` (range)
- Keywords: `end`, `last` for boundary access
- Multi-dimensional: `A[i,j]` for matrices, `A[i,_]` for row views
- Creates views without copying data

### Integer-Focused Math
- Rational arithmetic for exact computations  
- Greatest Common Divisor (GCD) algorithms with SIMD
- Linear Diophantine equation solving
- Unimodularization and normal form algorithms
- Polyhedral analysis utilities

### Memory Management
- Arena allocator for efficient bulk allocation
- RAII-based lifetime management
- Support for custom allocators (mimalloc, jemalloc)
- Trivial type optimizations

## Compiler Requirements
- **Clang 21+** required (GCC not supported due to use of `ext_vector_type`)
- **CMake 3.31+** required for C++23 module support
- Uses cutting-edge features: deducing this, C++23 modules
- Multiple sanitizer configurations for robust testing
- Native compilation support with `-march=native`

## Testing Strategy
- Comprehensive ut (micro test) framework in `test/` directory
- Architecture-specific testing (baseline, AVX2, AVX512)
- Sanitizer builds catch memory/undefined behavior issues
- Performance benchmarks in `benchmark/` directory
- Individual test executables for each test file

## Development Notes
- Library uses C++23 modules
- Expression templates require careful attention to lifetime management
- SIMD operations use Clang's `ext_vector_type` for portable vectorization
- Focus on integer operations distinguishes this from typical math libraries
- Used by LoopModels for polyhedral compilation analysis

## Important Development Considerations
- The build system automatically detects CPU features (AVX2, AVX512) and adjusts targets accordingly
- All builds use `-fno-exceptions` and `-fno-rtti` for performance
- Memory allocators can be configured (mimalloc, jemalloc) via CMake options
- Sanitizer builds help catch subtle bugs in template-heavy code
- Use ut framework for testing with individual test executables per test file

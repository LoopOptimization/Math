# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

### Primary Build Targets
- `make clang-san`: Build and test with Clang using AddressSanitizer and UndefinedBehaviorSanitizer using libc++
- `make clang-no-san`: Build and test with Clang without sanitizers using libc++
- `make clang-san-libstdcxx`: Build and test with Clang using AddressSanitizer and UndefinedBehaviorSanitizer using libstdc++
- `make clang-no-san-libstdcxx`: Build and test with Clang without sanitizers using libstdc++
- `make gcc-san`: Build and test with GCC using AddressSanitizer and UndefinedBehaviorSanitizer
- `make gcc-no-san`: Build and test with GCC without sanitizers

### Architecture-Specific Builds
- `make clang-avx512`: Build targeting AVX512 instructions
- `make clang-base-arch`: Build for baseline x86-64 architecture
- `make gcc-avx2`: Build targeting AVX2 instructions

### Other Builds
- `make clang-release`: Release build with Clang
- `make gcc-release`: Release build with GCC
- `make clang-nosimd`: Build without explicit SIMD operations
- `make clang-type`: Build with TypeSanitizer (Clang only)
- `make clang-bench`: Build benchmarks with Clang
- `make gcc-bench`: Build benchmarks with GCC

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

### Dual Build System
The library supports both traditional header-only compilation and modern C++23 modules:
- Module files (`.cxxm`) in `mod/` serve as module interfaces
- Module implementation files (`.cxx`) in `mod/` provide implementations
- Headers in `include/` provide the main API entry points
- CMake automatically detects and builds C++23 modules when using compatible compilers

### Core Components

**Math Module (`mod/Math/`)**
- `Array.cxx`: Core array types, expression templates, and operations
- `ArrayConcepts.cxx`: C++20 concepts for type safety and constraints  
- `StaticArrays.cxx`: Compile-time sized arrays for stack allocation
- `ManagedArray.cxx`: Runtime-sized arrays with dynamic memory management
- `LinearAlgebra.cxx`: Matrix operations, LU decomposition, etc.
- `Simplex.cxx`: Rational simplex solver for linear programming

**SIMD Module (`mod/SIMD/`)**
- `Vec.cxx`: SIMD vector types using GCC vector extensions
- `Intrin.cxx`: Cross-platform intrinsic wrappers
- Platform-specific optimizations for x86-64 (SSE, AVX, AVX512)
- Automatic SIMD width detection and vectorization

**Container Module (`mod/Containers/`)**
- `TinyVector.cxx`: Small vector optimizations
- `Storage.cxx`: Memory layout and allocation strategies
- `Tuple.cxx`: Multi-dimensional indexing utilities
- `BitSets.cxx`: Efficient bit manipulation containers
- `UnrolledList.cxx`: Unrolled linked list implementation
- `Flat.cxx`: Flattened container structures

**Allocation Module (`mod/Alloc/`)**
- `Arena.cxx`: Arena-based memory allocator
- `Mallocator.cxx`: Custom memory allocation strategies

**Utilities Module (`mod/Utilities/`)**
- `ArrayPrint.cxx`: Array printing and formatting utilities
- `MatrixStringParse.cxx`: Parsing matrices from string representations
- `Optional.cxx`: Optional type implementations
- `Valid.cxx`: Validation utilities

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
- **Clang 17+** or **GCC 13.2+** required for C++23 features
- Uses cutting-edge features: deducing this, C++23 modules
- Multiple sanitizer configurations for robust testing
- Native compilation support with `-march=native`

## Testing Strategy
- Comprehensive ut (micro test) framework in `test/` directory
- Cross-compiler testing (Clang/GCC) with multiple configurations
- Architecture-specific testing (baseline, AVX2, AVX512)
- Sanitizer builds catch memory/undefined behavior issues
- Performance benchmarks in `benchmark/` directory
- Individual test executables for each test file

## Development Notes
- Library is header-only by default but can build as C++23 modules
- Expression templates require careful attention to lifetime management
- SIMD operations are automatically vectorized where possible
- Focus on integer operations distinguishes this from typical math libraries
- Used by LoopModels for polyhedral compilation analysis

## Important Development Considerations
- The build system automatically detects CPU features (AVX2, AVX512) and adjusts targets accordingly
- All builds use `-fno-exceptions` and `-fno-rtti` for performance
- Memory allocators can be configured (mimalloc, jemalloc) via CMake options
- Cross-compiler testing with both Clang and GCC is essential due to C++23 features
- Sanitizer builds help catch subtle bugs in template-heavy code
- Use ut framework for testing with individual test executables per test file
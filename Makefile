HAVE_AVX512 := $(shell grep avx512 /proc/cpuinfo &> /dev/null; echo $$?)
HAVE_AVX2 := $(shell grep avx2 /proc/cpuinfo &> /dev/null; echo $$?)

ifeq ($(HAVE_AVX512),0)
all: clang-modules clang-no-san clang-san gcc-no-san gcc-san clang-nosimd clang-base-arch clang-release gcc-release gcc-avx2 clang-avx512
else ifeq ($(HAVE_AVX2),0)
all: clang-modules clang-no-san clang-san gcc-no-san gcc-san clang-nosimd clang-base-arch clang-release gcc-release gcc-avx2
else
all: clang-modules clang-no-san clang-san gcc-no-san gcc-san clang-nosimd clang-base-arch clang-release gcc-release
endif
#TODO: re-enable GCC once multidimensional indexing in `requires` is fixed:
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=111493

# `command -v` returns nothing if not found (and we redirect stderr)
NINJA := $(shell command -v ninja 2> /dev/null)
ifdef NINJA
    NINJAGEN := "-G Ninja"
else
    NINJAGEN := ""
endif


build-gcc/no-san/:
	CXXFLAGS="-Og" CXX=g++ cmake $(NINJAGEN) -S test -B build-gcc/no-san/ -DCMAKE_BUILD_TYPE=Debug

build-gcc/san/:
	CXXFLAGS="" CXX=g++ cmake $(NINJAGEN) -S test -B build-gcc/san/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined' -DPOLYMATHNOEXPLICITSIMDARRAY=OFF

build-clang/no-san/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/no-san/ -DCMAKE_BUILD_TYPE=Debug

build-clang/san/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/san/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined'
	
build-gcc/avx2/:
	CXXFLAGS="-Og -march=x86-64-v3" CXX=g++ cmake $(NINJAGEN) -S test -B build-gcc/avx2/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build-gcc/modules/:
	CXXFLAGS="" CXX=g++ cmake $(NINJAGEN) -S test -B build-gcc/modules/ -DCMAKE_BUILD_TYPE=Debug -DUSE_MODULES=ON

build-clang/base-arch/:
	CXXFLAGS="-Og -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/base-arch/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build-clang/avx512/:
	CXXFLAGS="-Og -march=x86-64-v4 -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/avx512/ -DCMAKE_BUILD_TYPE=Debug

build-clang/no-simd/:
	CXXFLAGS="-Og -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/no-simd/ -DCMAKE_BUILD_TYPE=Debug -DPOLYMATHNOEXPLICITSIMDARRAY=ON

build-clang/modules/:
	CXXFLAGS="-Og -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/modules/ -DCMAKE_BUILD_TYPE=Debug -DUSE_MODULES=ON

build-clang/bench/:
	CXXFLAGS="-Og -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S benchmark -B build-clang/bench/ -DCMAKE_BUILD_TYPE=Release

build-gcc/bench/:
	CXXFLAGS="-Og" CXX=g++ cmake $(NINJAGEN) -S benchmark -B build-gcc/bench/ -DCMAKE_BUILD_TYPE=Release

build-clang/type/:
	CXXFLAGS="-Og -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/type/ -DCMAKE_BUILD_TYPE=Debug -DUSE_TYPE_SANITIZER=ON

build-clang/release/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/release/ -DCMAKE_BUILD_TYPE=RELEASE

build-gcc/release/:
	CXXFLAGS="" CXX=g++ cmake $(NINJAGEN) -S test -B build-gcc/release/ -DCMAKE_BUILD_TYPE=RELEASE

gcc-no-san: build-gcc/no-san/
	cmake --build build-gcc/no-san/
	cmake --build build-gcc/no-san/ --target test

gcc-san: build-gcc/san/
	cmake --build build-gcc/san/
	cmake --build build-gcc/san/ --target test

clang-no-san: build-clang/no-san/
	cmake --build build-clang/no-san/
	cmake --build build-clang/no-san/ --target test

clang-san: build-clang/san/
	cmake --build build-clang/san/
	cmake --build build-clang/san/ --target test

gcc-avx2: build-gcc/avx2/
	cmake --build build-gcc/avx2/
	cmake --build build-gcc/avx2/ --target test

clang-base-arch: build-clang/base-arch/
	cmake --build build-clang/base-arch/
	cmake --build build-clang/base-arch/ --target test

clang-avx512: build-clang/avx512/
	cmake --build build-clang/avx512/
	cmake --build build-clang/avx512/ --target test

clang-nosimd: build-clang/no-simd/
	cmake --build build-clang/no-simd/
	cmake --build build-clang/no-simd/ --target test

clang-modules: build-clang/modules/
	cmake --build build-clang/modules/
	cmake --build build-clang/modules/ --target test

gcc-modules: build-gcc/modules/
	cmake --build build-gcc/modules/
	cmake --build build-gcc/modules/ --target test

clang-bench: build-clang/bench/
	cmake --build build-clang/bench

gcc-bench: build-gcc/bench/
	cmake --build build-gcc/bench

clang-type: build-clang/type/
	cmake --build build-clang/type
	TYSAN_OPTIONS=print_stacktrace=1 cmake --build build-clang/type --target test

clang-release: build-clang/release/
	cmake --build build-clang/release
	cmake --build build-clang/release --target test

gcc-release: build-gcc/release/
	cmake --build build-gcc/release
	cmake --build build-gcc/release --target test


clean:
	rm -rf build-clang build-gcc

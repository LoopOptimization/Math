HAVE_AVX512 := $(shell grep avx512 /proc/cpuinfo &> /dev/null; echo $$?)
HAVE_AVX2 := $(shell grep avx2 /proc/cpuinfo &> /dev/null; echo $$?)

ifeq ($(HAVE_AVX512),0)
all: clangmodules clangnosan clangsan gccnosan gccsan clangnosimd clangbasearch clangrelease gccrelease gccavx2 clangavx512
else ifeq ($(HAVE_AVX2),0)
all: clangmodules clangnosan clangsan gccnosan gccsan clangnosimd clangbasearch clangrelease gccrelease gccavx2
else
all: clangmodules clangnosan clangsan gccnosan gccsan clangnosimd clangbasearch clangrelease gccrelease
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


buildgcc/nosan/:
	CXXFLAGS="-Og" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/nosan/ -DCMAKE_BUILD_TYPE=Debug

buildgcc/test/:
	CXXFLAGS="" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/test/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined' -DPOLYMATHNOEXPLICITSIMDARRAY=OFF

buildclang/nosan/:
	CXXFLAGS="" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/nosan/ -DCMAKE_BUILD_TYPE=Debug

buildclang/test/:
	CXXFLAGS="" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/test/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined'
	
buildgcc/avx2/:
	CXXFLAGS="-Og -march=x86-64-v3" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/avx2/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

buildgcc/modules/:
	CXXFLAGS="" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/modules/ -DCMAKE_BUILD_TYPE=Debug -DUSE_MODULES=ON

buildclang/basearch/:
	CXXFLAGS="-Og" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/basearch/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

buildclang/avx512/:
	CXXFLAGS="-Og -march=x86-64-v4" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/avx512/ -DCMAKE_BUILD_TYPE=Debug

buildclang/nosimdarrayop/:
	CXXFLAGS="-Og" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/nosimdarrayop/ -DCMAKE_BUILD_TYPE=Debug -DPOLYMATHNOEXPLICITSIMDARRAY=ON

buildclang/modules/:
	CXXFLAGS="-Og" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/modules/ -DCMAKE_BUILD_TYPE=Debug -DUSE_MODULES=ON

buildclang/bench/:
	CXXFLAGS="-Og" CXX=clang++ cmake $(NINJAGEN) -S benchmark -B buildclang/bench/ -DCMAKE_BUILD_TYPE=Release

buildgcc/bench/:
	CXXFLAGS="-Og" CXX=g++ cmake $(NINJAGEN) -S benchmark -B buildgcc/bench/ -DCMAKE_BUILD_TYPE=Release

buildclang/type/:
	CXXFLAGS="-Og" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/type/ -DCMAKE_BUILD_TYPE=Debug -DUSE_TYPE_SANITIZER=ON

buildclang/release/:
	CXXFLAGS="" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/release/ -DCMAKE_BUILD_TYPE=RELEASE

buildgcc/release/:
	CXXFLAGS="" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/release/ -DCMAKE_BUILD_TYPE=RELEASE

gccnosan: buildgcc/nosan/
	cmake --build buildgcc/nosan/
	cmake --build buildgcc/nosan/ --target test

gccsan: buildgcc/test/
	cmake --build buildgcc/test/ 
	cmake --build buildgcc/test/ --target test

clangnosan: buildclang/nosan/
	cmake --build buildclang/nosan/
	cmake --build buildclang/nosan/ --target test

clangsan: buildclang/test/
	cmake --build buildclang/test/ 
	cmake --build buildclang/test/ --target test

gccavx2: buildgcc/avx2/
	cmake --build buildgcc/avx2/
	cmake --build buildgcc/avx2/ --target test

clangbasearch: buildclang/basearch/
	cmake --build buildclang/basearch/
	cmake --build buildclang/basearch/ --target test

clangavx512: buildclang/avx512/
	cmake --build buildclang/avx512/
	cmake --build buildclang/avx512/ --target test

clangnosimd: buildclang/nosimdarrayop/
	cmake --build buildclang/nosimdarrayop/
	cmake --build buildclang/nosimdarrayop/ --target test

clangmodules: buildclang/modules/
	cmake --build buildclang/modules/
	cmake --build buildclang/modules/ --target test

gccmodules: buildgcc/modules/
	cmake --build buildgcc/modules/
	cmake --build buildgcc/modules/ --target test

clangbench: buildclang/bench/
	cmake --build buildclang/bench

gccbench: buildgcc/bench/
	cmake --build buildgcc/bench

clangtype: buildclang/type/
	cmake --build buildclang/type
	TYSAN_OPTIONS=print_stacktrace=1 cmake --build buildclang/type --target test

clangrelease: buildclang/release/
	cmake --build buildclang/release
	cmake --build buildclang/release --target test

gccrelease: buildgcc/release/
	cmake --build buildgcc/release
	cmake --build buildgcc/release --target test


clean:
	rm -rf buildclang buildgcc

all: clangnosan clangsan gccnosan gccsan gccavx2 clangavx512 clangnosimd
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
	CXXFLAGS="" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/nosan/ -DCMAKE_BUILD_TYPE=Debug

buildgcc/test/:
	CXXFLAGS="-Og" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/test/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined' -DPOLYMATHNOEXPLICITSIMDARRAY=OFF

buildclang/nosan/:
	CXXFLAGS="" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/nosan/ -DCMAKE_BUILD_TYPE=Debug

buildclang/test/:
	CXXFLAGS="" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/test/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined'
	
buildgcc/avx2/:
	CXXFLAGS="-march=haswell" CXX=g++ cmake $(NINJAGEN) -S test -B buildgcc/avx2/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

buildclang/sse/:
	CXXFLAGS="" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/sse/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

buildclang/nosimdarrayop/:
	CXXFLAGS="" CXX=clang++ cmake $(NINJAGEN) -S test -B buildclang/nosimdarrayop/ -DCMAKE_BUILD_TYPE=Debug -DPOLYMATHNOEXPLICITSIMDARRAY=OFF


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

clangavx512: buildclang/sse/
	cmake --build buildclang/sse/
	cmake --build buildclang/sse/ --target test

clangnosimd: buildclang/nosimdarrayop/
	cmake --build buildclang/nosimdarrayop/
	cmake --build buildclang/nosimdarrayop/ --target test

clean:
	rm -rf buildclang #buildgcc

HAVE_AVX512 := $(shell grep avx512 /proc/cpuinfo &> /dev/null; echo $$?)
HAVE_AVX2 := $(shell grep avx2 /proc/cpuinfo &> /dev/null; echo $$?)

ifeq ($(HAVE_AVX512),0)
all: clang-no-san clang-san clang-no-san-libcpp clang-san-libcpp clang-base-arch clang-release clang-avx2 clang-avx512
else ifeq ($(HAVE_AVX2),0)
all: clang-no-san clang-san clang-no-san-libcpp clang-san-libcpp clang-base-arch clang-release clang-avx2
else
all: clang-no-san clang-san clang-no-san-libcpp clang-san-libcpp clang-base-arch clang-release
endif

# `command -v` returns nothing if not found (and we redirect stderr)
NINJA := $(shell command -v ninja 2> /dev/null)
ifdef NINJA
    NINJAGEN := "-G Ninja"
else
    NINJAGEN := ""
endif


build-clang/no-san-libcpp/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/no-san-libcpp/ -DCMAKE_BUILD_TYPE=Debug

build-clang/san-libcpp/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/san-libcpp/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined'

build-clang/no-san/:
	CXXFLAGS="-stdlib=libstdc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/no-san/ -DCMAKE_BUILD_TYPE=Debug

build-clang/san/:
	CXXFLAGS="-stdlib=libstdc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/san/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined'
	
build-clang/avx2/:
	CXXFLAGS="-march=x86-64-v3" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/avx2/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build-clang/base-arch/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/base-arch/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build-clang/avx512/:
	CXXFLAGS="-march=x86-64-v4 -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/avx512/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build-clang/bench/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S benchmark -B build-clang/bench/ -DCMAKE_BUILD_TYPE=Release

build-clang/type/:
	CXXFLAGS="-Og -stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/type/ -DCMAKE_BUILD_TYPE=Debug -DUSE_TYPE_SANITIZER=ON

build-clang/release/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/release/ -DCMAKE_BUILD_TYPE=RELEASE

clang-no-san-libcpp: build-clang/no-san-libcpp/
	cmake --build build-clang/no-san-libcpp/
	cmake --build build-clang/no-san-libcpp/ --target test

clang-san-libcpp: build-clang/san-libcpp/
	cmake --build build-clang/san-libcpp/
	cmake --build build-clang/san-libcpp/ --target test

clang-no-san: build-clang/no-san/
	cmake --build build-clang/no-san/
	cmake --build build-clang/no-san/ --target test

clang-san: build-clang/san/
	cmake --build build-clang/san/
	cmake --build build-clang/san/ --target test

clang-avx2: build-clang/avx2/
	cmake --build build-clang/avx2/
	cmake --build build-clang/avx2/ --target test

clang-base-arch: build-clang/base-arch/
	cmake --build build-clang/base-arch/
	cmake --build build-clang/base-arch/ --target test

clang-avx512: build-clang/avx512/
	cmake --build build-clang/avx512/
	cmake --build build-clang/avx512/ --target test

clang-bench: build-clang/bench/
	cmake --build build-clang/bench

clang-type: build-clang/type/
	cmake --build build-clang/type
	TYSAN_OPTIONS=print_stacktrace=1 cmake --build build-clang/type --target test

clang-release: build-clang/release/
	cmake --build build-clang/release
	cmake --build build-clang/release --target test


# Coverage builds
build-clang/coverage/:
	CXXFLAGS="-stdlib=libc++ -fprofile-instr-generate -fcoverage-mapping" CXX=clang++ cmake $(NINJAGEN) -S test -B build-clang/coverage/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_TEST_COVERAGE=ON

coverage-local: build-clang/coverage/
	cmake --build build-clang/coverage/
	LLVM_PROFILE_FILE="reports/%p.profraw" cmake --build build-clang/coverage/ --target test
	@echo "Generating coverage report..."
	@if command -v llvm-cov >/dev/null 2>&1; then \
		echo "Generating coverage report with llvm-cov..."; \
		llvm-profdata merge -sparse build-clang/coverage/reports/*.profraw -o build-clang/coverage/merged.profdata 2>/dev/null || true; \
		if [ -f build-clang/coverage/merged.profdata ]; then \
			echo "Finding test executables..."; \
			OBJECT_FLAGS=""; \
			for exe in build-clang/coverage/*test*; do \
				if [ -x "$$exe" ] && [ -f "$$exe" ]; then \
					OBJECT_FLAGS="$$OBJECT_FLAGS --object $$exe"; \
				fi; \
			done; \
			echo "Using object flags: $$OBJECT_FLAGS"; \
			if [ -n "$$OBJECT_FLAGS" ]; then \
				llvm-cov show $$OBJECT_FLAGS -instr-profile=build-clang/coverage/merged.profdata -format=html -output-dir=build-clang/coverage/html -ignore-filename-regex='extern|test'; \
				echo "Coverage report generated: build-clang/coverage/html/index.html"; \
			else \
				echo "No test executables found in build-clang/coverage/"; \
			fi; \
		else \
			echo "No coverage data found. Make sure ENABLE_TEST_COVERAGE=ON and tests were run."; \
		fi \
	else \
		echo "llvm-cov not found."; \
		echo "Install llvm-cov (usually in llvm package)"; \
		gcov build-clang/coverage/CMakeFiles/Math.dir/mod/**/*.gcno 2>/dev/null || echo "gcov also failed"; \
		echo "Basic coverage files may be in build-clang/coverage/"; \
	fi

coverage-clean:
	@echo "Cleaning coverage data..."
	@rm -f build-clang/coverage/CMakeFiles/Math.dir/mod/**/*.gcda
	@rm -f build-clang/coverage/CMakeFiles/*/test/*.gcda
	@rm -f build-clang/coverage/coverage.info
	@rm -rf build-clang/coverage/html

coverage-check: coverage-local
	@echo "Checking coverage thresholds..."
	@if command -v lcov >/dev/null 2>&1; then \
		lcov --summary build-clang/coverage/coverage.info | grep -E "lines.*: [0-9.]+%" | \
		awk '{print $$2}' | sed 's/%//' | \
		awk '{if($$1 < 85) {print "Coverage " $$1 "% below threshold of 85%"; exit 1} else {print "Coverage " $$1 "% meets threshold"}}'; \
	else \
		echo "lcov not available for threshold checking"; \
	fi

clean:
	rm -rf build-clang

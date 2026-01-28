HAVE_AVX512 := $(shell grep avx512 /proc/cpuinfo &> /dev/null; echo $$?)
HAVE_AVX2 := $(shell grep avx2 /proc/cpuinfo &> /dev/null; echo $$?)

ifeq ($(HAVE_AVX512),0)
all: no-san san no-san-libcpp san-libcpp base-arch release avx2 avx512
else ifeq ($(HAVE_AVX2),0)
all: no-san san no-san-libcpp san-libcpp base-arch release avx2
else
all: no-san san no-san-libcpp san-libcpp base-arch release
endif

build/no-san-libcxx/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake -G Ninja -S test -B build/no-san-libcxx/ -DCMAKE_BUILD_TYPE=Debug

build/san-libcxx/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake -G Ninja -S test -B build/san-libcxx/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined'

build/no-san/:
	CXXFLAGS="-stdlib=libstdc++" CXX=clang++ cmake -G Ninja -S test -B build/no-san/ -DCMAKE_BUILD_TYPE=Debug

build/san/:
	CXXFLAGS="-stdlib=libstdc++" CXX=clang++ cmake -G Ninja -S test -B build/san/ -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER='Address;Undefined'
	
build/avx2/:
	CXXFLAGS="-march=x86-64-v3" CXX=clang++ cmake -G Ninja -S test -B build/avx2/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build/base-arch/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake -G Ninja -S test -B build/base-arch/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build/avx512/:
	CXXFLAGS="-march=x86-64-v4 -stdlib=libc++" CXX=clang++ cmake -G Ninja -S test -B build/avx512/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_NATIVE_COMPILATION=OFF

build/bench/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake -G Ninja -S benchmark -B build/bench/ -DCMAKE_BUILD_TYPE=Release

build/type/:
	CXXFLAGS="-Og -stdlib=libc++" CXX=clang++ cmake -G Ninja -S test -B build/type/ -DCMAKE_BUILD_TYPE=Debug -DUSE_TYPE_SANITIZER=ON

build/release/:
	CXXFLAGS="-stdlib=libc++" CXX=clang++ cmake -G Ninja -S test -B build/release/ -DCMAKE_BUILD_TYPE=RELEASE

no-san-libcxx: build/no-san-libcxx/
	cmake --build build/no-san-libcxx/
	cmake --build build/no-san-libcxx/ --target test

san-libcxx: build/san-libcxx/
	cmake --build build/san-libcxx/
	cmake --build build/san-libcxx/ --target test

no-san: build/no-san/
	cmake --build build/no-san/
	cmake --build build/no-san/ --target test

san: build/san/
	cmake --build build/san/
	cmake --build build/san/ --target test

avx2: build/avx2/
	cmake --build build/avx2/
	cmake --build build/avx2/ --target test

base-arch: build/base-arch/
	cmake --build build/base-arch/
	cmake --build build/base-arch/ --target test

avx512: build/avx512/
	cmake --build build/avx512/
	cmake --build build/avx512/ --target test

bench: build/bench/
	cmake --build build/bench

type: build/type/
	cmake --build build/type
	TYSAN_OPTIONS=print_stacktrace=1 cmake --build build/type --target test

release: build/release/
	cmake --build build/release
	cmake --build build/release --target test


# Coverage builds
build/coverage/:
	CXXFLAGS="-stdlib=libc++ -fprofile-instr-generate -fcoverage-mapping" CXX=clang++ cmake -G Ninja -S test -B build/coverage/ -DCMAKE_BUILD_TYPE=Debug -DENABLE_TEST_COVERAGE=ON

coverage-local: build/coverage/
	cmake --build build/coverage/
	LLVM_PROFILE_FILE="reports/%p.profraw" cmake --build build/coverage/ --target test
	@echo "Generating coverage report..."
	@if command -v llvm-cov >/dev/null 2>&1; then \
		echo "Generating coverage report with llvm-cov..."; \
		llvm-profdata merge -sparse build/coverage/reports/*.profraw -o build/coverage/merged.profdata 2>/dev/null || true; \
		if [ -f build/coverage/merged.profdata ]; then \
			echo "Finding test executables..."; \
			OBJECT_FLAGS=""; \
			for exe in build/coverage/*test*; do \
				if [ -x "$$exe" ] && [ -f "$$exe" ]; then \
					OBJECT_FLAGS="$$OBJECT_FLAGS --object $$exe"; \
				fi; \
			done; \
			echo "Using object flags: $$OBJECT_FLAGS"; \
			if [ -n "$$OBJECT_FLAGS" ]; then \
				llvm-cov show $$OBJECT_FLAGS -instr-profile=build/coverage/merged.profdata -format=html -output-dir=build/coverage/html -ignore-filename-regex='extern|test'; \
				echo "Coverage report generated: build/coverage/html/index.html"; \
			else \
				echo "No test executables found in build/coverage/"; \
			fi; \
		else \
			echo "No coverage data found. Make sure ENABLE_TEST_COVERAGE=ON and tests were run."; \
		fi \
	else \
		echo "llvm-cov not found."; \
		echo "Install llvm-cov (usually in llvm package)"; \
		gcov build/coverage/CMakeFiles/Math.dir/mod/**/*.gcno 2>/dev/null || echo "gcov also failed"; \
		echo "Basic coverage files may be in build/coverage/"; \
	fi

coverage-clean:
	@echo "Cleaning coverage data..."
	@rm -f build/coverage/CMakeFiles/Math.dir/mod/**/*.gcda
	@rm -f build/coverage/CMakeFiles/*/test/*.gcda
	@rm -f build/coverage/coverage.info
	@rm -rf build/coverage/html

coverage-check: coverage-local
	@echo "Checking coverage thresholds..."
	@if command -v lcov >/dev/null 2>&1; then \
		lcov --summary build/coverage/coverage.info | grep -E "lines.*: [0-9.]+%" | \
		awk '{print $$2}' | sed 's/%//' | \
		awk '{if($$1 < 85) {print "Coverage " $$1 "% below threshold of 85%"; exit 1} else {print "Coverage " $$1 "% meets threshold"}}'; \
	else \
		echo "lcov not available for threshold checking"; \
	fi

clean:
	rm -rf build

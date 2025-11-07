cmake_minimum_required(VERSION 3.31 FATAL_ERROR)
# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(POLYMATHNOEXPLICITSIMDARRAY "No explicit SIMD for Array operations" OFF)
option(USE_TYPE_SANITIZER "Use type sanitizer" OFF)
option(ENABLE_NATIVE_COMPILATION "Compile with -march=native" ON)

# --- Native compilation ---

if(ENABLE_NATIVE_COMPILATION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    add_compile_options(-xhost)
    if(ENABLE_WIDE_VECTORS)
      add_compile_options(-qopt-zmm-usage=high)
    endif()
  else()
    add_compile_options(-march=native)
    if(ENABLE_WIDE_VECTORS)
      if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
          add_compile_options(-mprefer-vector-width=512)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
          add_compile_options(-mprefer-vector-width=512
                              -mtune-ctrl=avx512_move_by_pieces)
        endif()
      endif()
    endif()
  endif()
endif()
if(APPLE)
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:CXX>:-nostdinc++;-nostdlib++;-isystem;/opt/homebrew/opt/llvm/include/c++/v1>"
  )
  add_link_options(
    "$<$<COMPILE_LANGUAGE:CXX>:-L/opt/homebrew/opt/llvm/lib/c++;-Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++;-lc++>"
  )
endif()
add_compile_options(
  "$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-omit-frame-pointer;-fno-rtti>"
)
# add_compile_options("-save-temps") # for preprocessed source
if((USE_SANITIZER MATCHES "([Uu]ndefined)") OR (USE_SANITIZER MATCHES
                                                "([Aa]ddress);([Uu]ndefined)"))
  set(USE_MI_MALLOC
      OFF
      CACHE BOOL "No mimalloc with sanitizers")
  set(USE_JE_MALLOC
      OFF
      CACHE BOOL "No jemalloc with sanitizers")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fstack-protector-strong>")
  endif()
else()
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fstack-protector-strong>")
endif()
if(USE_TYPE_SANITIZER)
  set(USE_MI_MALLOC
      OFF
      CACHE BOOL "No mimalloc with sanitizers")
  set(USE_JE_MALLOC
      OFF
      CACHE BOOL "No jemalloc with sanitizers")
endif()
# set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
    # This specific value changes as experimental support evolves. See
    # `Help/dev/experimental.rst` in the CMake source corresponding to
    # your CMake build for the exact value to use.
    "d0edc3af-4c50-42ea-a356-e2862fe7a444")
# cmake master uses: "d0edc3af-4c50-42ea-a356-e2862fe7a444")

# # ---- Static Analyzers ---- find_program(CLANGTIDY clang-tidy) if(NOT
# CLANGTIDY MATCHES CLANGTIDY-NOTFOUND) set(CMAKE_CXX_CLANG_TIDY ${CLANGTIDY})
# endif() find_program(CPPCHECK cppcheck) if(NOT CPPCHECK MATCHES
# CPPCHECK-NOTFOUND) set(CMAKE_CXX_CPPCHECK "${CPPCHECK} --language=c++
# --enable=warning,style,performance,information") endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
else()
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
endif()
# Check if mold linker is available (only on non-macOS systems)
if(NOT APPLE)
  find_program(MOLD_LINKER mold)
  if(MOLD_LINKER)
    set(CMAKE_LINKER_TYPE MOLD)
    message(STATUS "Using mold linker: ${MOLD_LINKER}")
  else()
    message(STATUS "mold linker not found, using default linker")
  endif()
else()
  message(STATUS "macOS detected, using default linker (mold not supported)")
endif()

project(MathTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_MODULE_STD 1)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)

# ---- compile_commands.json ----
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(BOOST_UT_DISABLE_MODULE OFF)
add_subdirectory("${CMAKE_SOURCE_DIR}/.." "extern_build/math")
add_subdirectory("${CMAKE_SOURCE_DIR}/../extern/ut" "extern_build/ut")

cpmaddpackage("gh:TheLartians/Format.cmake@1.7.3")

# set(FOR_TESTING ON CACHE BOOL "Build for testing")
# add_subdirectory("${PROJECT_SOURCE_DIR}/.." "extern_build/math")

# ---- Create binary ----

# Get all test files
file(GLOB test_sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Create a test executable for each source file
foreach(test_source ${test_sources})
  get_filename_component(test_name ${test_source} NAME_WE)
  add_executable(${test_name} ${test_source})
  target_link_libraries(${test_name} PRIVATE Math ut_module)
  set_target_properties(
    ${test_name}
    PROPERTIES CXX_STANDARD 23
               CXX_VISIBILITY_PRESET hidden
               VISIBILITY_INLINES_HIDDEN ON)

  # Apply the same compiler options as the main target
  if((USE_SANITIZER MATCHES "([Uu]ndefined)")
     OR (USE_SANITIZER MATCHES "([Aa]ddress);([Uu]ndefined)"))
    target_compile_options(${test_name} PRIVATE -fsanitize-trap=all)
  elseif(USE_TYPE_SANITIZER)
    target_compile_options(${test_name} PRIVATE -fsanitize-trap=all
                                                -fsanitize=type)
    target_link_options(${test_name} PRIVATE -fsanitize=type)
  endif()

  if(POLYMATHNOEXPLICITSIMDARRAY)
    target_compile_definitions(${test_name} PRIVATE POLYMATHNOEXPLICITSIMDARRAY)
  endif()

  target_compile_options(${test_name} PUBLIC -Wall -Wpedantic -Wextra -Wshadow
                                             -Werror)

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(
      ${test_name} PRIVATE -ferror-limit=8 -fcolor-diagnostics)# -ftime-trace)
    # -fexperimental-new-constant-interpreter

    # TODO: try if this passes with clang 22
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(
      ${test_name}
      PRIVATE -fmax-errors=8 -fconcepts-diagnostics-depth=4
              -fno-semantic-interposition -fdiagnostics-color=always
              -Wno-comma-subscript -Wno-psabi)
  endif()

  target_compile_options(${test_name} PRIVATE -D_GLIBCXX_ASSERTIONS
                                              -ftemplate-backtrace-limit=0)

  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# ---- Add MathTests ----

include(CTest)
enable_testing()

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # LLVM/Clang coverage instrumentation
    target_compile_options(Math PUBLIC -O0 -g -fprofile-instr-generate
                                       -fcoverage-mapping)
    target_link_options(Math PUBLIC -fprofile-instr-generate)
    add_custom_target(
      clean_coverage
      COMMAND rm -f ${PROJECT_BINARY_DIR}/reports/*.profraw
      COMMAND rm -f ${PROJECT_BINARY_DIR}/CMakeFiles/*/*.profraw
      COMMAND
        rm -f
        ${PROJECT_BINARY_DIR}/_deps/math-build/CMakeFiles/Math.dir/mod/**/*.profraw
      COMMAND rm -f ${PROJECT_BINARY_DIR}/*.profraw
      COMMENT "Clearing out old coverage data.")
  else()
    # GCC gcov coverage
    target_compile_options(Math PUBLIC -O0 -g --coverage)
    target_link_options(Math PUBLIC --coverage)
    add_custom_target(
      clean_coverage
      COMMAND rm -f ${PROJECT_BINARY_DIR}/CMakeFiles/*/test/*.gcda
      COMMAND
        rm -f
        ${PROJECT_BINARY_DIR}/_deps/math-build/CMakeFiles/Math.dir/mod/**/*.gcda
      COMMENT "Clearing out old coverage data.")
  endif()
  add_dependencies(Math clean_coverage)
  message(STATUS "Test coverage enabled.")
endif()

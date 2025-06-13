cmake_minimum_required(VERSION 4.0.2 FATAL_ERROR)
# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)
option(POLYMATHNOEXPLICITSIMDARRAY "No explicit SIMD for Array operations" OFF)
option(ENABLE_LLD "Use lld as the linker" ON)
option(USE_TYPE_SANITIZER "Use type sanitizer" OFF)

add_compile_options(
  "$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-omit-frame-pointer>")
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
  set(ENABLE_RTTI
      ON
      CACHE BOOL "enable RTTI when using ubsan")
else()
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:CXX>:-fno-rtti;-fstack-protector-strong>")
endif()
if(USE_TYPE_SANITIZER)
  set(USE_MI_MALLOC
      OFF
      CACHE BOOL "No mimalloc with sanitizers")
  set(USE_JE_MALLOC
      OFF
      CACHE BOOL "No jemalloc with sanitizers")
endif()
if(USE_MODULES)
  # set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
  # "d0edc3af-4c50-42ea-a356-e2862fe7a444")
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
      # This specific value changes as experimental support evolves. See
      # `Help/dev/experimental.rst` in the CMake source corresponding to
      # your CMake build for the exact value to use.
      "d0edc3af-4c50-42ea-a356-e2862fe7a444")
  # cmake master uses: "d0edc3af-4c50-42ea-a356-e2862fe7a444")
endif()

# # ---- Static Analyzers ---- find_program(CLANGTIDY clang-tidy) if(NOT
# CLANGTIDY MATCHES CLANGTIDY-NOTFOUND) set(CMAKE_CXX_CLANG_TIDY ${CLANGTIDY})
# endif() find_program(CPPCHECK cppcheck) if(NOT CPPCHECK MATCHES
# CPPCHECK-NOTFOUND) set(CMAKE_CXX_CPPCHECK "${CPPCHECK} --language=c++
# --enable=warning,style,performance,information") endif()

project(MathTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
if(USE_MODULES)
  set(CMAKE_CXX_MODULE_STD 1)
endif()

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)

# ---- compile_commands.json ----
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

# include(CheckCXXSymbolExists)

# # Check for _LIBCPP_VERSION check_cxx_symbol_exists(_LIBCPP_VERSION
# "<version>" LIBCXX_FOUND) if(LIBCXX_FOUND) message(STATUS "libc++ is being
# used.") FetchContent_Declare( googletest NAMES GTest URL
# https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz #
# GIT_REPOSITORY https://github.com/google/googletest.git # GIT_TAG
# 52eb8108c5bdec04579160ae17225d66034bd723 # 1.17 ) else() message(STATUS
# "libstdc++ (or another standard library) is being used.")
# FetchContent_Declare( googletest FIND_PACKAGE_ARGS NAMES GTest URL
# https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz #
# GIT_REPOSITORY https://github.com/google/googletest.git # GIT_TAG
# 52eb8108c5bdec04579160ae17225d66034bd723 # 1.17 ) endif() using
# `FIND_PACKAGE_ARGS` causes it to try `find_package`
if(DEFINED CMAKE_CXX_FLAGS) # OR
  if(${CMAKE_CXX_FLAGS} STREQUAL "")
    set(has_libcpp -1)
  else()
    string(FIND ${CMAKE_CXX_FLAGS} "-stdlib=libc++" has_libcpp)
  endif()
else()
  set(has_libcpp -1)
endif()
if(has_libcpp GREATER_EQUAL 0)
  FetchContent_Declare(
    googletest
    NAMES GTest
    URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz
    # GIT_REPOSITORY https://github.com/google/googletest.git GIT_TAG
    # 52eb8108c5bdec04579160ae17225d66034bd723 # 1.17
  )
else()
  FetchContent_Declare(
    googletest
    FIND_PACKAGE_ARGS NAMES GTest
    URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz
    # GIT_REPOSITORY https://github.com/google/googletest.git GIT_TAG
    # 52eb8108c5bdec04579160ae17225d66034bd723 # 1.17
  )
endif()
FetchContent_MakeAvailable(googletest)

cpmaddpackage("gh:TheLartians/Format.cmake@1.7.3")

set(FOR_TESTING
    ON
    CACHE BOOL "Build for testing")
add_subdirectory("${PROJECT_SOURCE_DIR}/.." "extern_build/math")

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} Math GTest::gtest_main)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 23)

if((USE_SANITIZER MATCHES "([Uu]ndefined)") OR (USE_SANITIZER MATCHES
                                                "([Aa]ddress);([Uu]ndefined)"))
  target_compile_options(${PROJECT_NAME} PRIVATE -fsanitize-trap=all)
elseif(USE_TYPE_SANITIZER)
  target_compile_options(${PROJECT_NAME} PRIVATE -fsanitize-trap=all
                                                 -fsanitize=type)
  target_link_options(${PROJECT_NAME} PRIVATE -fsanitize=type)
endif()
if(POLYMATHNOEXPLICITSIMDARRAY)
  target_compile_definitions(${PROJECT_NAME}
                             PRIVATE POLYMATHNOEXPLICITSIMDARRAY)
endif()
if(NOT USE_MODULES)
  set_target_properties(${PROJECT_NAME} PROPERTIES CXX_SCAN_FOR_MODULES 0)
endif()

# enable compiler warnings
if(NOT TEST_INSTALLED_VERSION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES
                                              "GNU")
    # target_compile_options(Math PUBLIC -Wall -Wpedantic -Wextra -Wshadow
    # -Werror)
    target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wpedantic -Wextra
                                                  -Wshadow -Werror)
  elseif(MSVC)
    # target_compile_options(Math PUBLIC /W4 /WX)
    target_compile_definitions(${PROJECT_NAME}
                               PUBLIC DOCTEST_CONFIG_USE_STD_HEADERS)
  endif()
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(
    ${PROJECT_NAME} PRIVATE -ferror-limit=8 -fcolor-diagnostics -ftime-trace)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(
    ${PROJECT_NAME}
    PRIVATE -fmax-errors=8 -fconcepts-diagnostics-depth=4
            -fno-semantic-interposition -fdiagnostics-color=always
            -Wno-comma-subscript -Wno-psabi)
endif()
# target_compile_options(Math INTERFACE -D_GLIBCXX_ASSERTIONS)
target_compile_options(${PROJECT_NAME} PRIVATE -D_GLIBCXX_ASSERTIONS
                                               -ftemplate-backtrace-limit=0)
if(ENABLE_LLD)
  # TODO: use `CMAKE_LINKER_TYPE`/`LINKER_TYPE` once you've upgraded to cmake 29
  target_link_options(${PROJECT_NAME} PRIVATE -fuse-ld=lld)
endif()

# ---- Add MathTests ----

enable_testing()

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  # target_compile_options(Math PUBLIC -O0 -g --coverage)
  # target_link_options(Math PUBLIC --coverage)
  target_compile_options(${PROJECT_NAME} PUBLIC -O0 -g --coverage)
  target_link_options(${PROJECT_NAME} PUBLIC --coverage)
  add_custom_target(
    clean_coverage
    COMMAND rm -f ${PROJECT_BINARY_DIR}/CMakeFiles/MathTests.dir/*.gcda
    COMMAND
      rm -f
      ${PROJECT_BINARY_DIR}/_deps/math-build/CMakeFiles/Math.dir/lib/*.gcda
    COMMENT "Clearing out old coverage data.")
  add_dependencies(${PROJECT_NAME} clean_coverage)
  message(STATUS "Test coverage enabled.")
endif()

cmake_minimum_required(VERSION 4.0.2 FATAL_ERROR)

add_compile_options(
  "$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fomit-frame-pointer;-fno-rtti;-fno-plt;-fno-asynchronous-unwind-tables>"
)

set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
    # This specific value changes as experimental support evolves. See
    # `Help/dev/experimental.rst` in the CMake source corresponding to
    # your CMake build for the exact value to use.
    "d0edc3af-4c50-42ea-a356-e2862fe7a444")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
else()
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
endif()
project(MathBenchmarks LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_MODULE_STD 1)

option(ENABLE_NATIVE_COMPILATION "Compile with -march=native" ON)
option(ENABLE_WIDE_VECTORS "Compile with 512bit vectors if available" ON)
# option(ENABLE_OPENMP "Use OpenMP for a multithreading benchmark" OFF)
option(POLYMATHNOEXPLICITSIMDARRAY "No explicit SIMD for Array operations" OFF)

# --- Import tools ----

# include(../cmake/tools.cmake)

# ---- Dependencies ----

include(FetchContent)

# ---- compile_commands.json ----
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# if("${CMAKE_CXX_FLAGS}" STREQUAL "") set(has_libcpp -1) else() string(FIND
# ${CMAKE_CXX_FLAGS} "-stdlib=libc++" has_libcpp) endif() if(has_libcpp
# GREATER_EQUAL 0) main doesn't support -fno-exceptions or -fno-rtti
FetchContent_Declare(
  nanobench
  GIT_REPOSITORY https://github.com/Spartan322/nanobench
  GIT_TAG f49b052a68bc05d5497acf1f6c0c2d6ef65c6eef
  # GIT_REPOSITORY https://github.com/martinus/nanobench.git GIT_TAG v4.3.11
  GIT_SHALLOW TRUE)
# else() FetchContent_Declare( nanobench FIND_PACKAGE_ARGS GIT_REPOSITORY
# https://github.com/Spartan322/nanobench GIT_REPOSITORY
# https://github.com/martinus/nanobench.git GIT_TAG v4.3.11 GIT_SHALLOW TRUE)
# endif()

FetchContent_MakeAvailable(nanobench)

# FetchContent_Declare( Math GIT_REPOSITORY git@github.com:LoopModels/Math.git
# GIT_TAG origin/main ) FetchContent_MakeAvailable(Math)

# file(GLOB_RECURSE headers CONFIGURE_DEPENDS
# ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB benchmarks CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE module_interfaces CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cxxm)
file(GLOB_RECURSE module_implementations CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)

add_executable(${PROJECT_NAME} ${benchmarks})
target_sources(
  ${PROJECT_NAME}
  PUBLIC FILE_SET CXX_MODULES FILES ${module_interfaces}
  PRIVATE ${module_implementations})
# INTERFACE FILE_SET HEADERS FILES ${headers})

message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")
message(STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
add_subdirectory("${PROJECT_SOURCE_DIR}/.." "extern_build/math")
# if(ENABLE_OPENMP) find_package(OpenMP) target_link_libraries(${PROJECT_NAME}
# OpenMP::OpenMP_CXX) endif()

# --- Maybe Use SIMD Array Ops ---
if(POLYMATHNOEXPLICITSIMDARRAY)
  target_compile_definitions(${PROJECT_NAME} PUBLIC POLYMATHNOEXPLICITSIMDARRAY)
  message(STATUS "Do not use SIMD Array operations")
else()
  message(STATUS "Use SIMD Array operations")
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE Math nanobench)

if((CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES
                                               "IntelLLVM"))
  target_compile_options(${PROJECT_NAME} PRIVATE -ferror-limit=2
                                                 -fcolor-diagnostics)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(
    ${PROJECT_NAME}
    PRIVATE -fmax-errors=2 -fconcepts-diagnostics-depth=4
            -fno-semantic-interposition -fdiagnostics-color=always
            -fverbose-asm)
endif()

if(ENABLE_NATIVE_COMPILATION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    target_compile_options(${PROJECT_NAME} PRIVATE -xhost)
    if(ENABLE_WIDE_VECTORS)
      target_compile_options(${PROJECT_NAME} PRIVATE -qopt-zmm-usage=high)
    endif()
  else()
    target_compile_options(${PROJECT_NAME} PRIVATE -march=native
    )# -fno-unroll-loops)
    if(ENABLE_WIDE_VECTORS)
      check_cxx_compiler_flag("-mprefer-vector-width=512" VEC512)
      if(VEC512)
        target_compile_options(${PROJECT_NAME}
                               PRIVATE -mprefer-vector-width=512)
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
          if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(${PROJECT_NAME}
                                   PRIVATE -mtune-ctrl=avx512_move_by_pieces)
          endif()
        endif()
      endif()
    endif()
  endif()
endif()
set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CXX_STANDARD 23
             CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)
set_target_properties(
  ${PROJECT_NAME} PROPERTIES ENVIRONMENT
                             WORKING_DIRECTORY=${PROJECT_BINARY_DIR})

target_compile_options(
  ${PROJECT_NAME}
  PRIVATE -fstrict-aliasing
          -fstrict-overflow
          -fno-signed-zeros
          -fassociative-math
          -ffinite-math-only
          -funsafe-math-optimizations
          -fno-trapping-math
          -Wall
          -Wshadow
          -Wextra
          # -save-temps
          -Werror)
# if(ENABLE_OPENMP) if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
# target_compile_options(${PROJECT_NAME} PRIVATE -fiopenmp) else()
# target_compile_options(${PROJECT_NAME} PRIVATE -fopenmp) endif() else()
# target_compile_options(${PROJECT_NAME} PRIVATE -fopenmp-simd) endif()
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
  target_compile_options(${PROJECT_NAME} PRIVATE -masm=intel)
endif()

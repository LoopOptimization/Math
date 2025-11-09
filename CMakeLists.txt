cmake_minimum_required(VERSION 3.31 FATAL_ERROR)
cmake_policy(SET CMP0167 NEW) # new find boost

option(ENABLE_LLD "Use LLD? Off by default." OFF)
option(USE_MI_MALLOC "Use MI_MALLOC if available." ON)
option(USE_JE_MALLOC "Use jemallc if available and not using mi_malloc." ON)
option(ENABLE_WIDE_VECTORS "Compile with 512bit vectors if available" ON)
option(POLYMATHNOEXPLICITSIMDARRAY "No explicit SIMD for Array operations" OFF)
# option(FOR_TESTING "Build for testing" OFF)

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-rtti>")
set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
    # This specific value changes as experimental support evolves. See
    # `Help/dev/experimental.rst` in the CMake source corresponding to
    # your CMake build for the exact value to use.
    "d0edc3af-4c50-42ea-a356-e2862fe7a444")

# ---- Project ----

message(
  STATUS
    "CMAKE_EXPERIMENTAL_CXX_IMPORT_STD = ${CMAKE_EXPERIMENTAL_CXX_IMPORT_STD}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
else()
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
endif()

# Note: update this to your new project's name and version
project(
  Math
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_MODULE_STD ON)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ---- Add dependencies via CPM ----
# see https://github.com/TheLartians/CPM.cmake for more info

include(cmake/CPM.cmake)

# PackageProject.cmake will be used to make our target installable
cpmaddpackage("gh:TheLartians/PackageProject.cmake@1.8.0")

# ---- Add source files ----

# Note: globbing sources is considered bad practice as CMake's generators may
# not detect new files automatically. Keep that in mind when changing files, or
# explicitly mention them here.
file(GLOB_RECURSE headers CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx")
file(GLOB_RECURSE module_interfaces CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/mod/*.cxxm")
file(GLOB_RECURSE module_implementations CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/mod/*.cxx")

# ---- Create library ----
add_library(${PROJECT_NAME} OBJECT) # ${headers})
target_sources(
  ${PROJECT_NAME}
  PUBLIC FILE_SET CXX_MODULES FILES ${module_interfaces}
  PRIVATE ${module_implementations}
  INTERFACE FILE_SET HEADERS FILES ${headers})
target_include_directories(
  ${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>)
# --- Maybe Use SIMD Array Ops ---
if(POLYMATHNOEXPLICITSIMDARRAY)
  target_compile_definitions(${PROJECT_NAME}
                             PRIVATE POLYMATHNOEXPLICITSIMDARRAY)
endif()
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_23)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(
    ${PROJECT_NAME} PRIVATE -ferror-limit=8 -fcolor-diagnostics) # -ftime-trace)
  # -fexperimental-new-constant-interpreter)

  # TODO: try if this passes with clang 22
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(
    ${PROJECT_NAME}
    PRIVATE -fmax-errors=8 -fconcepts-diagnostics-depth=4
            -fno-semantic-interposition -fdiagnostics-color=always
            -Wno-comma-subscript -Wno-psabi)
endif()
target_compile_options(${PROJECT_NAME} PRIVATE -Wshadow -Werror) # -Wpedantic
                                                                 # -Wextra -Wall
# -Wno-reference-tu-local-entity-in-other-tu)
target_compile_features(
  ${PROJECT_NAME}
  PRIVATE cxx_std_23
  INTERFACE cxx_std_23)
set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CXX_STANDARD 23
             CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 23)

# being a cross-platform target, we enforce standards conformance on MSVC
# target_compile_options( ${PROJECT_NAME} INTERFACE
# "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")

# --- Allocator override ----
set(SEARCH_JE_MALLOC USE_JE_MALLOC)
if(USE_MI_MALLOC)
  find_package(mimalloc) # REQUIRED)
  if(mimalloc_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE mimalloc) # mimalloc-static
    target_compile_definitions(${PROJECT_NAME} PRIVATE USING_MIMALLOC)
    set(SEARCH_JE_MALLOC OFF)
  endif()
endif()
if(SEARCH_JE_MALLOC)
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(JEMALLOC jemalloc)
    if(jemalloc_FOUND)
      target_link_libraries(${PROJECT_NAME} PRIVATE ${JEMALLOC_LIBRARIES})
      target_include_directories(${PROJECT_NAME}
                                 PRIVATE ${JEMALLOC_INCLUDE_DIRS})
      target_compile_definitions(${PROJECT_NAME} PRIVATE USING_JEMALLOC)
    endif()
  endif()
endif()

# --- Boost Unordered ---

target_compile_definitions(${PROJECT_NAME} PUBLIC BOOST_EXCEPTION_DISABLE)

find_package(Boost 1.83.0 COMPONENTS headers)
if(Boost_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE Boost::headers)
else()
  include(FetchContent)
  FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.89.0/boost-1.89.0-cmake.tar.xz
        EXCLUDE_FROM_ALL)
  FetchContent_MakeAvailable(Boost)
  target_link_libraries(${PROJECT_NAME} PRIVATE Boost::unordered)
endif()
# Link Boost

# ---- Create an installable target ----
# this allows users to install and find the library via `find_package()`.

# the location where the project's version header will be placed should match
# the project's regular header paths string(TOLOWER ${PROJECT_NAME}/version.h
# VERSION_HEADER_LOCATION)

include(GNUInstallDirs)

# install( TARGETS ${PROJECT_NAME} EXPORT math-targets FILE_SET CXX_MODULES #
# There's currently no convention for this location DESTINATION
# ${CMAKE_INSTALL_LIBDIR}/cmake/math/mod FILE_SET HEADERS DESTINATION
# ${CMAKE_INSTALL_INCLUDEDIR} # Same as default, could be omitted INCLUDES
# DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}) install( EXPORT math-targets
# DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/math CXX_MODULES_DIRECTORY .)

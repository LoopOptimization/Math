cmake_minimum_required(VERSION 3.28.2)

option(ENABLE_LLD "Use LLD? Off by default." OFF)
option(USE_MI_MALLOC "Use MI_MALLOC if available." ON)
option(USE_JE_MALLOC "Use jemallc if available and not using mi_malloc." ON)
option(ENABLE_NATIVE_COMPILATION "Compile with -march=native" ON)
option(ENABLE_WIDE_VECTORS "Compile with 512bit vectors if available" ON)
option(ENABLE_RTTI "Use RTTI" OFF)
option(ENABLE_STACK_PROTECTOR "Use stack protection" OFF)
# ---- Project ----

# Note: update this to your new project's name and version
project(
  Math
  VERSION 1.0
  LANGUAGES CXX)

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
file(GLOB_RECURSE modules CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/mod/*.cxx")
message(STATUS "Headers: ${headers}")
message(STATUS "Modules: ${modules}")

# ---- Create library ----

# Note: for header-only libraries change all PUBLIC flags to INTERFACE and
# create an interface target:
add_library(${PROJECT_NAME} OBJECT) # ${headers})
# add_library(${PROJECT_NAME} INTERFACE ${headers})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 23)
target_sources(
  ${PROJECT_NAME}
  PUBLIC FILE_SET CXX_MODULES FILES ${modules}
  INTERFACE FILE_SET HEADERS FILES ${headers})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_23)

# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(
  ${PROJECT_NAME} INTERFACE "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")

target_include_directories(
  ${PROJECT_NAME}
  PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
          $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>)
target_compile_options(${PROJECT_NAME} PRIVATE -fno-exceptions)
if(NOT ENABLE_RTTI)
  target_compile_options(${PROJECT_NAME} PRIVATE -fno-rtti)
endif()
if(ENABLE_STACK_PROTECTOR)
  target_compile_options(${PROJECT_NAME} PRIVATE -fstack-protector-strong)
endif()

# --- Allocator override ----
set(SEARCH_JE_MALLOC USE_JE_MALLOC)
if(USE_MI_MALLOC)
  find_package(mimalloc) # REQUIRED)
  if(mimalloc_FOUND)
    target_link_libraries(${PROJECT_NAME} INTERFACE mimalloc) # mimalloc-static
    target_compile_definitions(${PROJECT_NAME} INTERFACE USING_MIMALLOC)
    set(SEARCH_JE_MALLOC OFF)
  endif()
endif()
if(SEARCH_JE_MALLOC)
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(JEMALLOC jemalloc)
    if(jemalloc_FOUND)
      target_link_libraries(${PROJECT_NAME} INTERFACE ${JEMALLOC_LIBRARIES})
      target_include_directories(${PROJECT_NAME}
                                 INTERFACE ${JEMALLOC_INCLUDE_DIRS})
      target_compile_definitions(${PROJECT_NAME} INTERFACE USING_JEMALLOC)
    endif()
  endif()
endif()

# --- Native compilation ---

if(ENABLE_NATIVE_COMPILATION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    target_compile_options(${PROJECT_NAME} INTERFACE -xhost)
    if(ENABLE_WIDE_VECTORS)
      target_compile_options(${PROJECT_NAME} INTERFACE -qopt-zmm-usage=high)
    endif()
  else()
    target_compile_options(${PROJECT_NAME} INTERFACE -march=native)
    if(ENABLE_WIDE_VECTORS)
      include(CheckCXXCompilerFlag)
      check_cxx_compiler_flag("-mprefer-vector-width=512" VEC512)
      if(VEC512)
        if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
          target_compile_options(${PROJECT_NAME}
                                 INTERFACE -mprefer-vector-width=512)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
          target_compile_options(
            ${PROJECT_NAME} INTERFACE -mprefer-vector-width=512
                                      -mtune-ctrl=avx512_move_by_pieces)
        endif()
      endif()
    endif()
  endif()
endif()
# ---- Create an installable target ----
# this allows users to install and find the library via `find_package()`.

# the location where the project's version header will be placed should match
# the project's regular header paths
string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)

include(GNUInstallDirs)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT math-targets
  FILE_SET CXX_MODULES
  # There's currently no convention for this location
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/my_package/src
  FILE_SET HEADERS
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Same as default, could be omitted
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(
  EXPORT math-targets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/math
  CXX_MODULES_DIRECTORY .)
# packageproject( NAME ${PROJECT_NAME} VERSION ${PROJECT_VERSION} NAMESPACE
# ${PROJECT_NAME} BINARY_DIR ${PROJECT_BINARY_DIR} INCLUDE_DIR
# ${PROJECT_SOURCE_DIR}/include INCLUDE_DESTINATION
# include/${PROJECT_NAME}-${PROJECT_VERSION} VERSION_HEADER
# "${VERSION_HEADER_LOCATION}" COMPATIBILITY SameMajorVersion)

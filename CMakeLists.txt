cmake_minimum_required(VERSION 4.0.2 FATAL_ERROR)

option(ENABLE_LLD "Use LLD? Off by default." OFF)
option(USE_MI_MALLOC "Use MI_MALLOC if available." ON)
option(USE_JE_MALLOC "Use jemallc if available and not using mi_malloc." ON)
option(ENABLE_RTTI "Use RTTI" OFF)
option(ENABLE_STACK_PROTECTOR "Use stack protection" OFF)
option(ENABLE_NATIVE_COMPILATION "Compile with -march=native" ON)
option(ENABLE_WIDE_VECTORS "Compile with 512bit vectors if available" ON)
option(POLYMATHNOEXPLICITSIMDARRAY "No explicit SIMD for Array operations" OFF)
option(USE_MODULES "Use C++ Modules" OFF)
option(FOR_TESTING "Build for testing" OFF)

if(ENABLE_STACK_PROTECTOR)
    if(ENABLE_RTTI)
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-rtti;-fstack-protector-strong>")
    else()
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-rtti;-fstack-protector-strong;-fno-rtti>")
    endif()
else()
    if(ENABLE_RTTI)
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-rtti>")
    else()
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fno-exceptions;-fno-rtti;-fstack-protector-strong;-fno-rtti>")
    endif()
endif()
if (USE_MODULES)
    # set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
    #     "d0edc3af-4c50-42ea-a356-e2862fe7a444")
    set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
      # This specific value changes as experimental support evolves. See
      # `Help/dev/experimental.rst` in the CMake source corresponding to
      # your CMake build for the exact value to use.
        "a9e1cf81-9932-4810-974b-6eccaf14e457")
      # cmake master uses:
        # "d0edc3af-4c50-42ea-a356-e2862fe7a444")
endif()
# ---- Project ----

message(STATUS "CMAKE_EXPERIMENTAL_CXX_IMPORT_STD = ${CMAKE_EXPERIMENTAL_CXX_IMPORT_STD}")
# Note: update this to your new project's name and version
project(
  Math
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
if (USE_MODULES)
    set(CMAKE_CXX_MODULE_STD 1)
else()
    set_target_properties(${PROJECT_NAME} PROPERTIES CXX_SCAN_FOR_MODULES 0)
endif()

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

if(USE_MODULES)
  # Note: globbing sources is considered bad practice as CMake's generators may
  # not detect new files automatically. Keep that in mind when changing files,
  # or explicitly mention them here.
  file(GLOB_RECURSE headers CONFIGURE_DEPENDS
       "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx")
  file(GLOB_RECURSE modules CONFIGURE_DEPENDS
       "${CMAKE_CURRENT_SOURCE_DIR}/mod/*.cxx")

    # ---- Create library ----
  add_library(${PROJECT_NAME} OBJECT) # ${headers})
  target_sources(
    ${PROJECT_NAME}
    PUBLIC FILE_SET CXX_MODULES FILES ${modules}
    INTERFACE FILE_SET HEADERS FILES ${headers})
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_MODULE)
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
  target_compile_definitions(${PROJECT_NAME} PRIVATE MATHTESTSTLDEFINITIONS)

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
else()
  file(GLOB_RECURSE headers CONFIGURE_DEPENDS
       "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx"
       "${CMAKE_CURRENT_SOURCE_DIR}/mod/*.cxx")

    # ---- Create library ----
  add_library(${PROJECT_NAME} INTERFACE ${headers})
  target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_23)
  target_include_directories(
    ${PROJECT_NAME}
    INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
              $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>)
  target_include_directories(
    ${PROJECT_NAME}
    INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/mod>
              $<INSTALL_INTERFACE:mod/${PROJECT_NAME}-${PROJECT_VERSION}>)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 23)
target_compile_features(${PROJECT_NAME}
  PRIVATE   cxx_std_23
  INTERFACE cxx_std_23)



# being a cross-platform target, we enforce standards conformance on MSVC
# target_compile_options( ${PROJECT_NAME} INTERFACE
# "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")

# --- Allocator override ----
set(SEARCH_JE_MALLOC USE_JE_MALLOC)
if(USE_MI_MALLOC)
  find_package(mimalloc) # REQUIRED)
  if(mimalloc_FOUND)
    if(USE_MODULES)
      target_link_libraries(${PROJECT_NAME} PRIVATE mimalloc) # mimalloc-static
      target_compile_definitions(${PROJECT_NAME} PRIVATE USING_MIMALLOC)
    else()
      target_link_libraries(${PROJECT_NAME} INTERFACE mimalloc
      )# mimalloc-static
      target_compile_definitions(${PROJECT_NAME} INTERFACE USING_MIMALLOC)
    endif()
    set(SEARCH_JE_MALLOC OFF)
  endif()
endif()
if(SEARCH_JE_MALLOC)
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(JEMALLOC jemalloc)
    if(jemalloc_FOUND)
      if(USE_MODULES)
        target_link_libraries(${PROJECT_NAME} PRIVATE ${JEMALLOC_LIBRARIES})
        target_include_directories(${PROJECT_NAME}
                                   PRIVATE ${JEMALLOC_INCLUDE_DIRS})
        target_compile_definitions(${PROJECT_NAME} PRIVATE USING_JEMALLOC)
      else()
        target_link_libraries(${PROJECT_NAME} INTERFACE ${JEMALLOC_LIBRARIES})
        target_include_directories(${PROJECT_NAME}
                                   INTERFACE ${JEMALLOC_INCLUDE_DIRS})
        target_compile_definitions(${PROJECT_NAME} INTERFACE USING_JEMALLOC)
      endif()
    endif()
  endif()
endif()

# --- Native compilation ---

if(ENABLE_NATIVE_COMPILATION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    if(USE_MODULES)
      target_compile_options(${PROJECT_NAME} PUBLIC -xhost)
      if(ENABLE_WIDE_VECTORS)
        target_compile_options(${PROJECT_NAME} PUBLIC -qopt-zmm-usage=high)
      endif()
    else()
      target_compile_options(${PROJECT_NAME} INTERFACE -xhost)
      if(ENABLE_WIDE_VECTORS)
        target_compile_options(${PROJECT_NAME} INTERFACE -qopt-zmm-usage=high)
      endif()
    endif()
  else()
    if(USE_MODULES)
      target_compile_options(${PROJECT_NAME} PUBLIC -march=native)
    else()
      target_compile_options(${PROJECT_NAME} INTERFACE -march=native)
    endif()
    if(ENABLE_WIDE_VECTORS)
      include(CheckCXXCompilerFlag)
      check_cxx_compiler_flag("-mprefer-vector-width=512" VEC512)
      if(VEC512)
        if(USE_MODULES)
          if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(${PROJECT_NAME}
                                   PUBLIC -mprefer-vector-width=512)
          elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            target_compile_options(
              ${PROJECT_NAME} PUBLIC -mprefer-vector-width=512
                                     -mtune-ctrl=avx512_move_by_pieces)
          endif()
        else()
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
endif()
# ---- Create an installable target ----
# this allows users to install and find the library via `find_package()`.

# the location where the project's version header will be placed should match
# the project's regular header paths string(TOLOWER ${PROJECT_NAME}/version.h
# VERSION_HEADER_LOCATION)

include(GNUInstallDirs)

if(USE_MODULES)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT math-targets
    FILE_SET CXX_MODULES
    # There's currently no convention for this location
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/math/mod
    FILE_SET HEADERS
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Same as default, could be omitted
    INCLUDES
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(
    EXPORT math-targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/math
    CXX_MODULES_DIRECTORY .)
else()
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT math-targets
    FILE_SET HEADERS
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # Same as default, could be omitted
    INCLUDES
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(EXPORT math-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/math)
endif()

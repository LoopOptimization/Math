#pragma once
#ifndef Macros_hxx_INCLUDED
#define Macros_hxx_INCLUDED

#ifdef __clang__

// #define NODEBUG
// #define TRIVIAL
#define NODEBUG [[gnu::nodebug]]
#define TRIVIAL [[gnu::nodebug, gnu::always_inline, gnu::flatten]]
// clang-format off
#define DEBUGABLE_START                                                       \
  _Pragma("clang attribute DoNotOptimize.push (__attribute((optnone)), apply_to = function)")
// clang-format on
#define DEBUGABLE_END _Pragma("clang attribute DoNotOptimize.pop")
#define DEBUGABLE [[clang::optnone]]

#elifdef __GNUC__

#define NODEBUG [[gnu::optimize("O3")]]
#define TRIVIAL [[gnu::always_inline, gnu::flatten, gnu::optimize("O3")]]
// clang-format off
#define DEBUGABLE_START                                                       \
  _Pragma("gcc attribute DoNotOptimize.push (__attribute((optimize(\"O0\"))), apply_to = function)")
// clang-format on
#define DEBUGABLE_END _Pragma("gcc attribute DoNotOptimize.pop")
#define DEBUGABLE [[gnu::optimize("O0")]]

#else

#define OPTIMIZE_TRIVIAL
#define DEBUGABLE_START
#define DEBUGABLE_END
#define DEBUGABLE

#endif

#ifndef NDEBUG
#define DEBUGUSED [[gnu::used]]
#else
#define DEBUGUSED
#endif

#endif // Macros_hxx_INCLUDED

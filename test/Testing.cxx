module;
#include <cstdio>
#include <cstdlib>
module Testing;
import CorePrint;
import std;

namespace testing {

namespace {
int tests_pass{};
int tests_fail{};
int asserts_pass{};
int asserts_fail{};
bool current_test_failed{};
bool atexit_registered{};
const char *fatal_file_{};
unsigned fatal_line_{};

void print_summary() {
  ::utils::print('\n');
  if (tests_fail) {
    ::utils::print("\033[31m"); // red
    ::utils::print("FAILED: ");
  } else {
    ::utils::print("\033[32m"); // green
    ::utils::print("ALL PASSED: ");
  }
  ::utils::print(tests_pass, " test(s) passed, ", tests_fail, " test(s) failed");
  ::utils::print(" (", asserts_pass, " assertion(s) passed, ", asserts_fail,
                 " assertion(s) failed)");
  ::utils::print("\033[0m\n");
  ::utils::flush();
  if (tests_fail) std::_Exit(1);
}
} // namespace

void report_pass() { ++asserts_pass; }

void report_fail(const char *file, unsigned line) {
  ++asserts_fail;
  current_test_failed = true;
  ::utils::print("\033[31m  FAIL\033[0m ");
  ::utils::print(std::string_view{file});
  ::utils::print(':');
  ::utils::print(int(line));
  ::utils::print('\n');
}

void set_fatal_location(const char *file, unsigned line) {
  fatal_file_ = file;
  fatal_line_ = line;
}

[[noreturn]] void do_fatal_abort() {
  if (fatal_file_) {
    ::utils::print("\033[31m  FATAL\033[0m ");
    ::utils::print(std::string_view{fatal_file_});
    ::utils::print(':');
    ::utils::print(int(fatal_line_));
    ::utils::print('\n');
  }
  print_summary();
  __builtin_trap();
}

auto should_run_test(std::string_view name) -> bool {
  auto f = cfg<>.filter;
  if (f.empty()) return true;
  return name.find(f) != std::string_view::npos;
}

void begin_test(std::string_view name) {
  if (!atexit_registered) {
    std::atexit(print_summary);
    atexit_registered = true;
  }
  current_test_failed = false;
  ::utils::print("Running \"");
  ::utils::print(name);
  ::utils::print("\"...");
  ::utils::flush();
}

void end_test(std::string_view name) {
  if (current_test_failed) {
    ::utils::print(" \033[31mFAILED\033[0m\n");
    ++tests_fail;
  } else {
    ::utils::print(" \033[32mPASSED\033[0m\n");
    ++tests_pass;
  }
  ::utils::flush();
}

} // namespace testing

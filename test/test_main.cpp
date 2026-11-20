#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "vendor/doctest.h"

#ifdef BLUEIMU_PIO_TEST_REPORT
// PlatformIO's doctest runner only counts test cases that appear in the console
// output, and passing cases print nothing by default (`--success` prints every
// assertion, which is very slow to parse). This listener prints one block per
// passing test case in the console reporter's format; failures are already
// reported by the console reporter itself.
#include <iostream>

namespace
{
   struct PioPassReporter : public doctest::IReporter
   {
      const doctest::TestCaseData* current = nullptr;

      explicit PioPassReporter(const doctest::ContextOptions&) {}

      void test_case_start(const doctest::TestCaseData& in) override { current = &in; }

      void test_case_end(const doctest::CurrentTestCaseStats& st) override
      {
         if (!current || !st.testCaseSuccess) return;
         std::cout << "===============================================================================\n"
                   << current->m_file.c_str() << ":" << current->m_line << ":\n";
         if (current->m_test_suite && current->m_test_suite[0])
            std::cout << "TEST SUITE: " << current->m_test_suite << "\n";
         std::cout << "TEST CASE:  " << current->m_name << "\n\n";
      }

      void report_query(const doctest::QueryData&) override {}
      void test_run_start() override {}
      void test_run_end(const doctest::TestRunStats&) override {}
      void test_case_reenter(const doctest::TestCaseData&) override {}
      void test_case_exception(const doctest::TestCaseException&) override {}
      void subcase_start(const doctest::SubcaseSignature&) override {}
      void subcase_end() override {}
      void log_assert(const doctest::AssertData&) override {}
      void log_message(const doctest::MessageData&) override {}
      void test_case_skipped(const doctest::TestCaseData&) override {}
   };
}

REGISTER_LISTENER("pio_pass", 1, PioPassReporter);
#endif

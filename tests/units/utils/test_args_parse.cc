#include <gtest/gtest.h>

#include "args_parse.h"

namespace utils {
namespace test {

class ArgParseTest : public ::testing::Test {
 protected:
  static char** begin(std::vector<const char*>& v) { return const_cast<char**>(v.data()); }
  static char** end(std::vector<const char*>& v) { return const_cast<char**>(v.data() + v.size()); }
};

TEST_F(ArgParseTest, GetCmdOption_Found) {
  std::vector<const char*> argv = {"prog", "--foo", "bar"};
  EXPECT_EQ(getCmdOption(begin(argv), end(argv), "--foo"), "bar");
}

TEST_F(ArgParseTest, GetCmdOption_NotFound) {
  std::vector<const char*> argv = {"prog", "--foo", "bar"};
  EXPECT_EQ(getCmdOption(begin(argv), end(argv), "--baz"), "");
}

TEST_F(ArgParseTest, GetCmdOption_OptionAtEnd) {
  std::vector<const char*> argv = {"prog", "--foo"};
  EXPECT_EQ(getCmdOption(begin(argv), end(argv), "--foo"), "");
}

TEST_F(ArgParseTest, CmdOptionExists_Present) {
  std::vector<const char*> argv = {"prog", "--verbose"};
  EXPECT_TRUE(cmdOptionExists(begin(argv), end(argv), "--verbose"));
}

TEST_F(ArgParseTest, CmdOptionExists_Absent) {
  std::vector<const char*> argv = {"prog", "--verbose"};
  EXPECT_FALSE(cmdOptionExists(begin(argv), end(argv), "--quiet"));
}

TEST_F(ArgParseTest, CmdOptionExists_EmptyList) {
  std::vector<const char*> argv;
  EXPECT_FALSE(cmdOptionExists(begin(argv), end(argv), "--verbose"));
}

}  // namespace test
}  // namespace utils

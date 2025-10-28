#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"

namespace utils
{
namespace test
{

// Wrapper class to test FIND_MAX macro
// since in GoogleTest, the test body is inside a class,
// but the class is non-copyable, so [=, *this] fails.
class FindMax1DHelper
{
 public:
  float find_max(const VECTOR_REAL_VIEW& array, int N) const
  {
    float result;
    FIND_MAX_1D(array, N, result);
    return result;
  }
};

class CommonMacrosTest : public ::testing::Test
{
 protected:
  FindMax1DHelper helper;
};

TEST_F(CommonMacrosTest, FindMax1D)
{
  constexpr int N = 5;
  auto array = allocateVector<VECTOR_REAL_VIEW>(N, "array");
  array(1) = 1.0f;
  array(2) = 3.5f;
  array(3) = -2.0f;
  array(4) = 7.2f;
  array(0) = 4.1f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 7.2f);
}

TEST_F(CommonMacrosTest, FindMax1DLargeArray)
{
  constexpr int N = 2000;
  auto array = allocateVector<VECTOR_REAL_VIEW>(N, "array_large");
  for (int i = 0; i < N; ++i)
  {
    array(i) = 3000.0f;
  }
  array(N / 2) = 3300.0f;  // One element is different
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 3300.0f);
}

TEST_F(CommonMacrosTest, FindMax1DAllNegative)
{
  constexpr int N = 4;
  auto array = allocateVector<VECTOR_REAL_VIEW>(N, "array_neg");
  array(0) = -5.0f;
  array(1) = -3.2f;
  array(2) = -7.8f;
  array(3) = -1.1f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, -1.1f);
}

TEST_F(CommonMacrosTest, FindMax1DAllEqual)
{
  constexpr int N = 3;
  auto array = allocateVector<VECTOR_REAL_VIEW>(N, "array_eq");
  array(0) = 2.5f;
  array(1) = 2.5f;
  array(2) = 2.5f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 2.5f);
}

TEST_F(CommonMacrosTest, FindMax1DSingleElement)
{
  constexpr int N = 1;
  auto array = allocateVector<VECTOR_REAL_VIEW>(N, "array_single");
  array(0) = 42.0f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 42.0f);
}

TEST_F(CommonMacrosTest, FindMax1DEmptyArrayThrows)
{
  constexpr int N = 0;
  auto array = allocateVector<VECTOR_REAL_VIEW>(N, "array_empty");
  EXPECT_THROW(helper.find_max(array, N), std::runtime_error);
  EXPECT_THROW(helper.find_max(array, 10), std::runtime_error);
}

}  // namespace test
}  // namespace utils
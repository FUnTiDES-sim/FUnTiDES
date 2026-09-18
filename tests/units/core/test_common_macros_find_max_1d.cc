#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"

namespace utils {
namespace test {

// Wrapper class to test FIND_MAX/FIND_MIN/SUM macros.
// GoogleTest test body lives inside a class, so [=, *this] fails — helper works around it.
class FindMax1DHelper {
 public:
  float find_max(const vectorReal& array, int N) const {
    float result;
    FIND_MAX_1D(array, N, result);
    return result;
  }

  float find_min(const vectorReal& array, int N) const {
    float result = array(0);
    FIND_MIN(array, N, result);
    return result;
  }

  float sum(const vectorReal& array, int N) const {
    float result = 0.0f;
    SUM(array, N, result);
    return result;
  }
};

class CommonMacrosTest : public ::testing::Test {
 protected:
  FindMax1DHelper helper;
};

TEST_F(CommonMacrosTest, FindMax1D) {
  constexpr int N = 5;
  auto array = allocateVector<vectorReal>(N, "array");
  array(1) = 1.0f;
  array(2) = 3.5f;
  array(3) = -2.0f;
  array(4) = 7.2f;
  array(0) = 4.1f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 7.2f);
}

TEST_F(CommonMacrosTest, FindMax1DLargeArray) {
  constexpr int N = 2000;
  auto array = allocateVector<vectorReal>(N, "array_large");
  for (int i = 0; i < N; ++i) {
    array(i) = 3000.0f;
  }
  array(N / 2) = 3300.0f;  // One element is different
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 3300.0f);
}

TEST_F(CommonMacrosTest, FindMax1DAllNegative) {
  constexpr int N = 4;
  auto array = allocateVector<vectorReal>(N, "array_neg");
  array(0) = -5.0f;
  array(1) = -3.2f;
  array(2) = -7.8f;
  array(3) = -1.1f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, -1.1f);
}

TEST_F(CommonMacrosTest, FindMax1DAllEqual) {
  constexpr int N = 3;
  auto array = allocateVector<vectorReal>(N, "array_eq");
  array(0) = 2.5f;
  array(1) = 2.5f;
  array(2) = 2.5f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 2.5f);
}

TEST_F(CommonMacrosTest, FindMax1DSingleElement) {
  constexpr int N = 1;
  auto array = allocateVector<vectorReal>(N, "array_single");
  array(0) = 42.0f;
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 42.0f);
}

TEST_F(CommonMacrosTest, FindMax1DEmptyArrayThrows) {
  constexpr int N = 0;
  auto array = allocateVector<vectorReal>(N, "array_empty");
  EXPECT_THROW(helper.find_max(array, N), std::runtime_error);
  EXPECT_THROW(helper.find_max(array, 10), std::runtime_error);
}

// ======================================================================
// FIND_MIN
// ======================================================================

TEST_F(CommonMacrosTest, FindMin_Basic) {
  constexpr int N = 5;
  auto array = allocateVector<vectorReal>(N, "array_min");
  array(0) = 4.1f;
  array(1) = 1.0f;
  array(2) = 3.5f;
  array(3) = -2.0f;
  array(4) = 7.2f;
  EXPECT_FLOAT_EQ(helper.find_min(array, N), -2.0f);
}

TEST_F(CommonMacrosTest, FindMin_AllPositive) {
  constexpr int N = 3;
  auto array = allocateVector<vectorReal>(N, "array_min_pos");
  array(0) = 5.0f;
  array(1) = 2.5f;
  array(2) = 8.0f;
  EXPECT_FLOAT_EQ(helper.find_min(array, N), 2.5f);
}

TEST_F(CommonMacrosTest, FindMin_SingleElement) {
  constexpr int N = 1;
  auto array = allocateVector<vectorReal>(N, "array_min_single");
  array(0) = -99.0f;
  EXPECT_FLOAT_EQ(helper.find_min(array, N), -99.0f);
}

// ======================================================================
// SUM
// ======================================================================

TEST_F(CommonMacrosTest, Sum_KnownValues) {
  constexpr int N = 4;
  auto array = allocateVector<vectorReal>(N, "array_sum");
  array(0) = 1.0f;
  array(1) = 2.0f;
  array(2) = 3.0f;
  array(3) = 4.0f;
  EXPECT_FLOAT_EQ(helper.sum(array, N), 10.0f);
}

TEST_F(CommonMacrosTest, Sum_AllZero) {
  constexpr int N = 5;
  auto array = allocateVector<vectorReal>(N, "array_sum_zero");
  for (int i = 0; i < N; ++i) array(i) = 0.0f;
  EXPECT_FLOAT_EQ(helper.sum(array, N), 0.0f);
}

TEST_F(CommonMacrosTest, Sum_NegativeValues) {
  constexpr int N = 3;
  auto array = allocateVector<vectorReal>(N, "array_sum_neg");
  array(0) = -1.0f;
  array(1) = -2.0f;
  array(2) = 3.0f;
  EXPECT_FLOAT_EQ(helper.sum(array, N), 0.0f);
}

}  // namespace test
}  // namespace utils
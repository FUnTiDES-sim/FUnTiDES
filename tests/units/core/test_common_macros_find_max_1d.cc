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
    float result = 0.0f;
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

namespace {

/// Device array holding the given values.
vectorReal makeArray(const char* label, std::initializer_list<float> values) {
  auto array = allocateVector<vectorReal>(values.size(), label);
  auto host = Kokkos::create_mirror_view(array);
  int i = 0;
  for (float v : values) host(i++) = v;
  Kokkos::deep_copy(array, host);
  return array;
}

}  // namespace

TEST_F(CommonMacrosTest, FindMax1D) {
  constexpr int N = 5;
  auto array = makeArray("array", {4.1f, 1.0f, 3.5f, -2.0f, 7.2f});
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 7.2f);
}

TEST_F(CommonMacrosTest, FindMax1DLargeArray) {
  constexpr int N = 2000;
  auto array = allocateVector<vectorReal>(N, "array_large");
  Kokkos::deep_copy(array, 3000.0f);
  Kokkos::deep_copy(Kokkos::subview(array, N / 2), 3300.0f);  // One element is different
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 3300.0f);
}

TEST_F(CommonMacrosTest, FindMax1DAllNegative) {
  constexpr int N = 4;
  auto array = makeArray("array_neg", {-5.0f, -3.2f, -7.8f, -1.1f});
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, -1.1f);
}

TEST_F(CommonMacrosTest, FindMax1DAllEqual) {
  constexpr int N = 3;
  auto array = makeArray("array_eq", {2.5f, 2.5f, 2.5f});
  auto result = helper.find_max(array, N);
  EXPECT_FLOAT_EQ(result, 2.5f);
}

TEST_F(CommonMacrosTest, FindMax1DSingleElement) {
  constexpr int N = 1;
  auto array = makeArray("array_single", {42.0f});
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
  auto array = makeArray("array_min", {4.1f, 1.0f, 3.5f, -2.0f, 7.2f});
  EXPECT_FLOAT_EQ(helper.find_min(array, N), -2.0f);
}

TEST_F(CommonMacrosTest, FindMin_AllPositive) {
  constexpr int N = 3;
  auto array = makeArray("array_min_pos", {5.0f, 2.5f, 8.0f});
  EXPECT_FLOAT_EQ(helper.find_min(array, N), 2.5f);
}

TEST_F(CommonMacrosTest, FindMin_SingleElement) {
  constexpr int N = 1;
  auto array = makeArray("array_min_single", {-99.0f});
  EXPECT_FLOAT_EQ(helper.find_min(array, N), -99.0f);
}

// ======================================================================
// SUM
// ======================================================================

TEST_F(CommonMacrosTest, Sum_KnownValues) {
  constexpr int N = 4;
  auto array = makeArray("array_sum", {1.0f, 2.0f, 3.0f, 4.0f});
  EXPECT_FLOAT_EQ(helper.sum(array, N), 10.0f);
}

TEST_F(CommonMacrosTest, Sum_AllZero) {
  constexpr int N = 5;
  auto array = allocateVector<vectorReal>(N, "array_sum_zero");
  Kokkos::deep_copy(array, 0.0f);
  EXPECT_FLOAT_EQ(helper.sum(array, N), 0.0f);
}

TEST_F(CommonMacrosTest, Sum_NegativeValues) {
  constexpr int N = 3;
  auto array = makeArray("array_sum_neg", {-1.0f, -2.0f, 3.0f});
  EXPECT_FLOAT_EQ(helper.sum(array, N), 0.0f);
}

}  // namespace test
}  // namespace utils
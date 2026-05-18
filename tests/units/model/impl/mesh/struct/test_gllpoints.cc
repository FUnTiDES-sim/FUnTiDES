#include <gtest/gtest.h>

#include <cmath>

#include "model_struct.h"

namespace {

// ============================================================================
// num_points
// ============================================================================

TEST(GLLPointsTest, NumPointsEqualsOrderPlusOne) {
  for (int order = 1; order <= MAX_GLL_ORDER; ++order) EXPECT_EQ(GLLPoints::num_points(order), order + 1);
}

// ============================================================================
// Out-of-range inputs
// ============================================================================

TEST(GLLPointsTest, InvalidOrderZeroReturnsZero) { EXPECT_FLOAT_EQ(GLLPoints::get(0, 0), 0.0f); }

TEST(GLLPointsTest, InvalidOrderTooHighReturnsZero) { EXPECT_FLOAT_EQ(GLLPoints::get(MAX_GLL_ORDER + 1, 0), 0.0f); }

TEST(GLLPointsTest, InvalidIndexNegativeReturnsZero) { EXPECT_FLOAT_EQ(GLLPoints::get(3, -1), 0.0f); }

TEST(GLLPointsTest, InvalidIndexTooHighReturnsZero) {
  for (int order = 1; order <= MAX_GLL_ORDER; ++order) EXPECT_FLOAT_EQ(GLLPoints::get(order, order + 1), 0.0f);
}

// ============================================================================
// Boundary points always ±1
// ============================================================================

TEST(GLLPointsTest, FirstPointAlwaysMinusOne) {
  for (int order = 1; order <= MAX_GLL_ORDER; ++order)
    EXPECT_FLOAT_EQ(GLLPoints::get(order, 0), -1.0f) << "order=" << order;
}

TEST(GLLPointsTest, LastPointAlwaysPlusOne) {
  for (int order = 1; order <= MAX_GLL_ORDER; ++order)
    EXPECT_FLOAT_EQ(GLLPoints::get(order, order), 1.0f) << "order=" << order;
}

// ============================================================================
// Symmetry: get(order, i) == -get(order, order - i)
// ============================================================================

TEST(GLLPointsTest, SymmetryAboutZero) {
  for (int order = 1; order <= MAX_GLL_ORDER; ++order) {
    for (int i = 0; i <= order; ++i) {
      float xi = GLLPoints::get(order, i);
      float xj = GLLPoints::get(order, order - i);
      EXPECT_NEAR(xi, -xj, 1e-6f) << "order=" << order << " i=" << i;
    }
  }
}

// ============================================================================
// Monotonically increasing
// ============================================================================

TEST(GLLPointsTest, MonotonicallyIncreasing) {
  for (int order = 1; order <= MAX_GLL_ORDER; ++order) {
    for (int i = 0; i < order; ++i) {
      float xi = GLLPoints::get(order, i);
      float xj = GLLPoints::get(order, i + 1);
      EXPECT_LT(xi, xj) << "order=" << order << " i=" << i;
    }
  }
}

// ============================================================================
// All points in [-1, 1]
// ============================================================================

TEST(GLLPointsTest, AllPointsInUnitInterval) {
  for (int order = 1; order <= MAX_GLL_ORDER; ++order) {
    for (int i = 0; i <= order; ++i) {
      float xi = GLLPoints::get(order, i);
      EXPECT_GE(xi, -1.0f) << "order=" << order << " i=" << i;
      EXPECT_LE(xi, 1.0f) << "order=" << order << " i=" << i;
    }
  }
}

// ============================================================================
// Known values
// ============================================================================

TEST(GLLPointsTest, KnownValuesOrder1) {
  EXPECT_FLOAT_EQ(GLLPoints::get(1, 0), -1.0f);
  EXPECT_FLOAT_EQ(GLLPoints::get(1, 1), 1.0f);
}

TEST(GLLPointsTest, KnownValuesOrder2) {
  EXPECT_FLOAT_EQ(GLLPoints::get(2, 0), -1.0f);
  EXPECT_FLOAT_EQ(GLLPoints::get(2, 1), 0.0f);
  EXPECT_FLOAT_EQ(GLLPoints::get(2, 2), 1.0f);
}

TEST(GLLPointsTest, KnownValuesOrder3) {
  EXPECT_FLOAT_EQ(GLLPoints::get(3, 0), -1.0f);
  EXPECT_NEAR(GLLPoints::get(3, 1), -0.4472135955f, 1e-6f);
  EXPECT_NEAR(GLLPoints::get(3, 2), 0.4472135955f, 1e-6f);
  EXPECT_FLOAT_EQ(GLLPoints::get(3, 3), 1.0f);
}

TEST(GLLPointsTest, KnownValuesOrder4InteriorZero) {
  // Order 4: interior point at i=2 is 0
  EXPECT_FLOAT_EQ(GLLPoints::get(4, 2), 0.0f);
}

TEST(GLLPointsTest, KnownValuesOrder6InteriorZero) {
  // Order 6: interior point at i=3 is 0
  EXPECT_FLOAT_EQ(GLLPoints::get(6, 3), 0.0f);
}

TEST(GLLPointsTest, KnownValuesOrder8InteriorZero) {
  // Order 8: interior point at i=4 is 0
  EXPECT_FLOAT_EQ(GLLPoints::get(8, 4), 0.0f);
}

}  // namespace

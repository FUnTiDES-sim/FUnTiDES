#include <gtest/gtest.h>

#include <cmath>

#include "source_and_receiver_utils.h"

namespace utils {
namespace test {

using SourceAndReceiverUtils::ComputeRHSWeights;

class RHSWeightsTest : public ::testing::Test {
 protected:
  void makeCornerCoords(float x0, float y0, float z0, float dx, float dy, float dz, real_t (&corners)[8][3]) {
    float cx[8] = {x0, x0 + dx, x0, x0 + dx, x0, x0 + dx, x0, x0 + dx};
    float cy[8] = {y0, y0, y0 + dy, y0 + dy, y0, y0, y0 + dy, y0 + dy};
    float cz[8] = {z0, z0, z0, z0, z0 + dz, z0 + dz, z0 + dz, z0 + dz};
    for (int i = 0; i < 8; ++i) {
      corners[i][0] = cx[i];
      corners[i][1] = cy[i];
      corners[i][2] = cz[i];
    }
  }
};

// ============================================================================
// Order 2 — center point — sum = 1, center node heaviest
// ============================================================================

TEST_F(RHSWeightsTest, ComputeRHSWeightsOrder2) {
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);
  std::array<float, 3> coord = {5.0f, 5.0f, 5.0f};
  arrayReal::host_mirror_type rhsWeights("rhs", 1, 27);
  ComputeRHSWeights<2>(corners, coord, rhsWeights);
  float sum = 0;
  for (int i = 0; i < 27; ++i) sum += rhsWeights(0, i);
  EXPECT_NEAR(sum, 1.0f, 1e-5f);
  EXPECT_GT(rhsWeights(0, 13), rhsWeights(0, 0));
}

// ============================================================================
// Order 1 — center point — all 8 weights = 0.125
// ============================================================================

TEST_F(RHSWeightsTest, ComputeRHSWeightsOrder1) {
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);
  std::array<float, 3> coord = {5.0f, 5.0f, 5.0f};
  arrayReal::host_mirror_type rhsWeights("rhs", 1, 8);
  ComputeRHSWeights<1>(corners, coord, rhsWeights);
  for (int i = 0; i < 8; ++i) EXPECT_NEAR(rhsWeights(0, i), 0.125f, 1e-5f);
}

// ============================================================================
// Off-center point — sum = 1 (partition of unity)
// ============================================================================

TEST_F(RHSWeightsTest, ComputeRHSWeightsSumToOne) {
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);
  std::array<float, 3> coord = {2.5f, 7.5f, 5.0f};
  arrayReal::host_mirror_type rhsWeights("rhs", 1, 27);
  ComputeRHSWeights<2>(corners, coord, rhsWeights);
  float sum = 0;
  for (int i = 0; i < 27; ++i) sum += rhsWeights(0, i);
  EXPECT_NEAR(sum, 1.0f, 1e-5f);
}

// ============================================================================
// Outside element — weights finite (Lagrange can be negative outside)
// ============================================================================

TEST_F(RHSWeightsTest, ComputeRHSWeightsOutsideElementFinite) {
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);
  std::array<float, 3> coord = {15.0f, 5.0f, 5.0f};
  arrayReal::host_mirror_type rhsWeights("rhs", 1, 8);
  ComputeRHSWeights<1>(corners, coord, rhsWeights);
  for (int i = 0; i < 8; ++i) EXPECT_TRUE(std::isfinite(rhsWeights(0, i)));
}

// ============================================================================
// Corner point — weight = 1 at node 0, 0 elsewhere
// ============================================================================

TEST_F(RHSWeightsTest, ComputeRHSWeightsAtCorner) {
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);
  std::array<float, 3> coord = {0.0f, 0.0f, 0.0f};
  arrayReal::host_mirror_type rhsWeights("rhs", 1, 8);
  ComputeRHSWeights<1>(corners, coord, rhsWeights);
  EXPECT_NEAR(rhsWeights(0, 0), 1.0f, 1e-5f);
  for (int i = 1; i < 8; ++i) EXPECT_NEAR(rhsWeights(0, i), 0.0f, 1e-5f);
}

// ============================================================================
// Face center (z=0) — nodes 0,1,2,3 each = 0.25, others = 0
// Corner ordering: node k = (k&1 ? x1:x0, k&2 ? y1:y0, k&4 ? z1:z0)
// z=0 face: nodes 0(x0,y0,0), 1(x1,y0,0), 2(x0,y1,0), 3(x1,y1,0)
// ============================================================================

TEST_F(RHSWeightsTest, ComputeRHSWeightsAtFaceCenter) {
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);
  std::array<float, 3> coord = {5.0f, 5.0f, 0.0f};
  arrayReal::host_mirror_type rhsWeights("rhs", 1, 8);
  ComputeRHSWeights<1>(corners, coord, rhsWeights);
  EXPECT_NEAR(rhsWeights(0, 0), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 1), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 2), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 3), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 4), 0.0f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 5), 0.0f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 6), 0.0f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 7), 0.0f, 1e-5f);
}

// ============================================================================
// Face center (y=0) — nodes 0,1,4,5 each = 0.25
// y=0 face: nodes 0(x0,0,z0), 1(x1,0,z0), 4(x0,0,z1), 5(x1,0,z1)
// ============================================================================

TEST_F(RHSWeightsTest, ComputeRHSWeightsAtYFaceCenter) {
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);
  std::array<float, 3> coord = {5.0f, 0.0f, 5.0f};
  arrayReal::host_mirror_type rhsWeights("rhs", 1, 8);
  ComputeRHSWeights<1>(corners, coord, rhsWeights);
  EXPECT_NEAR(rhsWeights(0, 0), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 1), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 4), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 5), 0.25f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 2), 0.0f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 3), 0.0f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 6), 0.0f, 1e-5f);
  EXPECT_NEAR(rhsWeights(0, 7), 0.0f, 1e-5f);
}

}  // namespace test
}  // namespace utils

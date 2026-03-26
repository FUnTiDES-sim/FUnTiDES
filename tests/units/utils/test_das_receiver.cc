#include <gtest/gtest.h>

#include <cmath>

#include "source_and_receiver_utils.h"

namespace utils
{
namespace test
{

using SourceAndReceiverUtils::ComputeDASVector;
using SourceAndReceiverUtils::ComputeDASWeightsForSample;
using SourceAndReceiverUtils::DASType;

class DASReceiverTest : public ::testing::Test
{
 protected:
  // Helper: build corner coordinates for a Cartesian element [x0,x1]^3
  void makeCornerCoords(float x0, float y0, float z0, float dx, float dy,
                        float dz, real_t (&corners)[8][3])
  {
    // Corner ordering: (0,0,0),(1,0,0),(0,1,0),(1,1,0),
    //                  (0,0,1),(1,0,1),(0,1,1),(1,1,1)
    float cx[8] = {x0, x0 + dx, x0, x0 + dx, x0, x0 + dx, x0, x0 + dx};
    float cy[8] = {y0, y0, y0 + dy, y0 + dy, y0, y0, y0 + dy, y0 + dy};
    float cz[8] = {z0, z0, z0, z0, z0 + dz, z0 + dz, z0 + dz, z0 + dz};
    for (int i = 0; i < 8; ++i)
    {
      corners[i][0] = cx[i];
      corners[i][1] = cy[i];
      corners[i][2] = cz[i];
    }
  }
};

// ---------------------------------------------------------------
// Test ComputeDASVector
// ---------------------------------------------------------------
TEST_F(DASReceiverTest, DASVectorHorizontalX)
{
  // dip=0, azimuth=0 → direction should be (1,0,0)
  auto v = ComputeDASVector(0.0f, 0.0f);
  EXPECT_NEAR(v[0], 1.0f, 1e-6f);
  EXPECT_NEAR(v[1], 0.0f, 1e-6f);
  EXPECT_NEAR(v[2], 0.0f, 1e-6f);
}

TEST_F(DASReceiverTest, DASVectorHorizontalY)
{
  // dip=0, azimuth=90 → direction should be (0,1,0)
  auto v = ComputeDASVector(0.0f, 90.0f);
  EXPECT_NEAR(v[0], 0.0f, 1e-5f);
  EXPECT_NEAR(v[1], 1.0f, 1e-5f);
  EXPECT_NEAR(v[2], 0.0f, 1e-5f);
}

TEST_F(DASReceiverTest, DASVectorVertical)
{
  // dip=90, azimuth=0 → direction should be (0,0,1)
  auto v = ComputeDASVector(90.0f, 0.0f);
  EXPECT_NEAR(v[0], 0.0f, 1e-5f);
  EXPECT_NEAR(v[1], 0.0f, 1e-5f);
  EXPECT_NEAR(v[2], 1.0f, 1e-5f);
}

TEST_F(DASReceiverTest, DASVectorIsUnitLength)
{
  auto v = ComputeDASVector(30.0f, 45.0f);
  float len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  EXPECT_NEAR(len, 1.0f, 1e-6f);
}

// ---------------------------------------------------------------
// Test ComputeDASWeightsForSample — Dipole mode
// ---------------------------------------------------------------
TEST_F(DASReceiverTest, DipoleWeightsSumToOne)
{
  // For dipole mode with integrationConstant=1, the weights should be
  // shape function values at the sample point. They should sum to 1.
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 100.0f, 100.0f, 100.0f, corners);

  // Sample at element center (50, 50, 50)
  std::array<float, 3> coord = {50.0f, 50.0f, 50.0f};
  std::array<float, 3> direction = {1.0f, 0.0f, 0.0f};

  constexpr int ORDER = 2;
  constexpr int npe =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;
  float weights[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kDipole, weights);

  float sum = 0;
  for (int i = 0; i < npe; ++i) sum += weights[i];
  EXPECT_NEAR(sum, 1.0f, 1e-5f);
}

TEST_F(DASReceiverTest, DipoleWeightsAtCorner)
{
  // Sample at corner (0,0,0) → only the corner basis function should be 1
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 100.0f, 100.0f, 100.0f, corners);

  std::array<float, 3> coord = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> direction = {1.0f, 0.0f, 0.0f};

  constexpr int ORDER = 1;
  constexpr int npe =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;
  float weights[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kDipole, weights);

  // Node 0 (a=0,b=0,c=0) should have weight 1, others 0
  EXPECT_NEAR(weights[0], 1.0f, 1e-5f);
  for (int i = 1; i < npe; ++i) EXPECT_NEAR(weights[i], 0.0f, 1e-5f);
}

TEST_F(DASReceiverTest, DipoleWeightsAllPositive)
{
  // For a point inside the element, all basis function values should
  // be non-negative (Lagrange on GLL points are not always positive,
  // but at the center they should be).
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 50.0f, 50.0f, 50.0f, corners);

  std::array<float, 3> coord = {25.0f, 25.0f, 25.0f};
  std::array<float, 3> direction = {1.0f, 0.0f, 0.0f};

  constexpr int ORDER = 1;
  constexpr int npe =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;
  float weights[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kDipole, weights);

  // For linear elements, all shape functions are positive inside the element
  for (int i = 0; i < npe; ++i) EXPECT_GE(weights[i], 0.0f);
}

// ---------------------------------------------------------------
// Test ComputeDASWeightsForSample — Strain Integration mode
// ---------------------------------------------------------------
TEST_F(DASReceiverTest, StrainWeightsGradientSumZero)
{
  // Gradient of shape functions should sum to zero (partition of unity).
  // So for strain integration at any point, sum of weights should be ~0.
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 100.0f, 100.0f, 100.0f, corners);

  std::array<float, 3> coord = {50.0f, 50.0f, 50.0f};
  std::array<float, 3> direction = {1.0f, 0.0f, 0.0f};

  constexpr int ORDER = 2;
  constexpr int npe =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;
  float weights[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kStrainIntegration, weights);

  float sum = 0;
  for (int i = 0; i < npe; ++i) sum += weights[i];
  EXPECT_NEAR(sum, 0.0f, 1e-4f);
}

TEST_F(DASReceiverTest, StrainWeightsLinearFieldGivesConstantStrain)
{
  // For a linear displacement field u_x = x, the strain du_x/dx = 1.
  // DAS weights dotted with node displacements should give ≈1.
  real_t corners[8][3];
  float x0 = 100.0f, y0 = 200.0f, z0 = 300.0f;
  float dx = 50.0f, dy = 50.0f, dz = 50.0f;
  makeCornerCoords(x0, y0, z0, dx, dy, dz, corners);

  // Sample at center of element
  std::array<float, 3> coord = {x0 + dx / 2, y0 + dy / 2, z0 + dz / 2};
  std::array<float, 3> direction = {1.0f, 0.0f, 0.0f};

  constexpr int ORDER = 2;
  using FEType =
      typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
  constexpr int npe = FEType::numNodes;
  float weights[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kStrainIntegration, weights);

  // Assign node displacements: u_x = x_node (linear in x)
  // Need to compute node positions. For a Cartesian element with ORDER=2,
  // nodes are at GLL points. For order 2 GLL points are {-1, 0, 1}.
  // In physical space, node (a,b,c) is at:
  //   x = x0 + dx*(a/ORDER), y = y0 + dy*(b/ORDER), z = z0 + dz*(c/ORDER)
  // (This is approximate for GLL but exact for equally spaced GLL at order 2)
  float strain = 0;
  for (int c = 0; c <= ORDER; ++c)
    for (int b = 0; b <= ORDER; ++b)
      for (int a = 0; a <= ORDER; ++a)
      {
        int nodeIdx = a + b * (ORDER + 1) + c * (ORDER + 1) * (ORDER + 1);
        // GLL points for order 2: xi = {-1, 0, 1} → physical x = x0 +
        // dx*(xi+1)/2 More precisely: node a of (ORDER+1) → xi_a. For ORDER=2:
        // xi_0=-1, xi_1=0, xi_2=1.
        float nodeX = x0 + dx * static_cast<float>(a) / ORDER;
        strain += weights[nodeIdx] * nodeX;
      }

  // du_x/dx = 1 everywhere, dotted with direction (1,0,0) gives strain_xx = 1
  EXPECT_NEAR(strain, 1.0f, 1e-3f);
}

TEST_F(DASReceiverTest, StrainWeightsOrder1)
{
  // Test strain weights for order 1 element
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 10.0f, 10.0f, 10.0f, corners);

  std::array<float, 3> coord = {5.0f, 5.0f, 5.0f};
  std::array<float, 3> direction = {0.0f, 1.0f, 0.0f};  // Y-direction

  constexpr int ORDER = 1;
  constexpr int npe =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;
  float weights[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kStrainIntegration, weights);

  // Gradients of shape functions in Y direction should sum to 0
  float sum = 0;
  for (int i = 0; i < npe; ++i) sum += weights[i];
  EXPECT_NEAR(sum, 0.0f, 1e-5f);

  // For linear elements: dN/dy for node (a,b,c):
  // dN/dy = N_a(x)*dN_b(y)/dy*N_c(z)
  // At center: N_a(0.5)=0.5, dN_0/dy=-1/dy, dN_1/dy=1/dy, N_c(0.5)=0.5
  // So weights should be ±0.25/dy = ±0.025
  float invDy = 1.0f / 10.0f;
  for (int c = 0; c < 2; ++c)
    for (int a = 0; a < 2; ++a)
    {
      int n0 = a + 0 * 2 + c * 4;  // b=0
      int n1 = a + 1 * 2 + c * 4;  // b=1
      EXPECT_NEAR(weights[n0], -0.25f * invDy, 1e-5f);
      EXPECT_NEAR(weights[n1], 0.25f * invDy, 1e-5f);
    }
}

TEST_F(DASReceiverTest, StrainWeightsOrder3)
{
  // Verify that order 3 strain weights have the right size and sum to 0
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 100.0f, 100.0f, 100.0f, corners);

  std::array<float, 3> coord = {50.0f, 50.0f, 50.0f};
  std::array<float, 3> direction = {0.0f, 0.0f, 1.0f};

  constexpr int ORDER = 3;
  constexpr int npe =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;
  static_assert(npe == 64, "Order 3 should have 64 nodes");
  float weights[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kStrainIntegration, weights);

  float sum = 0;
  for (int i = 0; i < npe; ++i) sum += weights[i];
  EXPECT_NEAR(sum, 0.0f, 1e-3f);
}

// ---------------------------------------------------------------
// Test integration constant scaling
// ---------------------------------------------------------------
TEST_F(DASReceiverTest, IntegrationConstantScaling)
{
  real_t corners[8][3];
  makeCornerCoords(0.0f, 0.0f, 0.0f, 100.0f, 100.0f, 100.0f, corners);

  std::array<float, 3> coord = {50.0f, 50.0f, 50.0f};
  std::array<float, 3> direction = {1.0f, 0.0f, 0.0f};

  constexpr int ORDER = 2;
  constexpr int npe =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;

  float weights1[npe] = {0};
  float weights2[npe] = {0};

  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 1.0f,
                                    DASType::kDipole, weights1);
  ComputeDASWeightsForSample<ORDER>(corners, coord, direction, 3.0f,
                                    DASType::kDipole, weights2);

  for (int i = 0; i < npe; ++i)
    EXPECT_NEAR(weights2[i], 3.0f * weights1[i], 1e-6f);
}

}  // namespace test
}  // namespace utils

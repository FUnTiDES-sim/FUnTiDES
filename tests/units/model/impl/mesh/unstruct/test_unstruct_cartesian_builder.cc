#include <gtest/gtest.h>

#include <memory>

// Assuming these are the headers for your classes
#include "cartesian_unstruct_builder.h"

namespace model
{
namespace test
{

class CartesianBuilderTest : public ::testing::Test
{
 protected:
  using FloatT = double;
  using ScalarT = int;

  // Helper to create basic parameters
  CartesianParams<FloatT, ScalarT> createParams(bool onNodes, int order)
  {
    CartesianParams<FloatT, ScalarT> p;
    p.ex = 2;
    p.ey = 1;
    p.ez = 1;
    p.lx = 10.0;
    p.ly = 5.0;
    p.lz = 5.0;
    p.order = order;
    p.isModelOnNodes = onNodes;
    p.isElastic = true;
    return p;
  }
};

// ============================================================================
// TEST 1: Optimized Corner Mode (Properties on Elements)
// ============================================================================
TEST_F(CartesianBuilderTest, OptimizedCornerMeshVerification)
{
  // We request Order 4, but since properties are on elements,
  // the builder should optimize to Order 1 geometry (8 corners).
  auto params = createParams(false, 4);
  CartesianUnstructBuilder<FloatT, ScalarT> builder(params);
  auto model = builder.getModel();

  // 1. Check Node Counts
  // For 2x1x1 elements, vertices are (2+1)*(1+1)*(1+1) = 3 * 2 * 2 = 12
  EXPECT_EQ(model->getNumberOfNodes(), 12);
  EXPECT_EQ(model->getNumberOfElements(), 2);

  // The internal order for the UNSTRUCT model should now be 1 (linear)
  EXPECT_EQ(model->getOrder(), 1);
  EXPECT_EQ(model->getNumberOfPointsPerElement(), 8);

  // 2. Check Coordinates (Should be uniform)
  // Element size dx = 10/2 = 5.0
  // Node at the middle interface (x=1, y=0, z=0) should be at X=5.0
  // In a 3x2x2 grid, global index of (1,0,0) is 1.
  EXPECT_NEAR(model->nodeCoord(1, 0), 5.0, 1e-7);

  // Node at the far end (x=2, y=1, z=1) is index (2 + 1*3 + 1*6) = 11
  EXPECT_NEAR(model->nodeCoord(11, 0), 10.0, 1e-7);
  EXPECT_NEAR(model->nodeCoord(11, 1), 5.0, 1e-7);

  // 3. Verify Properties are on Elements
  EXPECT_FALSE(model->isModelOnNodes());
  // Vp should be 1500 (default) on elements
  EXPECT_NEAR(model->getModelVpOnElement(0), 1500.0, 1e-7);
}

// ============================================================================
// TEST 2: High-Order GLL Mode (Properties on Nodes)
// ============================================================================
TEST_F(CartesianBuilderTest, HighOrderGLLMeshVerification)
{
  // Properties on nodes, Order 2 (Quadratic)
  auto params = createParams(true, 2);
  CartesianUnstructBuilder<FloatT, ScalarT> builder(params);
  auto model = builder.getModel();

  // 1. Check Node Counts
  // For 2x1x1 elements at Order 2:
  // nx = (2 * 2) + 1 = 5
  // ny = (1 * 2) + 1 = 3
  // nz = (1 * 2) + 1 = 3
  // Total nodes = 5 * 3 * 3 = 45
  EXPECT_EQ(model->getNumberOfNodes(), 45);
  EXPECT_EQ(model->getOrder(), 2);
  EXPECT_EQ(model->getNumberOfPointsPerElement(), 27);  // 3^3

  // 2. Check GLL Coordinate (Mid-node of element 0)
  // Element 0 spans X: [0, 5]. At order 2, the middle node is at 2.5.
  // Local indices for center: i=1, j=1, k=1.
  ScalarT midNode = model->globalNodeIndex(0, 1, 1, 1);
  EXPECT_NEAR(model->nodeCoord(midNode, 0), 2.5, 1e-7);

  // 3. Verify Properties are on Nodes
  EXPECT_TRUE(model->isModelOnNodes());
  EXPECT_NEAR(model->getModelVpOnNodes(0), 1500.0, 1e-7);
}

// ============================================================================
// TEST 3: Connectivity Sharing
// ============================================================================
TEST_F(CartesianBuilderTest, ElementConnectivitySharing)
{
  auto params = createParams(false, 1);  // Simple linear
  CartesianUnstructBuilder<FloatT, ScalarT> builder(params);
  auto model = builder.getModel();

  // Element 0 (left) and Element 1 (right) share a face.
  // Element 0's local right corners (i=1) must equal Element 1's local left
  // corners (i=0)

  // Check bottom-front shared corner
  ScalarT shared0 = model->globalNodeIndex(0, 1, 0, 0);
  ScalarT shared1 = model->globalNodeIndex(1, 0, 0, 0);
  EXPECT_EQ(shared0, shared1);
  EXPECT_EQ(shared0, 1);  // Lexicographical (1, 0, 0) in 3x2x2 grid

  // Check top-back shared corner
  ScalarT sharedUpper0 = model->globalNodeIndex(0, 1, 1, 1);
  ScalarT sharedUpper1 = model->globalNodeIndex(1, 0, 1, 1);
  EXPECT_EQ(sharedUpper0, sharedUpper1);
}

}  // namespace test
}  // namespace model

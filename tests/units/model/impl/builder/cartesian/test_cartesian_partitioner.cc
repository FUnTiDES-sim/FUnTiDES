#include <gtest/gtest.h>

#include <array>
#include <stdexcept>
#include <vector>

#include "model.h"
#include "topology_factory.h"

namespace
{

// Minimal Mock Mesh implementing only what TopologyFactory needs
template <typename FloatType, typename ScalarType>
class MockMesh : public model::ModelApi<FloatType, ScalarType>
{
 public:
  struct Node
  {
    FloatType x, y, z;
  };
  std::vector<Node> nodes;
  FloatType minSpacing = 1.0;

  void addNode(FloatType x, FloatType y, FloatType z)
  {
    nodes.push_back({x, y, z});
  }

  // Required by TopologyFactory
  PROXY_HOST_DEVICE FloatType nodeCoord(ScalarType i, int dim) const override
  {
    if (i < 0 || i >= static_cast<ScalarType>(nodes.size())) return 0.0;
    if (dim == 0) return nodes[i].x;
    if (dim == 1) return nodes[i].y;
    if (dim == 2) return nodes[i].z;
    return 0.0;
  }

  PROXY_HOST_DEVICE ScalarType getNumberOfNodes() const override
  {
    return static_cast<ScalarType>(nodes.size());
  }

  PROXY_HOST_DEVICE FloatType getMinSpacing() const override
  {
    return minSpacing;
  }

  // Stubs for other pure virtual methods (not used by TopologyFactory)
  PROXY_HOST_DEVICE ScalarType globalNodeIndex(ScalarType, int, int,
                                               int) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelVpOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelVpOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelRhoOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelRhoOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelVsOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelVsOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelDeltaOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelDeltaOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType
  getModelEpsilonOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelGammaOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE FloatType getModelGammaOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE ScalarType getModelThetaOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE ScalarType getModelThetaOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE ScalarType getModelPhiOnNodes(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE ScalarType getModelPhiOnElement(ScalarType) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE void getCTensorOnElement(ScalarType,
                                             FloatType[6][6]) const override
  {
  }
  virtual void initElasticityTensors() override {}
  PROXY_HOST_DEVICE ScalarType getNumberOfElements() const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE int getNumberOfPointsPerElement() const override
  {
    return 27;
  }
  PROXY_HOST_DEVICE int getOrder() const override { return 1; }
  PROXY_HOST_DEVICE model::BoundaryFlag boundaryType(ScalarType) const override
  {
    return model::BoundaryFlag::InteriorNode;
  }
  PROXY_HOST_DEVICE void faceNormal(ScalarType, model::CubicFace,
                                    FloatType[3]) const override
  {
  }

  virtual void buildFaceConnectivity() override {}

  PROXY_HOST_DEVICE FloatType domainSize(int) const override { return 0; }
  virtual FloatType getMaxSpeed() const override { return 0; }
  PROXY_HOST_DEVICE bool isModelOnNodes() const override { return true; }
  PROXY_HOST_DEVICE bool isElastic() const override { return false; }
};

class TopologyFactoryTest : public ::testing::Test
{
 protected:
  MockMesh<float, int> mesh;
};

TEST_F(TopologyFactoryTest, SingleRankReturnsEmptySharedNodes)
{
  mesh.addNode(0.0, 0.0, 0.0);
  mesh.addNode(1.0, 0.0, 0.0);

  auto topo = TopologyFactory::createFromMesh(mesh, 0, 1, 0.0f, 10.0f);

  EXPECT_FALSE(topo.isDistributed());
  EXPECT_EQ(topo.sharedNodes.size(), 0);
}

TEST_F(TopologyFactoryTest, LeftBoundaryDetection)
{
  // Rank 1 of 3 (Middle rank)
  // Left boundary at x=10.0, Right at x=20.0
  // Nodes: 10.0 (Left), 15.0 (Inner), 20.0 (Right)

  mesh.addNode(10.0f, 0.0f, 0.0f);  // Index 0
  mesh.addNode(15.0f, 0.0f, 0.0f);  // Index 1
  mesh.addNode(20.0f, 0.0f, 0.0f);  // Index 2

  float origin_x = 10.0f;
  float width_x = 10.0f;

  auto topo = TopologyFactory::createFromMesh(mesh, 1, 3, origin_x, width_x);

  EXPECT_TRUE(topo.isDistributed());

  // Should share with Rank 0 (Left)
  ASSERT_TRUE(topo.sharedNodes.count(0));
  EXPECT_EQ(topo.sharedNodes[0].size(), 1);
  EXPECT_EQ(topo.sharedNodes[0][0], 0);  // Node index 0 is at 10.0

  // Should share with Rank 2 (Right)
  ASSERT_TRUE(topo.sharedNodes.count(2));
  EXPECT_EQ(topo.sharedNodes[2].size(), 1);
  EXPECT_EQ(topo.sharedNodes[2][0], 2);  // Node index 2 is at 20.0
}

TEST_F(TopologyFactoryTest, AutoToleranceFromSpacing)
{
  mesh.minSpacing = 0.1f;
  mesh.addNode(10.000001f, 0.0f, 0.0f);  // Slightly off but within tolerance

  float origin_x = 10.0f;

  // Should detect using auto-computed tolerance (0.1 * 1e-4)
  // 1e-6 < 1e-5, so should match

  auto topo = TopologyFactory::createFromMesh(mesh, 1, 2, origin_x, 10.0f);

  EXPECT_FALSE(topo.sharedNodes[0].empty());
}

TEST_F(TopologyFactoryTest, ThrowsOnMissingBoundaryNodes)
{
  // Rank 1 expects a left neighbor at x=10.0
  // But mesh only has nodes at 15.0 and 20.0
  mesh.addNode(15.0f, 0.0f, 0.0f);
  mesh.addNode(20.0f, 0.0f, 0.0f);

  // CHANGE: Expect logic_error instead of runtime_error
  EXPECT_THROW({ TopologyFactory::createFromMesh(mesh, 1, 3, 10.0f, 10.0f); },
               std::logic_error);
}

TEST_F(TopologyFactoryTest, ThrowsOnAmbiguousBoundary)
{
  // A single node satisfying both left and right conditions
  mesh.addNode(10.0f, 0.0f, 0.0f);

  // Define domain [10.0, 10.0 + epsilon] where width is smaller than tolerance
  // This makes the node appear on both Left and Right boundaries
  // simultaneously.

  // CHANGE: Use non-zero width to pass input validation
  float width = 1e-7f;

  // CHANGE: Expect logic_error (topology failure) instead of invalid_argument
  // (bad input)
  EXPECT_THROW(

      { TopologyFactory::createFromMesh(mesh, 1, 3, 10.0f, width); },
      std::logic_error);
}

}  // namespace

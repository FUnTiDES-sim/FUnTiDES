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

  PROXY_HOST_DEVICE FloatType nodeCoord(ScalarType i, int dim) const override
  {
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
    if (i < 0 || i >= static_cast<ScalarType>(nodes.size())) return 0.0;
    if (dim == 0) return nodes[i].x;
    if (dim == 1) return nodes[i].y;
    if (dim == 2) return nodes[i].z;
#endif
    return 0.0;
  }

  PROXY_HOST_DEVICE ScalarType getNumberOfNodes() const override
  {
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
    return static_cast<ScalarType>(nodes.size());
#else
    return 0;
#endif
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
  PROXY_HOST_DEVICE FloatType getModelQpOnNodes(ScalarType) const override
  {
    return 1e9;
  }
  PROXY_HOST_DEVICE FloatType getModelQpOnElement(ScalarType) const override
  {
    return 1e9;
  }
  PROXY_HOST_DEVICE FloatType getModelQsOnNodes(ScalarType) const override
  {
    return 1e9;
  }
  PROXY_HOST_DEVICE FloatType getModelQsOnElement(ScalarType) const override
  {
    return 1e9;
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
  virtual void initElasticityTensors(
      model::AnisotropyType anisotropy_type) override
  {
  }
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

  PROXY_HOST_DEVICE bool isFreeSurface(ScalarType n) const override
  {
    return true;
  }

  PROXY_HOST_DEVICE ScalarType getGlobalNodeFromFace(ScalarType,
                                                     int) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE ScalarType getGlobalFace(ScalarType,
                                             model::CubicFace) const override
  {
    return 0;
  }
  PROXY_HOST_DEVICE bool isBoundaryFace(ScalarType) const override
  {
    return true;
  }
  PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const override { return 0; }

  PROXY_HOST_DEVICE FloatType domainSize(int) const override { return 0; }
  virtual FloatType getMaxSpeed() const override { return 0; }
  PROXY_HOST_DEVICE bool isModelOnNodes() const override { return true; }
  PROXY_HOST_DEVICE bool isElastic() const override { return false; }
  void setQualityFactors(FloatType qp, FloatType qs) override {}
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
  mesh.addNode(10.0f, 0.0f, 0.0f);  // Left
  mesh.addNode(15.0f, 0.0f, 0.0f);
  mesh.addNode(20.0f, 0.0f, 0.0f);  // Right

  auto topo = TopologyFactory::createFromMesh(mesh, 1, 3, 10.0f, 10.0f);

  EXPECT_TRUE(topo.isDistributed());
  ASSERT_TRUE(topo.sharedNodes.count(0));
  EXPECT_EQ(topo.sharedNodes[0].size(), 1);
  EXPECT_EQ(topo.sharedNodes[0][0], 0);

  ASSERT_TRUE(topo.sharedNodes.count(2));
  EXPECT_EQ(topo.sharedNodes[2].size(), 1);
  EXPECT_EQ(topo.sharedNodes[2][0], 2);
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
  mesh.addNode(15.0f, 0.0f, 0.0f);
  mesh.addNode(20.0f, 0.0f, 0.0f);

  // Expect logic_error (topology inconsistency)
  EXPECT_THROW(
      { TopologyFactory::createFromMesh(mesh, 1, 3, 10.0f, 10.0f); },
      std::logic_error);
}

TEST_F(TopologyFactoryTest, ThrowsOnAmbiguousBoundary)
{
  mesh.addNode(10.0f, 0.0f, 0.0f);

  // Width 1e-7 is < tolerance 1e-6, so node is detected on BOTH boundaries
  float width = 1e-7f;

  EXPECT_THROW(
      { TopologyFactory::createFromMesh(mesh, 1, 3, 10.0f, width); },
      std::logic_error);
}

}  // namespace

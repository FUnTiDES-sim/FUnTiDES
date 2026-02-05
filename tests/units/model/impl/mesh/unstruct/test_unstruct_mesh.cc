#include <gtest/gtest.h>

#include "model_unstruct.h"

namespace model
{
namespace
{

// ============================================================================
// Test Fixture
// ============================================================================

template <typename T>
class ModelUnstructTest : public ::testing::Test
{
 protected:
  using FloatType = typename T::FloatType;
  using ScalarType = typename T::ScalarType;
  using ModelType = ModelUnstruct<FloatType, ScalarType>;

  // Helper to create minimal test mesh
  ModelType createTestMesh()
  {
    ModelUnstructData<FloatType, ScalarType> data;
    data.order_ = 1;
    data.n_element_ = 1;
    data.n_node_ = 8;
    data.lx_ = data.ly_ = data.lz_ = 1.0;
    data.isModelOnNodes_ = false;
    data.isElastic_ = false;

    data.global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(1, 8);
    data.nodes_coords_x_ = allocateVector<VECTOR_REAL_VIEW>(8);
    data.nodes_coords_y_ = allocateVector<VECTOR_REAL_VIEW>(8);
    data.nodes_coords_z_ = allocateVector<VECTOR_REAL_VIEW>(8);
    data.model_vp_element_ = allocateVector<VECTOR_REAL_VIEW>(1);
    data.model_rho_element_ = allocateVector<VECTOR_REAL_VIEW>(1);
    data.boundaries_t_ = allocateVector<VECTOR_REAL_VIEW>(8);

    for (int i = 0; i < 8; ++i)
    {
      data.global_node_index_(0, i) = i;
      data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
      data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
      data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
    }
    data.model_vp_element_[0] = 1500.0;
    data.model_rho_element_[0] = 1.0;

    return ModelType(data);
  }
};

struct FloatInt
{
  using FloatType = float;
  using ScalarType = int;
};
struct DoubleInt
{
  using FloatType = double;
  using ScalarType = int;
};

using TestTypes = ::testing::Types<FloatInt, DoubleInt>;
TYPED_TEST_SUITE(ModelUnstructTest, TestTypes);

// ============================================================================
// Construction Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, DefaultConstructor)
{
  typename TestFixture::ModelType model;
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, ConstructFromData)
{
  auto model = this->createTestMesh();
  EXPECT_EQ(model.getNumberOfElements(), 1);
  EXPECT_EQ(model.getNumberOfNodes(), 8);
  EXPECT_EQ(model.getOrder(), 1);
}

// ============================================================================
// Indexing Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, GlobalNodeIndexInRange)
{
  auto model = this->createTestMesh();

  for (int k = 0; k <= 1; ++k)
    for (int j = 0; j <= 1; ++j)
      for (int i = 0; i <= 1; ++i)
      {
        auto node = model.globalNodeIndex(0, i, j, k);
        EXPECT_GE(node, 0);
        EXPECT_LT(node, 8);
      }
}

TYPED_TEST(ModelUnstructTest, VertexCoordsInUnitCube)
{
  auto model = this->createTestMesh();

  for (int i = 0; i < 8; ++i)
  {
    typename TestFixture::FloatType coords[3];
    model.vertexCoords(i, coords);

    EXPECT_GE(coords[0], 0.0);
    EXPECT_LE(coords[0], 1.0);
    EXPECT_GE(coords[1], 0.0);
    EXPECT_LE(coords[1], 1.0);
    EXPECT_GE(coords[2], 0.0);
    EXPECT_LE(coords[2], 1.0);
  }
}

// ============================================================================
// Material Property Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, VelocityAndDensityPositive)
{
  auto model = this->createTestMesh();

  EXPECT_GT(model.getModelVpOnElement(0), 0.0);
  EXPECT_GT(model.getModelRhoOnElement(0), 0.0);
}

TYPED_TEST(ModelUnstructTest, IsModelOnNodesConfig)
{
  auto model = this->createTestMesh();
  EXPECT_FALSE(model.isModelOnNodes());
}

TYPED_TEST(ModelUnstructTest, IsElasticConfig)
{
  auto model = this->createTestMesh();
  EXPECT_FALSE(model.isElastic());
}

// ============================================================================
// Geometry Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, DomainSizeCorrect)
{
  auto model = this->createTestMesh();

  EXPECT_FLOAT_EQ(model.domainSize(0), 1.0);
  EXPECT_FLOAT_EQ(model.domainSize(1), 1.0);
  EXPECT_FLOAT_EQ(model.domainSize(2), 1.0);
}

TYPED_TEST(ModelUnstructTest, MinSpacingPositive)
{
  auto model = this->createTestMesh();
  EXPECT_GT(model.getMinSpacing(), 0.0);
}

TYPED_TEST(ModelUnstructTest, MaxSpeedPositive)
{
  auto model = this->createTestMesh();
  EXPECT_GT(model.getMaxSpeed(), 0.0);
}

// ============================================================================
// Face Normal Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, FaceNormalNormalized)
{
  auto model = this->createTestMesh();

  typename TestFixture::FloatType normal[3];
  model.faceNormal(0, CubicFace::kXMinus, normal);

  typename TestFixture::FloatType norm = std::sqrt(
      normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

  EXPECT_NEAR(norm, 1.0, 1e-5);
}

TYPED_TEST(ModelUnstructTest, OppositeFacesHaveOppositeNormals)
{
  auto model = this->createTestMesh();

  typename TestFixture::FloatType n_xminus[3], n_xplus[3];
  model.faceNormal(0, CubicFace::kXMinus, n_xminus);
  model.faceNormal(0, CubicFace::kXPlus, n_xplus);

  // Should be opposite (dot product ≈ -1)
  typename TestFixture::FloatType dot = n_xminus[0] * n_xplus[0] +
                                        n_xminus[1] * n_xplus[1] +
                                        n_xminus[2] * n_xplus[2];

  EXPECT_NEAR(dot, -1.0, 1e-5);
}

// ============================================================================
// Integration with FaceConnectivity
// ============================================================================

TYPED_TEST(ModelUnstructTest, BuildFaceConnectivityDoesNotCrash)
{
  auto model = this->createTestMesh();
  model.buildFaceConnectivity();
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, FaceConnectivityMethodsAccessible)
{
  auto model = this->createTestMesh();
  model.buildFaceConnectivity();

  EXPECT_EQ(model.getNumberOfFaces(), 6);
  EXPECT_TRUE(model.isBoundaryFace(0));

  auto face = model.getGlobalFace(0, CubicFace::kXMinus);
  EXPECT_GE(face, 0);

  auto node = model.getGlobalNodeFromFace(face, 0);
  EXPECT_GE(node, 0);
  EXPECT_LT(node, 8);
}

}  // namespace
}  // namespace model
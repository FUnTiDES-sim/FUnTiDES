#include <gtest/gtest.h>

#include <cmath>
#include <set>

#include "model_unstruct.h"

namespace model {
namespace {

// ============================================================================
// Test Fixture
// ============================================================================

template <typename T>
class ModelUnstructTest : public ::testing::Test {
 protected:
  using FloatType = typename T::FloatType;
  using ScalarType = typename T::ScalarType;
  using ModelType = ModelUnstruct<FloatType, ScalarType>;

  static void fillGeometry(ModelUnstructData<FloatType, ScalarType>& data) {
    for (int i = 0; i < 8; ++i) {
      data.global_node_index_(0, i) = i;
      data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
      data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
      data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
    }
  }

  // Helper to create minimal test mesh (acoustic, element-based)
  ModelType createTestMesh() {
    ModelUnstructData<FloatType, ScalarType> data;
    data.order_ = 1;
    data.n_element_ = 1;
    data.n_node_ = 8;
    data.lx_ = data.ly_ = data.lz_ = 1.0;
    data.isModelOnNodes_ = false;
    data.isElastic_ = false;

    data.global_node_index_ = allocateArray2D<arrayInt>(1, 8);
    data.nodes_coords_x_ = allocateVector<vectorReal>(8);
    data.nodes_coords_y_ = allocateVector<vectorReal>(8);
    data.nodes_coords_z_ = allocateVector<vectorReal>(8);
    data.model_vp_element_ = allocateVector<vectorReal>(1);
    data.model_rho_element_ = allocateVector<vectorReal>(1);
    data.boundaries_t_ = allocateVector<vectorInt>(8);

    fillGeometry(data);
    data.model_vp_element_[0] = 1500.0;
    data.model_rho_element_[0] = 1.0;

    return ModelType(data);
  }

  // Per-node acoustic mesh — includes vs for setModelNodeProps
  ModelType createNodeMesh() {
    ModelUnstructData<FloatType, ScalarType> data;
    data.order_ = 1;
    data.n_element_ = 1;
    data.n_node_ = 8;
    data.lx_ = data.ly_ = data.lz_ = 1.0;
    data.isModelOnNodes_ = true;
    data.isElastic_ = false;

    data.global_node_index_ = allocateArray2D<arrayInt>(1, 8);
    data.nodes_coords_x_ = allocateVector<vectorReal>(8);
    data.nodes_coords_y_ = allocateVector<vectorReal>(8);
    data.nodes_coords_z_ = allocateVector<vectorReal>(8);
    data.model_vp_node_ = allocateVector<vectorReal>(8);
    data.model_rho_node_ = allocateVector<vectorReal>(8);
    data.model_vs_node_ = allocateVector<vectorReal>(8);
    data.boundaries_t_ = allocateVector<vectorInt>(8);

    fillGeometry(data);
    for (int i = 0; i < 8; ++i) {
      data.model_vp_node_[i] = 1500.0;
      data.model_rho_node_[i] = 1.0;
      data.model_vs_node_[i] = 0.0;
      // Mark node 7 (top corner, z=1) as Surface for boundary tests
      data.boundaries_t_[i] =
          (i == 7) ? static_cast<int>(BoundaryFlag::Surface) : static_cast<int>(BoundaryFlag::InteriorNode);
    }
    return ModelType(data);
  }

  // Per-element elastic mesh with all Thomsen parameters
  ModelType createElasticMesh() {
    ModelUnstructData<FloatType, ScalarType> data;
    data.order_ = 1;
    data.n_element_ = 1;
    data.n_node_ = 8;
    data.lx_ = data.ly_ = data.lz_ = 1.0;
    data.isModelOnNodes_ = false;
    data.isElastic_ = true;

    data.global_node_index_ = allocateArray2D<arrayInt>(1, 8);
    data.nodes_coords_x_ = allocateVector<vectorReal>(8);
    data.nodes_coords_y_ = allocateVector<vectorReal>(8);
    data.nodes_coords_z_ = allocateVector<vectorReal>(8);
    data.model_vp_element_ = allocateVector<vectorReal>(1);
    data.model_rho_element_ = allocateVector<vectorReal>(1);
    data.model_vs_element_ = allocateVector<vectorReal>(1);
    data.model_delta_element_ = allocateVector<vectorReal>(1);
    data.model_epsilon_element_ = allocateVector<vectorReal>(1);
    data.model_gamma_element_ = allocateVector<vectorReal>(1);
    data.model_theta_element_ = allocateVector<vectorReal>(1);
    data.model_phi_element_ = allocateVector<vectorReal>(1);
    data.boundaries_t_ = allocateVector<vectorInt>(8);

    fillGeometry(data);
    data.model_vp_element_[0] = 3000.0;
    data.model_rho_element_[0] = 2000.0;
    data.model_vs_element_[0] = 1500.0;
    data.model_delta_element_[0] = 0.0;
    data.model_epsilon_element_[0] = 0.0;
    data.model_gamma_element_[0] = 0.0;
    data.model_theta_element_[0] = 0.0;
    data.model_phi_element_[0] = 0.0;
    return ModelType(data);
  }
};

struct FloatInt {
  using FloatType = float;
  using ScalarType = int;
};
struct DoubleInt {
  using FloatType = double;
  using ScalarType = int;
};

using TestTypes = ::testing::Types<FloatInt, DoubleInt>;
TYPED_TEST_SUITE(ModelUnstructTest, TestTypes);

// ============================================================================
// Construction Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, DefaultConstructor) {
  typename TestFixture::ModelType model;
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, ConstructFromData) {
  auto model = this->createTestMesh();
  EXPECT_EQ(model.getNumberOfElements(), 1);
  EXPECT_EQ(model.getNumberOfNodes(), 8);
  EXPECT_EQ(model.getOrder(), 1);
}

// ============================================================================
// Indexing Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, GlobalNodeIndexInRange) {
  auto model = this->createTestMesh();

  for (int k = 0; k <= 1; ++k)
    for (int j = 0; j <= 1; ++j)
      for (int i = 0; i <= 1; ++i) {
        auto node = model.globalNodeIndex(0, i, j, k);
        EXPECT_GE(node, 0);
        EXPECT_LT(node, 8);
      }
}

TYPED_TEST(ModelUnstructTest, VertexCoordsInUnitCube) {
  auto model = this->createTestMesh();

  for (int i = 0; i < 8; ++i) {
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

TYPED_TEST(ModelUnstructTest, VelocityAndDensityPositive) {
  auto model = this->createTestMesh();

  EXPECT_GT(model.getModelVpOnElement(0), 0.0);
  EXPECT_GT(model.getModelRhoOnElement(0), 0.0);
}

TYPED_TEST(ModelUnstructTest, IsModelOnNodesConfig) {
  auto model = this->createTestMesh();
  EXPECT_FALSE(model.isModelOnNodes());
}

TYPED_TEST(ModelUnstructTest, IsElasticConfig) {
  auto model = this->createTestMesh();
  EXPECT_FALSE(model.isElastic());
}

// ============================================================================
// Geometry Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, DomainSizeCorrect) {
  auto model = this->createTestMesh();

  EXPECT_FLOAT_EQ(model.domainSize(0), 1.0);
  EXPECT_FLOAT_EQ(model.domainSize(1), 1.0);
  EXPECT_FLOAT_EQ(model.domainSize(2), 1.0);
}

TYPED_TEST(ModelUnstructTest, MinSpacingPositive) {
  auto model = this->createTestMesh();
  EXPECT_GT(model.getMinSpacing(), 0.0);
}

TYPED_TEST(ModelUnstructTest, MaxSpeedPositive) {
  auto model = this->createTestMesh();
  EXPECT_GT(model.getMaxSpeed(), 0.0);
}

// ============================================================================
// Face Normal Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, FaceNormalNormalized) {
  auto model = this->createTestMesh();

  typename TestFixture::FloatType normal[3];
  model.faceNormal(0, CubicFace::kXMinus, normal);

  typename TestFixture::FloatType norm =
      std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

  EXPECT_NEAR(norm, 1.0, 1e-5);
}

TYPED_TEST(ModelUnstructTest, OppositeFacesHaveOppositeNormals) {
  auto model = this->createTestMesh();

  typename TestFixture::FloatType n_xminus[3], n_xplus[3];
  model.faceNormal(0, CubicFace::kXMinus, n_xminus);
  model.faceNormal(0, CubicFace::kXPlus, n_xplus);

  typename TestFixture::FloatType dot = n_xminus[0] * n_xplus[0] + n_xminus[1] * n_xplus[1] + n_xminus[2] * n_xplus[2];

  EXPECT_NEAR(dot, -1.0, 1e-5);
}

// ============================================================================
// All 6 face normals
// ============================================================================

TYPED_TEST(ModelUnstructTest, AllFaceNormalsNormalized) {
  auto model = this->createTestMesh();
  typename TestFixture::FloatType n[3];
  for (int lf = 0; lf < 6; ++lf) {
    model.faceNormal(0, static_cast<CubicFace>(lf), n);
    typename TestFixture::FloatType norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    EXPECT_NEAR(norm, 1.0, 1e-5) << "face " << lf;
  }
}

TYPED_TEST(ModelUnstructTest, YFacesHaveOppositeNormals) {
  auto model = this->createTestMesh();
  typename TestFixture::FloatType ny[3], py[3];
  model.faceNormal(0, CubicFace::kYMinus, ny);
  model.faceNormal(0, CubicFace::kYPlus, py);
  typename TestFixture::FloatType dot = ny[0] * py[0] + ny[1] * py[1] + ny[2] * py[2];
  EXPECT_NEAR(dot, -1.0, 1e-5);
}

TYPED_TEST(ModelUnstructTest, ZFacesHaveOppositeNormals) {
  auto model = this->createTestMesh();
  typename TestFixture::FloatType nz[3], pz[3];
  model.faceNormal(0, CubicFace::kZMinus, nz);
  model.faceNormal(0, CubicFace::kZPlus, pz);
  typename TestFixture::FloatType dot = nz[0] * pz[0] + nz[1] * pz[1] + nz[2] * pz[2];
  EXPECT_NEAR(dot, -1.0, 1e-5);
}

// ============================================================================
// elementIndex / globalVertexIndex / getNumberOfPointsPerElement
// ============================================================================

TYPED_TEST(ModelUnstructTest, ElementIndexIsIdentity) {
  auto model = this->createTestMesh();
  for (int i = 0; i < 5; ++i) EXPECT_EQ(model.elementIndex(i), i);
}

TYPED_TEST(ModelUnstructTest, GlobalVertexIndexInRange) {
  auto model = this->createTestMesh();
  for (int i = 0; i <= 1; ++i)
    for (int j = 0; j <= 1; ++j)
      for (int k = 0; k <= 1; ++k) {
        auto v = model.globalVertexIndex(0, i, j, k);
        EXPECT_GE(v, 0);
        EXPECT_LT(v, 8);
      }
}

TYPED_TEST(ModelUnstructTest, NumberOfPointsPerElementOrder1) {
  auto model = this->createTestMesh();
  EXPECT_EQ(model.getNumberOfPointsPerElement(), 8);  // (1+1)^3
}

// ============================================================================
// nodeCoord / domainSize invalid dim
// ============================================================================

TYPED_TEST(ModelUnstructTest, NodeCoordInvalidDimReturnsMinusOne) {
  auto model = this->createTestMesh();
  EXPECT_FLOAT_EQ(model.nodeCoord(0, 3), -1.0f);
}

TYPED_TEST(ModelUnstructTest, DomainSizeInvalidDimReturnsMinusOne) {
  auto model = this->createTestMesh();
  EXPECT_FLOAT_EQ(model.domainSize(3), -1.0f);
}

// ============================================================================
// Per-node material getters
// ============================================================================

TYPED_TEST(ModelUnstructTest, VelocityPropertiesOnNodes) {
  auto model = this->createNodeMesh();
  EXPECT_FLOAT_EQ(model.getModelVpOnNodes(0), 1500.0f);
  EXPECT_FLOAT_EQ(model.getModelRhoOnNodes(0), 1.0f);
}

TYPED_TEST(ModelUnstructTest, MaxSpeedViaNodeArray) {
  auto model = this->createNodeMesh();
  EXPECT_FLOAT_EQ(model.getMaxSpeed(), 1500.0f);
}

TYPED_TEST(ModelUnstructTest, SetModelNodePropsUpdatesValues) {
  auto model = this->createNodeMesh();
  model.setModelNodeProps(0, 3000.0f, 1500.0f, 2000.0f);
  EXPECT_FLOAT_EQ(model.getModelVpOnNodes(0), 3000.0f);
  EXPECT_FLOAT_EQ(model.getModelVsOnNodes(0), 1500.0f);
  EXPECT_FLOAT_EQ(model.getModelRhoOnNodes(0), 2000.0f);
  // Other nodes unchanged
  EXPECT_FLOAT_EQ(model.getModelVpOnNodes(1), 1500.0f);
}

// ============================================================================
// Quality factors
// ============================================================================

TYPED_TEST(ModelUnstructTest, QualityFactorDefaultsToLargeValue) {
  auto model = this->createTestMesh();
  EXPECT_FLOAT_EQ(model.getModelQpOnElement(0), 1.0e9f);
  EXPECT_FLOAT_EQ(model.getModelQsOnElement(0), 1.0e9f);
  EXPECT_FLOAT_EQ(model.getModelQpOnNodes(0), 1.0e9f);
  EXPECT_FLOAT_EQ(model.getModelQsOnNodes(0), 1.0e9f);
}

TYPED_TEST(ModelUnstructTest, SetQualityFactorsFillsAllElements) {
  auto model = this->createTestMesh();
  model.setQualityFactors(100.0f, 50.0f);
  EXPECT_FLOAT_EQ(model.getModelQpOnElement(0), 100.0f);
  EXPECT_FLOAT_EQ(model.getModelQsOnElement(0), 50.0f);
}

// ============================================================================
// Boundary type / isFreeSurface
// ============================================================================

TYPED_TEST(ModelUnstructTest, BoundaryTypeWithEmptyBoundariesIsInterior) {
  ModelUnstructData<typename TestFixture::FloatType, typename TestFixture::ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 1;
  data.n_node_ = 8;
  data.lx_ = data.ly_ = data.lz_ = 1.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;
  data.global_node_index_ = allocateArray2D<arrayInt>(1, 8);
  data.nodes_coords_x_ = allocateVector<vectorReal>(8);
  data.nodes_coords_y_ = allocateVector<vectorReal>(8);
  data.nodes_coords_z_ = allocateVector<vectorReal>(8);
  data.model_vp_element_ = allocateVector<vectorReal>(1);
  data.model_rho_element_ = allocateVector<vectorReal>(1);
  // boundaries_t_ left empty on purpose
  TestFixture::fillGeometry(data);
  data.model_vp_element_[0] = 1500.0;
  data.model_rho_element_[0] = 1.0;
  typename TestFixture::ModelType model(data);
  EXPECT_EQ(model.boundaryType(0), BoundaryFlag::InteriorNode);
  EXPECT_FALSE(model.isFreeSurface(0));
}

TYPED_TEST(ModelUnstructTest, IsFreeSurfaceWithSurfaceNode) {
  auto model = this->createNodeMesh();
  // Node 7 marked as Surface in createNodeMesh
  EXPECT_TRUE(model.isFreeSurface(7));
  EXPECT_FALSE(model.isFreeSurface(0));
}

TYPED_TEST(ModelUnstructTest, BoundaryTypeWithSurfaceNode) {
  auto model = this->createNodeMesh();
  EXPECT_EQ(model.boundaryType(7), BoundaryFlag::Surface);
  EXPECT_EQ(model.boundaryType(0), BoundaryFlag::InteriorNode);
}

// ============================================================================
// buildFaceConnectivity + face queries
// ============================================================================

TYPED_TEST(ModelUnstructTest, BuildFaceConnectivitySingleElement) {
  auto model = this->createTestMesh();
  model.buildFaceConnectivity();
  EXPECT_EQ(model.getNumberOfFaces(), 6);
}

TYPED_TEST(ModelUnstructTest, BuildFaceConnectivityIdempotent) {
  auto model = this->createTestMesh();
  model.buildFaceConnectivity();
  auto n1 = model.getNumberOfFaces();
  model.buildFaceConnectivity();  // second call must not re-build
  EXPECT_EQ(model.getNumberOfFaces(), n1);
}

TYPED_TEST(ModelUnstructTest, GetGlobalFaceReturnsValidId) {
  auto model = this->createTestMesh();
  model.buildFaceConnectivity();
  std::set<int> face_ids;
  for (int lf = 0; lf < 6; ++lf) {
    auto fid = model.getGlobalFace(0, static_cast<CubicFace>(lf));
    EXPECT_GE(fid, 0);
    EXPECT_LT(fid, 6);
    face_ids.insert(fid);
  }
  EXPECT_EQ(face_ids.size(), 6u);
}

TYPED_TEST(ModelUnstructTest, AllFacesBoundaryInSingleElement) {
  auto model = this->createTestMesh();
  model.buildFaceConnectivity();
  for (int lf = 0; lf < 6; ++lf) {
    auto fid = model.getGlobalFace(0, static_cast<CubicFace>(lf));
    EXPECT_TRUE(model.isBoundaryFace(fid));
  }
}

TYPED_TEST(ModelUnstructTest, GetGlobalNodeFromFaceReturnsValidNode) {
  auto model = this->createTestMesh();
  model.buildFaceConnectivity();
  auto fid = model.getGlobalFace(0, CubicFace::kXMinus);
  int n_dofs = (1 + 1) * (1 + 1);  // order=1 → 4
  for (int dof = 0; dof < n_dofs; ++dof) {
    auto node = model.getGlobalNodeFromFace(fid, dof);
    EXPECT_GE(node, 0);
    EXPECT_LT(node, 8);
  }
}

// ============================================================================
// Elastic material properties
// ============================================================================

TYPED_TEST(ModelUnstructTest, ElasticMaterialVsAndThomsens) {
  auto model = this->createElasticMesh();
  EXPECT_FLOAT_EQ(model.getModelVsOnElement(0), 1500.0f);
  EXPECT_FLOAT_EQ(model.getModelDeltaOnElement(0), 0.0f);
  EXPECT_FLOAT_EQ(model.getModelEpsilonOnElement(0), 0.0f);
  EXPECT_FLOAT_EQ(model.getModelGammaOnElement(0), 0.0f);
  EXPECT_EQ(model.getModelThetaOnElement(0), 0);
  EXPECT_EQ(model.getModelPhiOnElement(0), 0);
}

// ============================================================================
// initElasticityTensors
// ============================================================================

TYPED_TEST(ModelUnstructTest, InitElasticityTensorsNonElasticReturnsEarly) {
  auto model = this->createTestMesh();
  model.initElasticityTensors(kIso);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, InitElasticityTensorsIsoReturnsEarly) {
  auto model = this->createElasticMesh();
  model.initElasticityTensors(kIso);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, InitElasticityTensorsVTIReturnsEarly) {
  auto model = this->createElasticMesh();
  model.initElasticityTensors(kVTI);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, InitElasticityTensorsTTISymmetric) {
  auto model = this->createElasticMesh();
  model.initElasticityTensors(kTTI);
  Kokkos::fence();
  typename TestFixture::FloatType C[6][6];
  model.getCTensorOnElement(0, C);
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) EXPECT_NEAR(C[i][j], C[j][i], 1e-4f) << "C[" << i << "][" << j << "]";
}

}  // namespace
}  // namespace model
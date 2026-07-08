#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <set>

#include "model_struct.h"

// Highest polynomial order compiled on this machine. Provided by CMake via
// MAX_SOLVER_ORDER (see tests/units/model/impl/mesh/struct/CMakeLists.txt).
// Fall back to the full implementation range when the build does not restrict
// it, so environments compiling all orders still exercise them.
#ifndef MAX_SOLVER_ORDER
#define MAX_SOLVER_ORDER 9
#endif

namespace model {
namespace {

// ============================================================================
// Type Wrapper Classes for Non-Type Template Parameters
// ============================================================================
// Google Test's TYPED_TEST_SUITE requires types, not non-type parameters.
// We wrap template instantiations in structs to make them compatible.

// Order 1 wrappers
struct FloatOrder1 {
  using FloatType = float;
  using ScalarType = int;
  static constexpr int Order = 1;
};

struct DoubleOrder1 {
  using FloatType = double;
  using ScalarType = int;
  static constexpr int Order = 1;
};

// Order 2 wrappers
struct FloatOrder2 {
  using FloatType = float;
  using ScalarType = int;
  static constexpr int Order = 2;
};

struct DoubleOrder2 {
  using FloatType = double;
  using ScalarType = int;
  static constexpr int Order = 2;
};

// Order 3 wrappers
struct FloatOrder3 {
  using FloatType = float;
  using ScalarType = int;
  static constexpr int Order = 3;
};

struct DoubleOrder3 {
  using FloatType = double;
  using ScalarType = int;
  static constexpr int Order = 3;
};

// Order 4 wrappers
#if MAX_SOLVER_ORDER >= 4
struct FloatOrder4 {
  using FloatType = float;
  using ScalarType = int;
  static constexpr int Order = 4;
};

struct DoubleOrder4 {
  using FloatType = double;
  using ScalarType = int;
  static constexpr int Order = 4;
};
#endif

// Order 5 wrappers
#if MAX_SOLVER_ORDER >= 5
struct FloatOrder5 {
  using FloatType = float;
  using ScalarType = int;
  static constexpr int Order = 5;
};

struct DoubleOrder5 {
  using FloatType = double;
  using ScalarType = int;
  static constexpr int Order = 5;
};
#endif

// ============================================================================
// Test Fixture
// ============================================================================

template <typename TypeWrapper>
class ModelStructTest : public ::testing::Test {
 protected:
  using FloatType = typename TypeWrapper::FloatType;
  using ScalarType = typename TypeWrapper::ScalarType;
  static constexpr int Order = TypeWrapper::Order;

  using ModelStructType = ModelStruct<FloatType, ScalarType, Order>;
  using ModelStructDataType = ModelStructData<FloatType, ScalarType>;

  void SetUp() override {
    // Standard structured mesh: 10x10x10 elements
    // With Order=1: 11x11x11 nodes
    // With Order=2: 21x21x21 nodes
    // With Order=3: 31x31x31 nodes
    data_.ex_ = 10;
    data_.ey_ = 10;
    data_.ez_ = 10;
    data_.dx_ = 100.0;
    data_.dy_ = 100.0;
    data_.dz_ = 100.0;
    data_.ox_ = 0.0;
    data_.oy_ = 0.0;
    data_.oz_ = 0.0;
    data_.isModelOnNodes_ = true;
    data_.isElastic_ = false;

    model_ = std::make_unique<ModelStructType>(data_);
  }

  ModelStructDataType data_;
  std::unique_ptr<ModelStructType> model_;
};

// Register type wrappers for typed tests. Orders 4 and 5 are only included
// when they are compiled on the current machine (see MAX_SOLVER_ORDER above).
#if MAX_SOLVER_ORDER >= 5
using TypeWrappers = ::testing::Types<FloatOrder1, FloatOrder2, FloatOrder3, FloatOrder4, FloatOrder5, DoubleOrder1,
                                      DoubleOrder2, DoubleOrder3, DoubleOrder4, DoubleOrder5>;
#elif MAX_SOLVER_ORDER == 4
using TypeWrappers = ::testing::Types<FloatOrder1, FloatOrder2, FloatOrder3, FloatOrder4, DoubleOrder1, DoubleOrder2,
                                      DoubleOrder3, DoubleOrder4>;
#else
using TypeWrappers =
    ::testing::Types<FloatOrder1, FloatOrder2, FloatOrder3, DoubleOrder1, DoubleOrder2, DoubleOrder3>;
#endif

TYPED_TEST_SUITE(ModelStructTest, TypeWrappers);

// ============================================================================
// Constructor and Initialization Tests
// ============================================================================

TYPED_TEST(ModelStructTest, DefaultConstructor) {
  using ModelStructType = typename TestFixture::ModelStructType;
  ModelStructType model;
  // Default constructor should compile and not crash
  SUCCEED();
}

TYPED_TEST(ModelStructTest, ConstructorFromData) {
  auto& model = *this->model_;
  EXPECT_EQ(model.getNumberOfElements(), 1000);  // 10*10*10
}

TYPED_TEST(ModelStructTest, CopyConstructor) {
  auto& original = *this->model_;
  auto copy = original;
  EXPECT_EQ(copy.getNumberOfElements(), original.getNumberOfElements());
  EXPECT_EQ(copy.getNumberOfNodes(), original.getNumberOfNodes());
}

TYPED_TEST(ModelStructTest, AssignmentOperator) {
  auto original = *this->model_;
  auto assigned = original;
  EXPECT_EQ(assigned.getNumberOfElements(), original.getNumberOfElements());
}

// ============================================================================
// Mesh Topology Tests
// ============================================================================

TYPED_TEST(ModelStructTest, ElementIndexConversion) {
  auto& model = *this->model_;

  // Test conversion from linear to 3D indices and back
  for (int linear = 0; linear < 1000; ++linear) {
    auto idx = model.elementIndex(linear);
    EXPECT_GE(idx[0], 0);
    EXPECT_LT(idx[0], 10);
    EXPECT_GE(idx[1], 0);
    EXPECT_LT(idx[1], 10);
    EXPECT_GE(idx[2], 0);
    EXPECT_LT(idx[2], 10);
  }
}

TYPED_TEST(ModelStructTest, ElementIndexBoundaries) {
  auto& model = *this->model_;

  // First element
  auto first = model.elementIndex(0);
  EXPECT_EQ(first[0], 0);
  EXPECT_EQ(first[1], 0);
  EXPECT_EQ(first[2], 0);

  // Last element
  auto last = model.elementIndex(999);
  EXPECT_EQ(last[0], 9);
  EXPECT_EQ(last[1], 9);
  EXPECT_EQ(last[2], 9);
}

TYPED_TEST(ModelStructTest, GlobalVertexIndex) {
  auto& model = *this->model_;

  auto elem_idx = std::array<int, 3>{5, 5, 5};
  for (int i = 0; i <= 1; ++i) {
    for (int j = 0; j <= 1; ++j) {
      for (int k = 0; k <= 1; ++k) {
        auto vertex = model.globalVertexIndex(elem_idx, i, j, k);
        EXPECT_GE(vertex[0], 5);
        EXPECT_LE(vertex[0], 6);
        EXPECT_GE(vertex[1], 5);
        EXPECT_LE(vertex[1], 6);
        EXPECT_GE(vertex[2], 5);
        EXPECT_LE(vertex[2], 6);
      }
    }
  }
}

// ============================================================================
// Coordinate System Tests
// ============================================================================

TYPED_TEST(ModelStructTest, VertexCoordinates) {
  auto& model = *this->model_;
  typename TestFixture::FloatType coords[3];

  // Test all 8 vertices of the first element (cube from (0,0,0) to (10,10,10))
  // Each element has size 10x10x10 since domain is 100x100x100 with 10 elements
  // per dimension

  // Vertex 0: (0, 0, 0)
  auto dof_global0 = std::array<int, 3>{0, 0, 0};
  model.vertexCoords(dof_global0, coords);
  EXPECT_EQ(coords[0], 0.0);
  EXPECT_EQ(coords[1], 0.0);
  EXPECT_EQ(coords[2], 0.0);

  // Vertex 1: (1, 0, 0) -> (10, 0, 0)
  auto dof_global1 = std::array<int, 3>{1, 0, 0};
  model.vertexCoords(dof_global1, coords);
  EXPECT_EQ(coords[0], 10.0);
  EXPECT_EQ(coords[1], 0.0);
  EXPECT_EQ(coords[2], 0.0);

  // Vertex 2: (0, 1, 0) -> (0, 10, 0)
  auto dof_global2 = std::array<int, 3>{0, 1, 0};
  model.vertexCoords(dof_global2, coords);
  EXPECT_EQ(coords[0], 0.0);
  EXPECT_EQ(coords[1], 10.0);
  EXPECT_EQ(coords[2], 0.0);

  // Vertex 3: (1, 1, 0) -> (10, 10, 0)
  auto dof_global3 = std::array<int, 3>{1, 1, 0};
  model.vertexCoords(dof_global3, coords);
  EXPECT_EQ(coords[0], 10.0);
  EXPECT_EQ(coords[1], 10.0);
  EXPECT_EQ(coords[2], 0.0);

  // Vertex 4: (0, 0, 1) -> (0, 0, 10)
  auto dof_global4 = std::array<int, 3>{0, 0, 1};
  model.vertexCoords(dof_global4, coords);
  EXPECT_EQ(coords[0], 0.0);
  EXPECT_EQ(coords[1], 0.0);
  EXPECT_EQ(coords[2], 10.0);

  // Vertex 5: (1, 0, 1) -> (10, 0, 10)
  auto dof_global5 = std::array<int, 3>{1, 0, 1};
  model.vertexCoords(dof_global5, coords);
  EXPECT_EQ(coords[0], 10.0);
  EXPECT_EQ(coords[1], 0.0);
  EXPECT_EQ(coords[2], 10.0);

  // Vertex 6: (0, 1, 1) -> (0, 10, 10)
  auto dof_global6 = std::array<int, 3>{0, 1, 1};
  model.vertexCoords(dof_global6, coords);
  EXPECT_EQ(coords[0], 0.0);
  EXPECT_EQ(coords[1], 10.0);
  EXPECT_EQ(coords[2], 10.0);

  // Vertex 7: (1, 1, 1) -> (10, 10, 10)
  auto dof_global7 = std::array<int, 3>{1, 1, 1};
  model.vertexCoords(dof_global7, coords);
  EXPECT_EQ(coords[0], 10.0);
  EXPECT_EQ(coords[1], 10.0);
  EXPECT_EQ(coords[2], 10.0);
}

TYPED_TEST(ModelStructTest, NodeCoordAtOrigin) {
  auto& model = *this->model_;

  auto x = model.nodeCoord(0, 0);
  auto y = model.nodeCoord(0, 1);
  auto z = model.nodeCoord(0, 2);

  EXPECT_NEAR(x, 0.0, 1e-6);
  EXPECT_NEAR(y, 0.0, 1e-6);
  EXPECT_NEAR(z, 0.0, 1e-6);
}

TYPED_TEST(ModelStructTest, NodeCoordBoundaryX) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;

  // Number of nodes in x direction
  auto nx = Order * 10 + 1;

  // Last node in x direction, first in y and z
  auto last_node_idx = nx - 1;
  auto x = model.nodeCoord(last_node_idx, 0);

  EXPECT_NEAR(x, 100.0, 1e-6);
}

TYPED_TEST(ModelStructTest, NodeCoordMonotonicallyIncreasing) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;

  auto nx = Order * 10 + 1;

  // Sample nodes along x-direction
  auto prev_x = model.nodeCoord(0, 0);
  for (auto i = 1; i < nx; ++i) {
    auto curr_x = model.nodeCoord(i, 0);
    EXPECT_GT(curr_x, prev_x);
    prev_x = curr_x;
  }
}

// ============================================================================
// Global Node Index Tests
// ============================================================================

TYPED_TEST(ModelStructTest, GlobalNodeIndexFirstElement) {
  auto& model = *this->model_;

  // First element at origin
  auto node_idx = model.globalNodeIndex(0, 0, 0, 0);
  EXPECT_EQ(node_idx, 0);
}

TYPED_TEST(ModelStructTest, GlobalNodeIndexConsistency) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;

  // All nodes in first element should be within valid range
  for (int i = 0; i <= Order; ++i) {
    for (int j = 0; j <= Order; ++j) {
      for (int k = 0; k <= Order; ++k) {
        auto node_idx = model.globalNodeIndex(0, i, j, k);
        EXPECT_GE(node_idx, 0);
        EXPECT_LT(node_idx, model.getNumberOfNodes());
      }
    }
  }
}

// ============================================================================
// Element and Node Count Tests
// ============================================================================

TYPED_TEST(ModelStructTest, NumberOfElements) {
  auto& model = *this->model_;
  EXPECT_EQ(model.getNumberOfElements(), 1000);
}

TYPED_TEST(ModelStructTest, NumberOfNodes) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;

  auto expected_nodes = (Order * 10 + 1) * (Order * 10 + 1) * (Order * 10 + 1);
  EXPECT_EQ(model.getNumberOfNodes(), expected_nodes);
}

TYPED_TEST(ModelStructTest, NumberOfPointsPerElement) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;

  int expected_points = (Order + 1) * (Order + 1) * (Order + 1);
  EXPECT_EQ(model.getNumberOfPointsPerElement(), expected_points);
}

TYPED_TEST(ModelStructTest, GetOrder) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;

  EXPECT_EQ(model.getOrder(), Order);
}

// ============================================================================
// Material Property Tests
// ============================================================================

TYPED_TEST(ModelStructTest, VelocityPropertiesOnNodes) {
  auto& model = *this->model_;

  // Test Vp on first few nodes
  for (int i = 0; i < 10; ++i) {
    auto vp = model.getModelVpOnNodes(i);
    EXPECT_EQ(vp, 1500.0);
  }
}

TYPED_TEST(ModelStructTest, VelocityPropertiesOnElements) {
  auto& model = *this->model_;

  // Test on first few elements
  for (int i = 0; i < 10; ++i) {
    auto vp = model.getModelVpOnElement(i);
    auto vs = model.getModelVsOnElement(i);
    auto rho = model.getModelRhoOnElement(i);

    EXPECT_EQ(vp, 1500.0);
    EXPECT_EQ(vs, 755.0);
    EXPECT_EQ(rho, 1.0);
  }
}

TYPED_TEST(ModelStructTest, AnisotropicAnglesOnElements) {
  auto& model = *this->model_;

  // Test on first element
  auto theta = model.getModelThetaOnElement(0);
  auto phi = model.getModelPhiOnElement(0);

  EXPECT_EQ(theta, 0.0);
  EXPECT_EQ(phi, 0.0);
}

TYPED_TEST(ModelStructTest, DensityPropertiesConsistent) {
  auto& model = *this->model_;

  auto rho_node = model.getModelRhoOnNodes(0);
  auto rho_elem = model.getModelRhoOnElement(0);

  EXPECT_EQ(rho_node, rho_elem);
  EXPECT_EQ(rho_node, 1.0);
}

// ============================================================================
// Domain Properties Tests
// ============================================================================

TYPED_TEST(ModelStructTest, DomainSize) {
  auto& model = *this->model_;

  auto lx = model.domainSize(0);
  auto ly = model.domainSize(1);
  auto lz = model.domainSize(2);

  EXPECT_EQ(lx, 100.0);
  EXPECT_EQ(ly, 100.0);
  EXPECT_EQ(lz, 100.0);
}

TYPED_TEST(ModelStructTest, DomainSizeInvalidDimension) {
  auto& model = *this->model_;

  // Invalid dimension should return -1
  EXPECT_EQ(model.domainSize(3), -1);
  EXPECT_EQ(model.domainSize(-1), -1);
}

TYPED_TEST(ModelStructTest, MinSpacingPositive) {
  auto& model = *this->model_;

  auto min_spacing = model.getMinSpacing();
  EXPECT_GT(min_spacing, 0.0);
}

TYPED_TEST(ModelStructTest, MinSpacingDependsOnOrder) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;

  auto min_spacing = model.getMinSpacing();

  if constexpr (Order == 1) {
    EXPECT_EQ(min_spacing, 10.0);  // min(10, 10, 10)
  } else if constexpr (Order == 2) {
    EXPECT_EQ(min_spacing, 5.0);  // min(10, 10, 10) / 2
  }
  // Higher orders have specific coefficients
}

TYPED_TEST(ModelStructTest, MaxSpeedPositive) {
  auto& model = *this->model_;

  auto max_speed = model.getMaxSpeed();
  EXPECT_GT(max_speed, 0.0);
}

// ============================================================================
// Model Configuration Tests
// ============================================================================

TYPED_TEST(ModelStructTest, IsModelOnNodes) {
  auto& model = *this->model_;
  EXPECT_TRUE(model.isModelOnNodes());
}

TYPED_TEST(ModelStructTest, IsElastic) {
  auto& model = *this->model_;
  EXPECT_FALSE(model.isElastic());
}

TYPED_TEST(ModelStructTest, IsElasticInitialization) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;

  ModelStructDataType elastic_data;
  elastic_data.ex_ = 5;
  elastic_data.ey_ = 5;
  elastic_data.ez_ = 5;
  elastic_data.dx_ = 50.0;
  elastic_data.dy_ = 50.0;
  elastic_data.dz_ = 50.0;
  elastic_data.isElastic_ = true;
  elastic_data.isModelOnNodes_ = true;

  ModelStructType elastic_model(elastic_data);
  EXPECT_TRUE(elastic_model.isElastic());
}

// ============================================================================
// Boundary Tests
// ============================================================================

TYPED_TEST(ModelStructTest, FaceNormal) {
  auto& model = *this->model_;
  typename TestFixture::FloatType normal[3] = {0.0, 0.0, 0.0};

  // Test XPlus face normal
  model.faceNormal(0, CubicFace::kXPlus, normal);
  EXPECT_NEAR(normal[0], 1.0, 1e-6);
  EXPECT_NEAR(normal[1], 0.0, 1e-6);
  EXPECT_NEAR(normal[2], 0.0, 1e-6);

  // Test XMinus face normal
  model.faceNormal(0, CubicFace::kXMinus, normal);
  EXPECT_NEAR(normal[0], -1.0, 1e-6);
  EXPECT_NEAR(normal[1], 0.0, 1e-6);
  EXPECT_NEAR(normal[2], 0.0, 1e-6);

  // Test YPlus face normal
  model.faceNormal(0, CubicFace::kYPlus, normal);
  EXPECT_NEAR(normal[0], 0.0, 1e-6);
  EXPECT_NEAR(normal[1], 1.0, 1e-6);
  EXPECT_NEAR(normal[2], 0.0, 1e-6);

  // Test YMinus face normal
  model.faceNormal(0, CubicFace::kYMinus, normal);
  EXPECT_NEAR(normal[0], 0.0, 1e-6);
  EXPECT_NEAR(normal[1], -1.0, 1e-6);
  EXPECT_NEAR(normal[2], 0.0, 1e-6);

  // Test ZPlus face normal
  model.faceNormal(0, CubicFace::kZPlus, normal);
  EXPECT_NEAR(normal[0], 0.0, 1e-6);
  EXPECT_NEAR(normal[1], 0.0, 1e-6);
  EXPECT_NEAR(normal[2], 1.0, 1e-6);

  // Test ZMinus face normal
  model.faceNormal(0, CubicFace::kZMinus, normal);
  EXPECT_NEAR(normal[0], 0.0, 1e-6);
  EXPECT_NEAR(normal[1], 0.0, 1e-6);
  EXPECT_NEAR(normal[2], -1.0, 1e-6);
}

// ============================================================================
// Elasticity Tensor Tests - UPDATED
// ============================================================================

TYPED_TEST(ModelStructTest, InitElasticityTensorsNonElastic) {
  auto& model = *this->model_;

  // Should return early if not elastic
  model.initElasticityTensors(model::kIso);
  SUCCEED();
}

TYPED_TEST(ModelStructTest, InitElasticityTensorsElastic) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;
  using FloatType = typename TestFixture::FloatType;

  ModelStructDataType elastic_data;
  elastic_data.ex_ = 2;
  elastic_data.ey_ = 2;
  elastic_data.ez_ = 2;
  elastic_data.dx_ = 20.0;
  elastic_data.dy_ = 20.0;
  elastic_data.dz_ = 20.0;
  elastic_data.isElastic_ = true;
  elastic_data.isModelOnNodes_ = false;  // CHANGED: Use elements for TTI

  ModelStructType elastic_model(elastic_data);
  elastic_model.initElasticityTensors(model::kTTI);  // CHANGED: Use TTI instead of kIso

  // Verify tensors were created
  FloatType C[6][6];
  elastic_model.getCTensorOnElement(0, C);

  // Verify symmetry of C tensor (Voigt notation)
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      EXPECT_NEAR(C[i][j], C[j][i],
                  1e-6);  // CHANGED: Relaxed tolerance for float
    }
  }
}

TYPED_TEST(ModelStructTest, GetCTensorOnElement) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;
  using FloatType = typename TestFixture::FloatType;

  ModelStructDataType elastic_data;
  elastic_data.ex_ = 3;
  elastic_data.ey_ = 3;
  elastic_data.ez_ = 3;
  elastic_data.dx_ = 30.0;
  elastic_data.dy_ = 30.0;
  elastic_data.dz_ = 30.0;
  elastic_data.isElastic_ = true;
  elastic_data.isModelOnNodes_ = false;  // CHANGED: Use elements for TTI

  ModelStructType elastic_model(elastic_data);
  elastic_model.initElasticityTensors(model::kTTI);  // CHANGED: Use TTI instead of kIso

  FloatType C[6][6];
  elastic_model.getCTensorOnElement(0, C);

  // Verify tensor values are reasonable (non-zero)
  bool has_nonzero = false;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      if (C[i][j] != 0.0) {
        has_nonzero = true;
      }
    }
  }
  EXPECT_TRUE(has_nonzero);
}

// ============================================================================
// Large Mesh Tests
// ============================================================================

TYPED_TEST(ModelStructTest, LargeMeshHandling) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;

  ModelStructDataType large_data;
  large_data.ex_ = 50;
  large_data.ey_ = 50;
  large_data.ez_ = 50;
  large_data.dx_ = 500.0;
  large_data.dy_ = 500.0;
  large_data.dz_ = 500.0;
  large_data.isModelOnNodes_ = true;
  large_data.isElastic_ = false;

  ModelStructType large_model(large_data);

  EXPECT_EQ(large_model.getNumberOfElements(), 125000);
  EXPECT_GT(large_model.getNumberOfNodes(), 125000);
}

TYPED_TEST(ModelStructTest, UniformMeshProperties) {
  auto& model = *this->model_;

  // In a uniform mesh, element spacing should be consistent
  EXPECT_EQ(model.domainSize(0) / 10.0, 10.0);
  EXPECT_EQ(model.domainSize(1) / 10.0, 10.0);
  EXPECT_EQ(model.domainSize(2) / 10.0, 10.0);
}

// ============================================================================
// Documentation Tests
// ============================================================================

TYPED_TEST(ModelStructTest, StructuredVsUnstructuredComparison) { SUCCEED(); }

// ============================================================================
// Face Connectivity Delegate Tests
// ============================================================================

TYPED_TEST(ModelStructTest, BuildFaceConnectivityAndGetNumberOfFaces) {
  auto& model = *this->model_;
  model.buildFaceConnectivity();
  // 10x10x10 mesh:
  // X-faces: (10+1)*10*10 = 1100
  // Y-faces: 10*(10+1)*10 = 1100
  // Z-faces: 10*10*(10+1) = 1100
  // Total: 3300
  EXPECT_EQ(model.getNumberOfFaces(), 3300);
}

TYPED_TEST(ModelStructTest, GetGlobalFaceReturnsValidId) {
  auto& model = *this->model_;
  // All 6 faces of element 0 must be in [0, nFaces)
  auto n_faces = model.getNumberOfFaces();
  for (int lf = 0; lf < 6; ++lf) {
    auto face_id = model.getGlobalFace(0, static_cast<CubicFace>(lf));
    EXPECT_GE(face_id, 0);
    EXPECT_LT(face_id, n_faces);
  }
}

TYPED_TEST(ModelStructTest, GetGlobalFaceUniquenessForElement0) {
  auto& model = *this->model_;
  std::set<int> face_ids;
  for (int lf = 0; lf < 6; ++lf) face_ids.insert(model.getGlobalFace(0, static_cast<CubicFace>(lf)));
  EXPECT_EQ(face_ids.size(), 6u);
}

TYPED_TEST(ModelStructTest, GetGlobalNodeFromFaceReturnsValidNode) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;
  int n_dofs = (Order + 1) * (Order + 1);
  auto face0 = model.getGlobalFace(0, CubicFace::kXMinus);
  for (int dof = 0; dof < n_dofs; ++dof) {
    auto node = model.getGlobalNodeFromFace(face0, dof);
    EXPECT_GE(node, 0);
    EXPECT_LT(node, model.getNumberOfNodes());
  }
}

TYPED_TEST(ModelStructTest, IsBoundaryFaceWithoutBoundariesT) {
  auto& model = *this->model_;
  // No boundaries_t_ set → delegates to FaceConnectivityStruct
  // Face kXMinus of element 0 must be boundary (i=0)
  auto face_xminus_e0 = model.getGlobalFace(0, CubicFace::kXMinus);
  EXPECT_TRUE(model.isBoundaryFace(face_xminus_e0));
  // Shared X-face between elements 0 and 1 must NOT be boundary
  auto shared = model.getGlobalFace(0, CubicFace::kXPlus);
  EXPECT_FALSE(model.isBoundaryFace(shared));
}

TYPED_TEST(ModelStructTest, ElemOwnerAndNeighborConsistency) {
  auto& model = *this->model_;
  // Shared face between elem 0 and elem 1 along X
  auto shared = model.getGlobalFace(0, CubicFace::kXPlus);
  auto owner = model.elemOwner(shared);
  auto neighbor = model.elemNeighbor(shared);
  EXPECT_GE(owner, 0);
  EXPECT_GE(neighbor, 0);
  EXPECT_NE(owner, neighbor);
}

TYPED_TEST(ModelStructTest, ElemNeighborBoundaryReturnsMinusOne) {
  auto& model = *this->model_;
  auto boundary_face = model.getGlobalFace(0, CubicFace::kXMinus);
  EXPECT_EQ(model.elemNeighbor(boundary_face), -1);
}

TYPED_TEST(ModelStructTest, LocalFaceOwnerAndNeighborOpposite) {
  auto& model = *this->model_;
  auto shared = model.getGlobalFace(0, CubicFace::kXPlus);
  int owner_lf = model.localFaceOwner(shared);
  int neigh_lf = model.localFaceNeighbor(shared);
  EXPECT_NE(neigh_lf, -1);
  EXPECT_EQ(owner_lf ^ 1, neigh_lf);
}

TYPED_TEST(ModelStructTest, LocalFaceNeighborBoundaryReturnsMinusOne) {
  auto& model = *this->model_;
  auto boundary_face = model.getGlobalFace(0, CubicFace::kXMinus);
  EXPECT_EQ(model.localFaceNeighbor(boundary_face), -1);
}

// ============================================================================
// Boundary / isFreeSurface / boundaryType
// ============================================================================

TYPED_TEST(ModelStructTest, IsFreeSurfaceWithEmptyBoundariesFalse) {
  auto& model = *this->model_;
  // No boundaries_t_ → always false
  EXPECT_FALSE(model.isFreeSurface(0));
  EXPECT_FALSE(model.isFreeSurface(model.getNumberOfNodes() - 1));
}

TYPED_TEST(ModelStructTest, BoundaryTypeWithEmptyBoundariesIsInterior) {
  auto& model = *this->model_;
  EXPECT_EQ(model.boundaryType(0), BoundaryFlag::InteriorNode);
  EXPECT_EQ(model.boundaryType(model.getNumberOfNodes() - 1), BoundaryFlag::InteriorNode);
}

// ============================================================================
// setQualityFactors / getModelQp / getModelQs
// ============================================================================

TYPED_TEST(ModelStructTest, QualityFactorDefaultsToLargeValue) {
  auto& model = *this->model_;
  EXPECT_FLOAT_EQ(model.getModelQpOnNodes(0), 1.0e9f);
  EXPECT_FLOAT_EQ(model.getModelQsOnNodes(0), 1.0e9f);
  EXPECT_FLOAT_EQ(model.getModelQpOnElement(0), 1.0e9f);
  EXPECT_FLOAT_EQ(model.getModelQsOnElement(0), 1.0e9f);
}

TYPED_TEST(ModelStructTest, SetQualityFactorsFillsAllElements) {
  auto& model = *this->model_;
  using FloatType = typename TestFixture::FloatType;
  FloatType qp = 100.0f, qs = 50.0f;
  model.setQualityFactors(qp, qs);
  for (int e = 0; e < 10; ++e) {
    EXPECT_FLOAT_EQ(model.getModelQpOnElement(e), qp);
    EXPECT_FLOAT_EQ(model.getModelQsOnElement(e), qs);
  }
}

// ============================================================================
// Anisotropy defaults (Delta, Epsilon, Gamma, Theta, Phi)
// ============================================================================

TYPED_TEST(ModelStructTest, AnisotropyDefaultsZeroOnNodes) {
  auto& model = *this->model_;
  EXPECT_FLOAT_EQ(model.getModelDeltaOnNodes(0), 0.0f);
  EXPECT_FLOAT_EQ(model.getModelEpsilonOnNodes(0), 0.0f);
  EXPECT_FLOAT_EQ(model.getModelGammaOnNodes(0), 0.0f);
  EXPECT_EQ(model.getModelThetaOnNodes(0), 0);
  EXPECT_EQ(model.getModelPhiOnNodes(0), 0);
}

TYPED_TEST(ModelStructTest, AnisotropyDefaultsZeroOnElements) {
  auto& model = *this->model_;
  EXPECT_FLOAT_EQ(model.getModelDeltaOnElement(0), 0.0f);
  EXPECT_FLOAT_EQ(model.getModelEpsilonOnElement(0), 0.0f);
  EXPECT_FLOAT_EQ(model.getModelGammaOnElement(0), 0.0f);
  EXPECT_EQ(model.getModelThetaOnElement(0), 0);
  EXPECT_EQ(model.getModelPhiOnElement(0), 0);
}

// ============================================================================
// Populated material arrays (per-node and per-element)
// ============================================================================

TYPED_TEST(ModelStructTest, PerElementMaterialArraysUsedWhenPopulated) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;

  ModelStructDataType data;
  data.ex_ = 2;
  data.ey_ = 2;
  data.ez_ = 2;
  data.dx_ = 20.0;
  data.dy_ = 20.0;
  data.dz_ = 20.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = true;

  int n_elem = 8;
  data.model_vp_element_ = allocateVector<vectorReal>(n_elem, "vp_e");
  data.model_vs_element_ = allocateVector<vectorReal>(n_elem, "vs_e");
  data.model_rho_element_ = allocateVector<vectorReal>(n_elem, "rho_e");
  for (int e = 0; e < n_elem; ++e) {
    data.model_vp_element_(e) = 3000.0f;
    data.model_vs_element_(e) = 1500.0f;
    data.model_rho_element_(e) = 2000.0f;
  }

  ModelStructType model(data);
  EXPECT_FLOAT_EQ(model.getModelVpOnElement(0), 3000.0f);
  EXPECT_FLOAT_EQ(model.getModelVsOnElement(0), 1500.0f);
  EXPECT_FLOAT_EQ(model.getModelRhoOnElement(0), 2000.0f);
}

TYPED_TEST(ModelStructTest, PerNodeMaterialArraysUsedWhenPopulated) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;
  using FloatType = typename TestFixture::FloatType;

  ModelStructDataType data;
  data.ex_ = 2;
  data.ey_ = 2;
  data.ez_ = 2;
  data.dx_ = 20.0;
  data.dy_ = 20.0;
  data.dz_ = 20.0;
  data.isModelOnNodes_ = true;
  data.isElastic_ = true;

  constexpr int Order = TestFixture::Order;
  int n_nodes = (Order * 2 + 1) * (Order * 2 + 1) * (Order * 2 + 1);
  data.model_vp_node_ = allocateVector<vectorReal>(n_nodes, "vp_n");
  data.model_vs_node_ = allocateVector<vectorReal>(n_nodes, "vs_n");
  data.model_rho_node_ = allocateVector<vectorReal>(n_nodes, "rho_n");
  for (int n = 0; n < n_nodes; ++n) {
    data.model_vp_node_(n) = 2500.0f;
    data.model_vs_node_(n) = 1200.0f;
    data.model_rho_node_(n) = 1800.0f;
  }

  ModelStructType model(data);
  EXPECT_FLOAT_EQ(model.getModelVpOnNodes(0), 2500.0f);
  EXPECT_FLOAT_EQ(model.getModelVsOnNodes(0), 1200.0f);
  EXPECT_FLOAT_EQ(model.getModelRhoOnNodes(0), 1800.0f);
}

TYPED_TEST(ModelStructTest, SetModelNodePropsUpdatesValues) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;

  ModelStructDataType data;
  data.ex_ = 2;
  data.ey_ = 2;
  data.ez_ = 2;
  data.dx_ = 20.0;
  data.dy_ = 20.0;
  data.dz_ = 20.0;
  data.isModelOnNodes_ = true;
  data.isElastic_ = true;

  constexpr int Order = TestFixture::Order;
  int n_nodes = (Order * 2 + 1) * (Order * 2 + 1) * (Order * 2 + 1);
  data.model_vp_node_ = allocateVector<vectorReal>(n_nodes, "vp_n");
  data.model_vs_node_ = allocateVector<vectorReal>(n_nodes, "vs_n");
  data.model_rho_node_ = allocateVector<vectorReal>(n_nodes, "rho_n");
  for (int n = 0; n < n_nodes; ++n) {
    data.model_vp_node_(n) = 1500.0f;
    data.model_vs_node_(n) = 755.0f;
    data.model_rho_node_(n) = 1.0f;
  }

  ModelStructType model(data);
  model.setModelNodeProps(0, 3000.0f, 1500.0f, 2000.0f);
  EXPECT_FLOAT_EQ(model.getModelVpOnNodes(0), 3000.0f);
  EXPECT_FLOAT_EQ(model.getModelVsOnNodes(0), 1500.0f);
  EXPECT_FLOAT_EQ(model.getModelRhoOnNodes(0), 2000.0f);
  // Other nodes unchanged
  EXPECT_FLOAT_EQ(model.getModelVpOnNodes(1), 1500.0f);
}

// ============================================================================
// getMinSpacing — all orders
// ============================================================================

TYPED_TEST(ModelStructTest, MinSpacingForAllOrders) {
  auto& model = *this->model_;
  constexpr int Order = TestFixture::Order;
  using FloatType = typename TestFixture::FloatType;

  auto min_spacing = model.getMinSpacing();
  EXPECT_GT(min_spacing, 0.0f);

  // h = lx/ex = 100/10 = 10 for all orders
  FloatType h = 10.0f;
  // NOTE: braces are required around each EXPECT here — nvc++ mis-parses an
  // un-braced gtest macro used as the body of `if constexpr`.
  if constexpr (Order == 1) {
    EXPECT_NEAR(min_spacing, h, 1e-4f);
  }
  if constexpr (Order == 2) {
    EXPECT_NEAR(min_spacing, h * 0.5000000000f, 1e-4f);
  }
  if constexpr (Order == 3) {
    EXPECT_NEAR(min_spacing, h * 0.2763932023f, 1e-4f);
  }
#if MAX_SOLVER_ORDER >= 4
  if constexpr (Order == 4) {
    EXPECT_NEAR(min_spacing, h * 0.1726731646f, 1e-4f);
  }
#endif
#if MAX_SOLVER_ORDER >= 5
  if constexpr (Order == 5) {
    EXPECT_NEAR(min_spacing, h * 0.1174723380f, 1e-4f);
  }
#endif
}

// ============================================================================
// initElasticityTensors — kVTI returns early
// ============================================================================

TYPED_TEST(ModelStructTest, InitElasticityTensorsVTIReturnsEarly) {
  using ModelStructDataType = typename TestFixture::ModelStructDataType;
  using ModelStructType = typename TestFixture::ModelStructType;

  ModelStructDataType data;
  data.ex_ = 2;
  data.ey_ = 2;
  data.ez_ = 2;
  data.dx_ = 20.0;
  data.dy_ = 20.0;
  data.dz_ = 20.0;
  data.isElastic_ = true;
  data.isModelOnNodes_ = false;

  ModelStructType model(data);
  // kVTI returns early without allocating the C-tensor array
  EXPECT_NO_THROW(model.initElasticityTensors(model::kVTI));
}

// ============================================================================
// getMaxSpeed
// ============================================================================

TYPED_TEST(ModelStructTest, MaxSpeedReturns1500) {
  auto& model = *this->model_;
  EXPECT_FLOAT_EQ(model.getMaxSpeed(), 1500.0f);
}

}  // namespace
}  // namespace model

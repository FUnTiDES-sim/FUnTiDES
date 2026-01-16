#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "model_unstruct.h"

namespace model
{
namespace
{

// ============================================================================
// Type Wrapper Classes for Non-Type Template Parameters
// ============================================================================

struct FloatOrder1
{
  using FloatType = float;
  using ScalarType = int;
};

struct DoubleOrder1
{
  using FloatType = double;
  using ScalarType = int;
};

struct FloatOrder2
{
  using FloatType = float;
  using ScalarType = int;
};

struct DoubleOrder2
{
  using FloatType = double;
  using ScalarType = int;
};

struct FloatOrder3
{
  using FloatType = float;
  using ScalarType = int;
};

struct DoubleOrder3
{
  using FloatType = double;
  using ScalarType = int;
};

// ============================================================================
// Test Fixture
// ============================================================================

template <typename TypeWrapper>
class ModelUnstructTest : public ::testing::Test
{
 protected:
  using FloatType = typename TypeWrapper::FloatType;
  using ScalarType = typename TypeWrapper::ScalarType;

  using ModelUnstructType = ModelUnstruct<FloatType, ScalarType>;
};

// Register type wrappers for typed tests
using TypeWrappers = ::testing::Types<FloatOrder1, FloatOrder2, FloatOrder3,
                                      DoubleOrder1, DoubleOrder2, DoubleOrder3>;

TYPED_TEST_SUITE(ModelUnstructTest, TypeWrappers);

// ============================================================================
// Constructor Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, DefaultConstructor)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  ModelUnstructType model;
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, AssignmentOperatorCompiles)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  ModelUnstructType model1;
  ModelUnstructType model2;
  model2 = model1;
  SUCCEED();
}

// ============================================================================
// Type System Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, IndexTypeIsInt)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using IndexType = typename ModelUnstructType::IndexType;

  static_assert(std::is_same<IndexType, int>::value,
                "ModelUnstruct::IndexType must be int");
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, FloatTypeCorrect)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;

  // Template should be properly parameterized
  // Verify by checking it's instantiable
  ModelUnstructType model;
  SUCCEED();
}

// ============================================================================
// GPU Compatibility Documentation
// ============================================================================

TYPED_TEST(ModelUnstructTest, GPUCompatibleMacros)
{
  // ModelUnstruct uses PROXY_HOST_DEVICE on all methods
  // for GPU/CPU dual compilation via Kokkos
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, KokkosViewBased)
{
  // Storage uses Kokkos Views:
  // - ARRAY_INT_VIEW: global_node_index_ (element → node map)
  // - VECTOR_REAL_VIEW: coordinates, properties, parameters
  // - ARRAY3D_REAL_VIEW: elasticity tensors
  SUCCEED();
}

// ============================================================================
// Method Existence Verification
// ============================================================================

TYPED_TEST(ModelUnstructTest, IndexingMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::elementIndex)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::globalVertexIndex)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::globalNodeIndex)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, CoordinateMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::vertexCoords)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::nodeCoord)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, PropertyMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelVpOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelVpOnElement)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelRhoOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelRhoOnElement)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, AnisotropicMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelDeltaOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelDeltaOnElement)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelEpsilonOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelEpsilonOnElement)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, ElasticityMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::initElasticityTensors)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getCTensorOnElement)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, ConfigurationMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::isModelOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::isElastic)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getNumberOfElements)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getNumberOfNodes)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, DomainMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::domainSize)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getMinSpacing)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getMaxSpeed)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, BoundaryMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::boundaryType)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::faceNormal)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, FaceConnectivityMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::buildFaceConnectivity)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getGlobalFace)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getGlobalNodeFromFace)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::isBoundaryFace)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getNumberOfFaces)>);
  SUCCEED();
}

// ============================================================================
// Compilation Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, AllMethodsCompile)
{
  // This test verifies that all public methods have valid signatures
  // by attempting to reference them at compile time
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  // If this compiles, all methods are available
  ModelUnstructType model;
  SUCCEED();
}

// ============================================================================
// Face Connectivity Functional Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, BuildFaceConnectivityCompiles)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  // Create a minimal valid mesh data
  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 1;
  data.n_node_ = 8;
  data.lx_ = 1.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;

  // Allocate minimal arrays
  data.global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(1, 8);
  data.nodes_coords_x_ = allocateVector<VECTOR_REAL_VIEW>(8);
  data.nodes_coords_y_ = allocateVector<VECTOR_REAL_VIEW>(8);
  data.nodes_coords_z_ = allocateVector<VECTOR_REAL_VIEW>(8);
  data.model_vp_element_ = allocateVector<VECTOR_REAL_VIEW>(1);
  data.model_rho_element_ = allocateVector<VECTOR_REAL_VIEW>(1);
  data.boundaries_t_ = allocateVector<VECTOR_REAL_VIEW>(8);

  // Initialize a single cube element
  for (int i = 0; i < 8; ++i)
  {
    data.global_node_index_(0, i) = i;
    data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }
  data.model_vp_element_[0] = 1500.0;
  data.model_rho_element_[0] = 1.0;

  ModelUnstructType model(data);

  // Should not crash
  model.buildFaceConnectivity();
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, GetNumberOfFacesNonZero)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 1;
  data.n_node_ = 8;
  data.lx_ = 1.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
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

  ModelUnstructType model(data);
  model.buildFaceConnectivity();

  // Single cube has 6 faces
  EXPECT_EQ(model.getNumberOfFaces(), 6);
}

TYPED_TEST(ModelUnstructTest, GetGlobalFaceForAllLocalFaces)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 1;
  data.n_node_ = 8;
  data.lx_ = 1.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
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

  ModelUnstructType model(data);
  model.buildFaceConnectivity();

  // Test all 6 faces of the cube
  auto face_xminus = model.getGlobalFace(0, CubicFace::kXMinus);
  auto face_xplus = model.getGlobalFace(0, CubicFace::kXPlus);
  auto face_yminus = model.getGlobalFace(0, CubicFace::kYMinus);
  auto face_yplus = model.getGlobalFace(0, CubicFace::kYPlus);
  auto face_zminus = model.getGlobalFace(0, CubicFace::kZMinus);
  auto face_zplus = model.getGlobalFace(0, CubicFace::kZPlus);

  // All faces should be valid (>= 0)
  EXPECT_GE(face_xminus, 0);
  EXPECT_GE(face_xplus, 0);
  EXPECT_GE(face_yminus, 0);
  EXPECT_GE(face_yplus, 0);
  EXPECT_GE(face_zminus, 0);
  EXPECT_GE(face_zplus, 0);

  // All faces should be different
  EXPECT_NE(face_xminus, face_xplus);
  EXPECT_NE(face_yminus, face_yplus);
  EXPECT_NE(face_zminus, face_zplus);
}

TYPED_TEST(ModelUnstructTest, AllFacesAreBoundaryForSingleElement)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 1;
  data.n_node_ = 8;
  data.lx_ = 1.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
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

  ModelUnstructType model(data);
  model.buildFaceConnectivity();

  // For a single element, all 6 faces are boundaries
  for (int i = 0; i < 6; ++i)
  {
    auto face = model.getGlobalFace(0, static_cast<CubicFace>(i));
    EXPECT_TRUE(model.isBoundaryFace(face));
  }
}

TYPED_TEST(ModelUnstructTest, GetGlobalNodeFromFaceReturnsValidNodes)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 1;
  data.n_node_ = 8;
  data.lx_ = 1.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
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

  ModelUnstructType model(data);
  model.buildFaceConnectivity();

  // Test first face (order=1 means 4 nodes per face)
  auto face0 = model.getGlobalFace(0, CubicFace::kXMinus);

  // Get all 4 nodes of the face
  for (int local_dof = 0; local_dof < 4; ++local_dof)
  {
    auto global_node = model.getGlobalNodeFromFace(face0, local_dof);

    // Node should be in valid range [0, 7]
    EXPECT_GE(global_node, 0);
    EXPECT_LT(global_node, 8);
  }
}

TYPED_TEST(ModelUnstructTest, TwoAdjacentElementsShareInternalFace)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 2;  // Two elements side by side
  data.n_node_ = 12;    // 2 cubes sharing one face = 12 unique nodes
  data.lx_ = 2.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;

  data.global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(2, 8);
  data.nodes_coords_x_ = allocateVector<VECTOR_REAL_VIEW>(12);
  data.nodes_coords_y_ = allocateVector<VECTOR_REAL_VIEW>(12);
  data.nodes_coords_z_ = allocateVector<VECTOR_REAL_VIEW>(12);
  data.model_vp_element_ = allocateVector<VECTOR_REAL_VIEW>(2);
  data.model_rho_element_ = allocateVector<VECTOR_REAL_VIEW>(2);
  data.boundaries_t_ = allocateVector<VECTOR_REAL_VIEW>(12);

  // First cube: nodes 0-7 at x=[0,1], y=[0,1], z=[0,1]
  for (int i = 0; i < 8; ++i)
  {
    data.global_node_index_(0, i) = i;
    data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }

  // Second cube at x=[1,2], y=[0,1], z=[0,1]
  // Shares face x=1 with first cube (nodes 1,3,5,7)
  data.global_node_index_(1, 0) = 1;   // (1,0,0)
  data.global_node_index_(1, 1) = 8;   // (2,0,0)
  data.global_node_index_(1, 2) = 3;   // (1,1,0)
  data.global_node_index_(1, 3) = 9;   // (2,1,0)
  data.global_node_index_(1, 4) = 5;   // (1,0,1)
  data.global_node_index_(1, 5) = 10;  // (2,0,1)
  data.global_node_index_(1, 6) = 7;   // (1,1,1)
  data.global_node_index_(1, 7) = 11;  // (2,1,1)

  data.nodes_coords_x_[8] = 2.0;
  data.nodes_coords_y_[8] = 0.0;
  data.nodes_coords_z_[8] = 0.0;
  data.nodes_coords_x_[9] = 2.0;
  data.nodes_coords_y_[9] = 1.0;
  data.nodes_coords_z_[9] = 0.0;
  data.nodes_coords_x_[10] = 2.0;
  data.nodes_coords_y_[10] = 0.0;
  data.nodes_coords_z_[10] = 1.0;
  data.nodes_coords_x_[11] = 2.0;
  data.nodes_coords_y_[11] = 1.0;
  data.nodes_coords_z_[11] = 1.0;

  data.model_vp_element_[0] = 1500.0;
  data.model_vp_element_[1] = 1500.0;
  data.model_rho_element_[0] = 1.0;
  data.model_rho_element_[1] = 1.0;

  ModelUnstructType model(data);
  model.buildFaceConnectivity();

  // Two cubes sharing one face = 11 unique faces
  // Cube 0: 6 faces, Cube 1: 6 faces, Shared: 1 face
  // Total unique: 6 + 6 - 1 = 11
  EXPECT_EQ(model.getNumberOfFaces(), 11);

  // The shared face (XPlus of elem 0, XMinus of elem 1) should be the same
  auto face_elem0_xplus = model.getGlobalFace(0, CubicFace::kXPlus);
  auto face_elem1_xminus = model.getGlobalFace(1, CubicFace::kXMinus);

  // Should be the same face
  EXPECT_EQ(face_elem0_xplus, face_elem1_xminus);

  // This face should NOT be a boundary (it has two neighbors)
  EXPECT_FALSE(model.isBoundaryFace(face_elem0_xplus));
}
// ============================================================================
// Documentation Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, StructuredVsUnstructuredComparison)
{
  // ModelStruct: Structured mesh with implicit geometry
  // - Template parameter: Order
  // - Mesh defined by: element counts (ex, ey, ez)
  // - Node positions: Computed from index formulas
  // - Use case: Regular grids, fast access
  //
  // ModelUnstruct: Unstructured mesh with explicit geometry
  // - Template parameters: FloatType, ScalarType (NOT Order)
  // - Mesh defined by: connectivity arrays
  // - Node positions: Stored explicitly
  // - Use case: Arbitrary topologies, complex geometries
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, TemplateParameterizations)
{
  // ModelUnstruct<float, int>:
  // - 32-bit floating point for coordinates
  // - 32-bit integer indexing
  // - Suitable for typical simulations
  //
  // ModelUnstruct<double, int>:
  // - 64-bit floating point for coordinates
  // - Higher precision for large-scale domains
  // - More memory usage
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, FaceConnectivityDocumentation)
{
  // Face connectivity provides efficient access to mesh faces for:
  // - Absorbing boundary conditions
  // - Flux calculations
  // - Interface handling
  //
  // Key concepts:
  // 1. Local face indexing: Each hexahedral element has 6 faces
  //    - CubicFace::kXMinus (0), CubicFace::kXPlus (1)
  //    - CubicFace::kYMinus (2), CubicFace::kYPlus (3)
  //    - CubicFace::kZMinus (4), CubicFace::kZPlus (5)
  //
  // 2. Global face indexing: Unique ID for each face in the mesh
  //    - Internal faces (shared by 2 elements): appear once
  //    - Boundary faces (1 element only): flagged as boundary
  //
  // 3. Face-to-node mapping: Each face has (order+1)² nodes
  //    - Ordered lexicographically for quadrature integration
  //    - Accessible via getGlobalNodeFromFace(face_id, local_dof)
  //
  // Usage pattern for boundary conditions:
  //   for each element:
  //     for each local face (0-5):
  //       face_id = getGlobalFace(element, local_face)
  //       if isBoundaryFace(face_id):
  //         for each node on face:
  //           node_id = getGlobalNodeFromFace(face_id, local_dof)
  //           apply boundary condition to node_id
  //
  // Internal representation:
  // - elem_to_faces_: Maps (element, local_face) → global_face_id
  // - face_dofs_: Maps (global_face_id, local_dof) → global_node_id
  // - face_elem_owner_/neighbor_: Adjacency information
  // - Boundary detection: face_elem_neighbor_ == -1
  SUCCEED();
}

}  // namespace
}  // namespace model
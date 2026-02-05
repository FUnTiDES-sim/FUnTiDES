#include <gtest/gtest.h>

#include "face_connectivity_unstruct.h"
#include "model_unstruct.h"

namespace model
{
namespace
{

// ============================================================================
// Helper: Create simple test mesh
// ============================================================================

template <typename FloatType, typename ScalarType>
ModelUnstruct<FloatType, ScalarType> createSingleCubeMesh()
{
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

  return ModelUnstruct<FloatType, ScalarType>(data);
}

template <typename FloatType, typename ScalarType>
ModelUnstruct<FloatType, ScalarType> createTwoAdjacentCubes()
{
  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 2;
  data.n_node_ = 12;
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

  // First cube: nodes 0-7 at x=[0,1]
  for (int i = 0; i < 8; ++i)
  {
    data.global_node_index_(0, i) = i;
    data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }

  // Second cube at x=[1,2], shares face x=1
  data.global_node_index_(1, 0) = 1;
  data.global_node_index_(1, 1) = 8;
  data.global_node_index_(1, 2) = 3;
  data.global_node_index_(1, 3) = 9;
  data.global_node_index_(1, 4) = 5;
  data.global_node_index_(1, 5) = 10;
  data.global_node_index_(1, 6) = 7;
  data.global_node_index_(1, 7) = 11;

  for (int i = 8; i < 12; ++i)
  {
    data.nodes_coords_x_[i] = 2.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }

  data.model_vp_element_[0] = 1500.0;
  data.model_vp_element_[1] = 1500.0;
  data.model_rho_element_[0] = 1.0;
  data.model_rho_element_[1] = 1.0;

  return ModelUnstruct<FloatType, ScalarType>(data);
}

// ============================================================================
// Test Fixture
// ============================================================================

template <typename T>
class FaceConnectivityTest : public ::testing::Test
{
 protected:
  using FloatType = typename T::FloatType;
  using ScalarType = typename T::ScalarType;
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
TYPED_TEST_SUITE(FaceConnectivityTest, TestTypes);

// ============================================================================
// Build Tests
// ============================================================================

TYPED_TEST(FaceConnectivityTest, BuildFromSingleCube)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createSingleCubeMesh<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  EXPECT_EQ(fc.getNumberOfFaces(), 6);
  EXPECT_EQ(fc.getDofsPerFace(), 4);  // order=1 → (1+1)² = 4
}

TYPED_TEST(FaceConnectivityTest, BuildFromTwoAdjacentCubes)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createTwoAdjacentCubes<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  // 2 cubes × 6 faces - 1 shared = 11 unique faces
  EXPECT_EQ(fc.getNumberOfFaces(), 11);
}

// ============================================================================
// Face Mapping Tests
// ============================================================================

TYPED_TEST(FaceConnectivityTest, GetGlobalFaceReturnsValidIds)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createSingleCubeMesh<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  // All 6 local faces should map to valid global faces
  for (int lf = 0; lf < 6; ++lf)
  {
    auto face_id = fc.getGlobalFace(0, static_cast<CubicFace>(lf));
    EXPECT_GE(face_id, 0);
    EXPECT_LT(face_id, 6);
  }
}

TYPED_TEST(FaceConnectivityTest, GetGlobalFaceUniqueness)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createSingleCubeMesh<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  std::set<ScalarType> face_ids;
  for (int lf = 0; lf < 6; ++lf)
  {
    auto face_id = fc.getGlobalFace(0, static_cast<CubicFace>(lf));
    face_ids.insert(face_id);
  }

  // All 6 faces should be different
  EXPECT_EQ(face_ids.size(), 6);
}

// ============================================================================
// Boundary Detection Tests
// ============================================================================

TYPED_TEST(FaceConnectivityTest, SingleCubeAllFacesAreBoundary)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createSingleCubeMesh<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  // Single element → all faces are boundaries
  for (ScalarType f = 0; f < fc.getNumberOfFaces(); ++f)
  {
    EXPECT_TRUE(fc.isBoundaryFace(f));
  }
}

TYPED_TEST(FaceConnectivityTest, TwoCubesSharedFaceNotBoundary)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createTwoAdjacentCubes<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  auto face_elem0_xplus = fc.getGlobalFace(0, CubicFace::kXPlus);
  auto face_elem1_xminus = fc.getGlobalFace(1, CubicFace::kXMinus);

  // Should be the same face (shared)
  EXPECT_EQ(face_elem0_xplus, face_elem1_xminus);

  // Should NOT be a boundary
  EXPECT_FALSE(fc.isBoundaryFace(face_elem0_xplus));
}

TYPED_TEST(FaceConnectivityTest, TwoCubesCountBoundaryFaces)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createTwoAdjacentCubes<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  int boundary_count = 0;
  for (ScalarType f = 0; f < fc.getNumberOfFaces(); ++f)
  {
    if (fc.isBoundaryFace(f)) boundary_count++;
  }

  // 2 cubes: 10 external faces, 1 internal
  EXPECT_EQ(boundary_count, 10);
}

// ============================================================================
// Node-on-Face Tests
// ============================================================================

TYPED_TEST(FaceConnectivityTest, GetGlobalNodeFromFaceReturnsValidNodes)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createSingleCubeMesh<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  auto face0 = fc.getGlobalFace(0, CubicFace::kXMinus);

  // Order=1 → 4 nodes per face
  for (int dof = 0; dof < 4; ++dof)
  {
    auto node = fc.getGlobalNodeFromFace(face0, dof);
    EXPECT_GE(node, 0);
    EXPECT_LT(node, 8);
  }
}

TYPED_TEST(FaceConnectivityTest, FaceNodesAreUnique)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createSingleCubeMesh<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  auto face0 = fc.getGlobalFace(0, CubicFace::kXMinus);

  std::set<ScalarType> nodes;
  for (int dof = 0; dof < 4; ++dof)
  {
    nodes.insert(fc.getGlobalNodeFromFace(face0, dof));
  }

  // All 4 nodes should be different
  EXPECT_EQ(nodes.size(), 4);
}

// ============================================================================
// Adjacency Tests
// ============================================================================

TYPED_TEST(FaceConnectivityTest, InternalFaceHasTwoOwners)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createTwoAdjacentCubes<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  auto shared_face = fc.getGlobalFace(0, CubicFace::kXPlus);

  auto owner = fc.elemOwner(shared_face);
  auto neighbor = fc.elemNeighbor(shared_face);

  EXPECT_GE(owner, 0);
  EXPECT_GE(neighbor, 0);
  EXPECT_NE(owner, neighbor);
}

TYPED_TEST(FaceConnectivityTest, BoundaryFaceHasNoNeighbor)
{
  using FloatType = typename TestFixture::FloatType;
  using ScalarType = typename TestFixture::ScalarType;

  auto mesh = createSingleCubeMesh<FloatType, ScalarType>();
  auto fc = FaceConnectivity<FloatType, ScalarType>::build(mesh);

  auto boundary_face = fc.getGlobalFace(0, CubicFace::kXMinus);

  EXPECT_GE(fc.elemOwner(boundary_face), 0);
  EXPECT_EQ(fc.elemNeighbor(boundary_face), -1);  // No neighbor
}

}  // namespace
}  // namespace model
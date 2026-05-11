#include <gtest/gtest.h>

#include "face_connectivity_unstruct.h"
#include "model_unstruct.h"
#include "test_face_connectivity_helpers.h"

namespace model {
namespace {

template <typename FloatType, typename ScalarType>
ModelUnstruct<FloatType, ScalarType> createSingleCube() {
  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 1;
  data.n_node_ = 8;
  data.lx_ = 1.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;
  data.global_node_index_ = allocateArray2D<arrayInt>(1, 8);
  data.nodes_coords_x_ = allocateVector<vectorReal>(8);
  data.nodes_coords_y_ = allocateVector<vectorReal>(8);
  data.nodes_coords_z_ = allocateVector<vectorReal>(8);
  data.model_vp_element_ = allocateVector<vectorReal>(1);
  data.model_rho_element_ = allocateVector<vectorReal>(1);
  data.boundaries_t_ = allocateVector<vectorInt>(8);
  for (int i = 0; i < 8; ++i) {
    data.global_node_index_(0, i) = i;
    data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }
  data.model_vp_element_[0] = 1500.0;
  data.model_rho_element_[0] = 1.0;
  return ModelUnstruct<FloatType, ScalarType>(data);
}

/// Two Q1 cubes sharing an X-face where element 1 has its local j and k axes
/// swapped relative to element 0. This produces a non-trivial face permutation
/// perm = {0, 2, 1, 3} on the shared face (swap of the two tangential DOF
/// directions).
template <typename FloatType, typename ScalarType>
ModelUnstruct<FloatType, ScalarType> createTwoAdjacentCubesRotated() {
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
  data.boundaries_t_ = allocateVector<VECTOR_INT_VIEW>(12);
  // Elem 0: standard cube [0,1]^3, nodes 0-7
  for (int i = 0; i < 8; ++i) {
    data.global_node_index_(0, i) = i;
    data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }
  // Elem 1: local j-axis points in global z-direction, local k-axis in global
  // y-direction (90-degree rotation of the tangential axes on the shared face).
  // local (i,j,k) → global node mapping:
  //   (0,0,0)→1  (1,0,0)→8   (0,1,0)→5  (1,1,0)→10
  //   (0,0,1)→3  (1,0,1)→9   (0,1,1)→7  (1,1,1)→11
  data.global_node_index_(1, 0) = 1;
  data.global_node_index_(1, 1) = 8;
  data.global_node_index_(1, 2) = 5;
  data.global_node_index_(1, 3) = 10;
  data.global_node_index_(1, 4) = 3;
  data.global_node_index_(1, 5) = 9;
  data.global_node_index_(1, 6) = 7;
  data.global_node_index_(1, 7) = 11;
  // New nodes at x=2
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
  data.model_vp_element_[0] = data.model_vp_element_[1] = 1500.0;
  data.model_rho_element_[0] = data.model_rho_element_[1] = 1.0;
  return ModelUnstruct<FloatType, ScalarType>(data);
}

template <typename FloatType, typename ScalarType>
ModelUnstruct<FloatType, ScalarType> createTwoAdjacentCubes() {
  ModelUnstructData<FloatType, ScalarType> data;
  data.order_ = 1;
  data.n_element_ = 2;
  data.n_node_ = 12;
  data.lx_ = 2.0;
  data.ly_ = 1.0;
  data.lz_ = 1.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;
  data.global_node_index_ = allocateArray2D<arrayInt>(2, 8);
  data.nodes_coords_x_ = allocateVector<vectorReal>(12);
  data.nodes_coords_y_ = allocateVector<vectorReal>(12);
  data.nodes_coords_z_ = allocateVector<vectorReal>(12);
  data.model_vp_element_ = allocateVector<vectorReal>(2);
  data.model_rho_element_ = allocateVector<vectorReal>(2);
  data.boundaries_t_ = allocateVector<vectorInt>(12);
  for (int i = 0; i < 8; ++i) {
    data.global_node_index_(0, i) = i;
    data.nodes_coords_x_[i] = (i & 1) ? 1.0 : 0.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }
  data.global_node_index_(1, 0) = 1;
  data.global_node_index_(1, 1) = 8;
  data.global_node_index_(1, 2) = 3;
  data.global_node_index_(1, 3) = 9;
  data.global_node_index_(1, 4) = 5;
  data.global_node_index_(1, 5) = 10;
  data.global_node_index_(1, 6) = 7;
  data.global_node_index_(1, 7) = 11;
  for (int i = 8; i < 12; ++i) {
    data.nodes_coords_x_[i] = 2.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }
  data.model_vp_element_[0] = data.model_vp_element_[1] = 1500.0;
  data.model_rho_element_[0] = data.model_rho_element_[1] = 1.0;
  return ModelUnstruct<FloatType, ScalarType>(data);
}

template <typename T>
class FaceConnectivityUnstructTest : public ::testing::Test {
 protected:
  using FloatType = typename T::FloatType;
  using ScalarType = typename T::ScalarType;
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
TYPED_TEST_SUITE(FaceConnectivityUnstructTest, TestTypes);

TYPED_TEST(FaceConnectivityUnstructTest, BuildFromSingleCube) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createSingleCube<F, S>());
  testSingleCubeFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, BuildFromTwoAdjacentCubes) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createTwoAdjacentCubes<F, S>());
  testTwoAdjacentCubesFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, GetGlobalFaceReturnsValidIds) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createSingleCube<F, S>());
  testGlobalFaceIds<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, GetGlobalFaceUniqueness) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createSingleCube<F, S>());
  testGlobalFaceUniqueness<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, SingleCubeAllFacesAreBoundary) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createSingleCube<F, S>());
  testSingleCubeAllBoundary<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, TwoCubesSharedFaceNotBoundary) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createTwoAdjacentCubes<F, S>());
  testSharedFaceNotBoundary<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, TwoCubesCountBoundaryFaces) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createTwoAdjacentCubes<F, S>());
  testBoundaryFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, FaceNodesValid) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createSingleCube<F, S>());
  testFaceNodes<F, S>(fc, 8);
}

TYPED_TEST(FaceConnectivityUnstructTest, FaceNodesUnique) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createSingleCube<F, S>());
  testFaceNodesUnique<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, InternalFaceHasTwoOwners) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createTwoAdjacentCubes<F, S>());
  testInternalFaceOwners<F, S>(fc);
}

TYPED_TEST(FaceConnectivityUnstructTest, BoundaryFaceHasNoNeighbor) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(createSingleCube<F, S>());
  testBoundaryFaceNoNeighbor<F, S>(fc);
}

// ============================================================================
// faceLocalToElemLocal — free function, order=1
// ============================================================================

TEST(FaceLocalToElemLocalTest, AllFacesOrder1) {
  // order=1: n=2, ndofs_per_face=4
  // face_dof_2d layout: u = dof%2, v = dof/2
  // kXMinus: 0  + u*2 + v*4
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXMinus, 0, 1), 0);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXMinus, 1, 1), 2);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXMinus, 2, 1), 4);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXMinus, 3, 1), 6);
  // kXPlus: 1 + u*2 + v*4
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXPlus, 0, 1), 1);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXPlus, 1, 1), 3);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXPlus, 2, 1), 5);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kXPlus, 3, 1), 7);
  // kYMinus: u + v*4
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYMinus, 0, 1), 0);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYMinus, 1, 1), 1);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYMinus, 2, 1), 4);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYMinus, 3, 1), 5);
  // kYPlus: u + 2 + v*4
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYPlus, 0, 1), 2);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYPlus, 1, 1), 3);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYPlus, 2, 1), 6);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kYPlus, 3, 1), 7);
  // kZMinus: u + v*2
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZMinus, 0, 1), 0);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZMinus, 1, 1), 1);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZMinus, 2, 1), 2);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZMinus, 3, 1), 3);
  // kZPlus: u + v*2 + 4
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZPlus, 0, 1), 4);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZPlus, 1, 1), 5);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZPlus, 2, 1), 6);
  EXPECT_EQ(faceLocalToElemLocal(CubicFace::kZPlus, 3, 1), 7);
}

TEST(FaceLocalToElemLocalTest, CornersSpanAllElemDofs) {
  // Every element DOF must be reachable from exactly one (face, face_dof) pair
  // among the 6 faces — otherwise the mapping has gaps or collisions.
  std::set<int> reached;
  for (int lf = 0; lf < 6; ++lf)
    for (int d = 0; d < 4; ++d) reached.insert(faceLocalToElemLocal(static_cast<CubicFace>(lf), d, 1));
  EXPECT_EQ(static_cast<int>(reached.size()), 8);
}

// ============================================================================
// getNeighborFaceDof — aligned and rotated meshes
// ============================================================================

TYPED_TEST(FaceConnectivityUnstructTest, NeighborFaceDofAlignedIsIdentity) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto mesh = createTwoAdjacentCubes<F, S>();
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(mesh);

  auto shared_face = fc.getGlobalFace(0, CubicFace::kXPlus);
  ASSERT_FALSE(fc.isBoundaryFace(shared_face));

  const int ndofs = fc.getDofsPerFace();
  for (int i = 0; i < ndofs; ++i)
    EXPECT_EQ(fc.getNeighborFaceDof(shared_face, i), i) << "Aligned mesh: expected identity permutation at DOF " << i;
}

TYPED_TEST(FaceConnectivityUnstructTest, NeighborFaceDofRotatedIsCorrect) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto mesh = createTwoAdjacentCubesRotated<F, S>();
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(mesh);

  auto shared_face = fc.getGlobalFace(0, CubicFace::kXPlus);
  ASSERT_FALSE(fc.isBoundaryFace(shared_face));

  // Elem 1 has j↔z and k↔y swapped → perm = {0,2,1,3}
  EXPECT_EQ(fc.getNeighborFaceDof(shared_face, 0), 0);
  EXPECT_EQ(fc.getNeighborFaceDof(shared_face, 1), 2);
  EXPECT_EQ(fc.getNeighborFaceDof(shared_face, 2), 1);
  EXPECT_EQ(fc.getNeighborFaceDof(shared_face, 3), 3);
}

TYPED_TEST(FaceConnectivityUnstructTest, NeighborFaceDofNodesMatch) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto mesh = createTwoAdjacentCubesRotated<F, S>();
  FaceConnectivityUnstruct<F, S> fc;
  fc.build(mesh);

  auto shared_face = fc.getGlobalFace(0, CubicFace::kXPlus);
  ASSERT_FALSE(fc.isBoundaryFace(shared_face));
  testNeighborFaceDofCorrectness<F, S>(fc, mesh, shared_face);
}

}  // namespace
}  // namespace model
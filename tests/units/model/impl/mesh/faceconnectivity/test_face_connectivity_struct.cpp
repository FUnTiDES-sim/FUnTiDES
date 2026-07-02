#include <gtest/gtest.h>

#include "face_connectivity_struct.h"
#include "model_struct.h"
#include "test_face_connectivity_helpers.h"

namespace model {
namespace {

template <typename FloatType, typename ScalarType>
ModelStruct<FloatType, ScalarType, 1> createSingleCube() {
  ModelStructData<FloatType, ScalarType> data;
  data.ex_ = 1;
  data.ey_ = 1;
  data.ez_ = 1;
  data.dx_ = 1.0;
  data.dy_ = 1.0;
  data.dz_ = 1.0;
  data.ox_ = 0.0;
  data.oy_ = 0.0;
  data.oz_ = 0.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;
  return ModelStruct<FloatType, ScalarType, 1>(data);
}

template <typename FloatType, typename ScalarType>
ModelStruct<FloatType, ScalarType, 1> createTwoAdjacentCubes() {
  ModelStructData<FloatType, ScalarType> data;
  data.ex_ = 2;
  data.ey_ = 1;
  data.ez_ = 1;
  data.dx_ = 2.0;
  data.dy_ = 1.0;
  data.dz_ = 1.0;
  data.ox_ = 0.0;
  data.oy_ = 0.0;
  data.oz_ = 0.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;
  return ModelStruct<FloatType, ScalarType, 1>(data);
}

template <typename T>
class FaceConnectivityStructTest : public ::testing::Test {
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
TYPED_TEST_SUITE(FaceConnectivityStructTest, TestTypes);

TYPED_TEST(FaceConnectivityStructTest, BuildFromSingleCube) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testSingleCubeFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, BuildFromTwoAdjacentCubes) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testTwoAdjacentCubesFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, GetGlobalFaceReturnsValidIds) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testGlobalFaceIds<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, GetGlobalFaceUniqueness) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testGlobalFaceUniqueness<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, SingleCubeAllFacesAreBoundary) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testSingleCubeAllBoundary<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, TwoCubesSharedFaceNotBoundary) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testSharedFaceNotBoundary<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, TwoCubesCountBoundaryFaces) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testBoundaryFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, FaceNodesValid) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testFaceNodes<F, S>(fc, 8);
}

TYPED_TEST(FaceConnectivityStructTest, FaceNodesUnique) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testFaceNodesUnique<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, InternalFaceHasTwoOwners) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testInternalFaceOwners<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, BoundaryFaceHasNoNeighbor) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testBoundaryFaceNoNeighbor<F, S>(fc);
}

// ============================================================================
// localFaceOwner / localFaceNeighbor
// ============================================================================

TYPED_TEST(FaceConnectivityStructTest, LocalFaceOwnerXFaces) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // 2x1x1: X-face at i=0 → kXMinus owner, i=1 (interior) → kXMinus owner,
  // i=2 (right boundary) → kXPlus owner
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  // face 0: i=0 → boundary, owner sees kXMinus
  EXPECT_EQ(fc.localFaceOwner(0), static_cast<int>(CubicFace::kXMinus));
  // face 1: i=1 → interior, owner (i<ex) sees kXMinus
  EXPECT_EQ(fc.localFaceOwner(1), static_cast<int>(CubicFace::kXMinus));
  // face 2: i=2 → boundary, owner sees kXPlus
  EXPECT_EQ(fc.localFaceOwner(2), static_cast<int>(CubicFace::kXPlus));
}

TYPED_TEST(FaceConnectivityStructTest, LocalFaceOwnerYFaces) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // 1x2x1: Y-faces start at offset_y = (1+1)*2*1 = 4... wait let me compute:
  // ex=1,ey=2,ez=1: offset_y = (1+1)*2*1 = 4, offset_z = 4 + 1*(2+1)*1 = 7
  // Y-face at j=0 → kYMinus owner, j=2 → kYPlus owner
  FaceConnectivityStruct<F, S> fc(1, 2, 1, 1);
  S offset_y = (1 + 1) * 2 * 1;  // 4
  // j=0: first Y-face → kYMinus
  EXPECT_EQ(fc.localFaceOwner(offset_y), static_cast<int>(CubicFace::kYMinus));
  // j=1 (interior): kYMinus
  EXPECT_EQ(fc.localFaceOwner(offset_y + 1), static_cast<int>(CubicFace::kYMinus));
  // j=2 (last): kYPlus
  EXPECT_EQ(fc.localFaceOwner(offset_y + 2), static_cast<int>(CubicFace::kYPlus));
}

TYPED_TEST(FaceConnectivityStructTest, LocalFaceOwnerZFaces) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // ex=1,ey=1,ez=2
  // offset_y = (1+1)*1*2 = 4, offset_z = 4 + 1*(1+1)*2 = 8
  FaceConnectivityStruct<F, S> fc(1, 1, 2, 1);
  S offset_y = (1 + 1) * 1 * 2;             // 4
  S offset_z = offset_y + 1 * (1 + 1) * 2;  // 8
  // k=0: kZMinus
  EXPECT_EQ(fc.localFaceOwner(offset_z), static_cast<int>(CubicFace::kZMinus));
  // k=1 (interior): kZMinus
  EXPECT_EQ(fc.localFaceOwner(offset_z + 1), static_cast<int>(CubicFace::kZMinus));
  // k=2 (last): kZPlus
  EXPECT_EQ(fc.localFaceOwner(offset_z + 2), static_cast<int>(CubicFace::kZPlus));
}

TYPED_TEST(FaceConnectivityStructTest, LocalFaceNeighborBoundaryReturnsMinusOne) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  for (S f = 0; f < fc.getNumberOfFaces(); ++f) EXPECT_EQ(fc.localFaceNeighbor(f), -1);
}

TYPED_TEST(FaceConnectivityStructTest, LocalFaceNeighborInternalIsOpposite) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // Shared X-face between elem 0 and elem 1 in a 2x1x1 mesh
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  auto shared = fc.getGlobalFace(0, CubicFace::kXPlus);
  int owner_lf = fc.localFaceOwner(shared);
  int neigh_lf = fc.localFaceNeighbor(shared);
  EXPECT_NE(neigh_lf, -1);
  // XOR of opposite face pair must flip bit 0
  EXPECT_EQ(owner_lf ^ 1, neigh_lf);
}

// ============================================================================
// getNeighborFaceDof — identity permutation for structured meshes
// ============================================================================

TYPED_TEST(FaceConnectivityStructTest, NeighborFaceDofIdentityOrder1) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto mesh = createTwoAdjacentCubes<F, S>();
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  auto shared = fc.getGlobalFace(0, CubicFace::kXPlus);
  testNeighborFaceDofCorrectness<F, S>(fc, mesh, shared);
}

TYPED_TEST(FaceConnectivityStructTest, NeighborFaceDofIdentityOrder2) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  ModelStructData<F, S> data;
  data.ex_ = 2;
  data.ey_ = 1;
  data.ez_ = 1;
  data.dx_ = 2.0;
  data.dy_ = 1.0;
  data.dz_ = 1.0;
  data.isModelOnNodes_ = false;
  data.isElastic_ = false;
  auto mesh = ModelStruct<F, S, 2>(data);
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 2);
  auto shared = fc.getGlobalFace(0, CubicFace::kXPlus);
  testNeighborFaceDofCorrectness<F, S>(fc, mesh, shared);
}

// ============================================================================
// getGlobalNodeFromFace — Y-face and Z-face branches
// ============================================================================

TYPED_TEST(FaceConnectivityStructTest, GetGlobalNodeFromFaceYFaceValid) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // 1x2x1 mesh, order 1 → 8 nodes
  FaceConnectivityStruct<F, S> fc(1, 2, 1, 1);
  S offset_y = (1 + 1) * 2 * 1;  // 4
  // First Y-face (j=0)
  for (int dof = 0; dof < 4; ++dof) {
    auto node = fc.getGlobalNodeFromFace(offset_y, dof);
    EXPECT_GE(node, 0);
    EXPECT_LT(node, 2 * 3 * 2);  // nx=2, ny=3, nz=2
  }
}

TYPED_TEST(FaceConnectivityStructTest, GetGlobalNodeFromFaceYFaceUnique) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 2, 1, 1);
  S offset_y = (1 + 1) * 2 * 1;  // 4
  std::set<S> nodes;
  for (int dof = 0; dof < 4; ++dof) nodes.insert(fc.getGlobalNodeFromFace(offset_y, dof));
  EXPECT_EQ(nodes.size(), 4u);
}

TYPED_TEST(FaceConnectivityStructTest, GetGlobalNodeFromFaceZFaceValid) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // 1x1x2 mesh, order 1 → 8 nodes
  FaceConnectivityStruct<F, S> fc(1, 1, 2, 1);
  S offset_y = (1 + 1) * 1 * 2;             // 4
  S offset_z = offset_y + 1 * (1 + 1) * 2;  // 8
  // First Z-face (k=0)
  for (int dof = 0; dof < 4; ++dof) {
    auto node = fc.getGlobalNodeFromFace(offset_z, dof);
    EXPECT_GE(node, 0);
    EXPECT_LT(node, 2 * 2 * 3);  // nx=2, ny=2, nz=3
  }
}

TYPED_TEST(FaceConnectivityStructTest, GetGlobalNodeFromFaceZFaceUnique) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 2, 1);
  S offset_y = (1 + 1) * 1 * 2;
  S offset_z = offset_y + 1 * (1 + 1) * 2;
  std::set<S> nodes;
  for (int dof = 0; dof < 4; ++dof) nodes.insert(fc.getGlobalNodeFromFace(offset_z, dof));
  EXPECT_EQ(nodes.size(), 4u);
}

// ============================================================================
// getDofsPerFace — higher orders
// ============================================================================

TYPED_TEST(FaceConnectivityStructTest, DofsPerFaceOrder2) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 2);
  EXPECT_EQ(fc.getDofsPerFace(), 9);  // (2+1)^2
}

TYPED_TEST(FaceConnectivityStructTest, DofsPerFaceOrder3) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 3);
  EXPECT_EQ(fc.getDofsPerFace(), 16);  // (3+1)^2
}

// ============================================================================
// 3D mesh Y and Z boundary face detection
// ============================================================================

TYPED_TEST(FaceConnectivityStructTest, ThreeDMeshYBoundaryFaces) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // 2x2x2 mesh: Y-faces start at offset_y = (2+1)*2*2 = 12
  // Y boundary faces: j=0 or j=ey=2
  FaceConnectivityStruct<F, S> fc(2, 2, 2, 1);
  S offset_y = (2 + 1) * 2 * 2;             // 12
  S offset_z = offset_y + 2 * (2 + 1) * 2;  // 12+20=32? let me compute: 2*(2+1)*2=12 → offset_z=24
  // j=0: first row of Y-faces, face = offset_y + 0 + 0*2 + k*2*3 for i=0,k=0 → offset_y+0 = 12
  EXPECT_TRUE(fc.isBoundaryFace(offset_y));       // j=0, boundary
  EXPECT_FALSE(fc.isBoundaryFace(offset_y + 2));  // j=1, interior (i=0, j=1, k=0)
  // j=2 (last): i=0, j=2, k=0 → local = 0 + 2*2 + 0 = 4 → face = offset_y + 4
  EXPECT_TRUE(fc.isBoundaryFace(offset_y + 4));  // j=2, boundary
}

TYPED_TEST(FaceConnectivityStructTest, ThreeDMeshZBoundaryFaces) {
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  // 2x2x2 mesh
  FaceConnectivityStruct<F, S> fc(2, 2, 2, 1);
  S offset_y = (2 + 1) * 2 * 2;             // 12
  S offset_z = offset_y + 2 * (2 + 1) * 2;  // 12 + 12 = 24
  // k=0: face = offset_z + 0 = 24 → boundary
  EXPECT_TRUE(fc.isBoundaryFace(offset_z));
  // k=1 (interior): local = i + j*ex + 1*ex*ey = 0 + 0 + 4 = 4 → face = offset_z+4
  EXPECT_FALSE(fc.isBoundaryFace(offset_z + 4));
  // k=2 (last): local = 0 + 0 + 2*4 = 8 → face = offset_z+8
  EXPECT_TRUE(fc.isBoundaryFace(offset_z + 8));
}

}  // namespace
}  // namespace model
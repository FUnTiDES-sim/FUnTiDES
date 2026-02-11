#include <gtest/gtest.h>

#include "face_connectivity_unstruct.h"
#include "model_unstruct.h"
#include "test_face_connectivity_helpers.h"

namespace model
{
namespace
{

template <typename FloatType, typename ScalarType>
ModelUnstruct<FloatType, ScalarType> createSingleCube()
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
  for (int i = 0; i < 8; ++i)
  {
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
  for (int i = 8; i < 12; ++i)
  {
    data.nodes_coords_x_[i] = 2.0;
    data.nodes_coords_y_[i] = (i & 2) ? 1.0 : 0.0;
    data.nodes_coords_z_[i] = (i & 4) ? 1.0 : 0.0;
  }
  data.model_vp_element_[0] = data.model_vp_element_[1] = 1500.0;
  data.model_rho_element_[0] = data.model_rho_element_[1] = 1.0;
  return ModelUnstruct<FloatType, ScalarType>(data);
}

template <typename T>
class FaceConnectivityUnstructTest : public ::testing::Test
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
TYPED_TEST_SUITE(FaceConnectivityUnstructTest, TestTypes);

TYPED_TEST(FaceConnectivityUnstructTest, BuildFromSingleCube)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc = FaceConnectivityUnstruct<F, S>().build(createSingleCube<F, S>());
  testSingleCubeFaceCount<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, BuildFromTwoAdjacentCubes)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc =
      FaceConnectivityUnstruct<F, S>().build(createTwoAdjacentCubes<F, S>());
  testTwoAdjacentCubesFaceCount<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, GetGlobalFaceReturnsValidIds)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc = FaceConnectivityUnstruct<F, S>().build(createSingleCube<F, S>());
  testGlobalFaceIds<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, GetGlobalFaceUniqueness)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc = FaceConnectivityUnstruct<F, S>().build(createSingleCube<F, S>());
  testGlobalFaceUniqueness<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, SingleCubeAllFacesAreBoundary)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc = FaceConnectivityUnstruct<F, S>().build(createSingleCube<F, S>());
  testSingleCubeAllBoundary<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, TwoCubesSharedFaceNotBoundary)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc =
      FaceConnectivityUnstruct<F, S>().build(createTwoAdjacentCubes<F, S>());
  testSharedFaceNotBoundary<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, TwoCubesCountBoundaryFaces)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc =
      FaceConnectivityUnstruct<F, S>().build(createTwoAdjacentCubes<F, S>());
  testBoundaryFaceCount<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, FaceNodesValid)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc = FaceConnectivityUnstruct<F, S>().build(createSingleCube<F, S>());
  testFaceNodes<F, S>(fc, 8);
}
TYPED_TEST(FaceConnectivityUnstructTest, FaceNodesUnique)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc = FaceConnectivityUnstruct<F, S>().build(createSingleCube<F, S>());
  testFaceNodesUnique<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, InternalFaceHasTwoOwners)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc =
      FaceConnectivityUnstruct<F, S>().build(createTwoAdjacentCubes<F, S>());
  testInternalFaceOwners<F, S>(fc);
}
TYPED_TEST(FaceConnectivityUnstructTest, BoundaryFaceHasNoNeighbor)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  auto fc = FaceConnectivityUnstruct<F, S>().build(createSingleCube<F, S>());
  testBoundaryFaceNoNeighbor<F, S>(fc);
}

}  // namespace
}  // namespace model
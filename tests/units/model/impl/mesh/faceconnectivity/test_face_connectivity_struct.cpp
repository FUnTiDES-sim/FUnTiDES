#include <gtest/gtest.h>

#include "face_connectivity_struct.h"
#include "model_struct.h"
#include "test_face_connectivity_helpers.h"

namespace model
{
namespace
{

template <typename FloatType, typename ScalarType>
ModelStruct<FloatType, ScalarType, 1> createSingleCube()
{
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
ModelStruct<FloatType, ScalarType, 1> createTwoAdjacentCubes()
{
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
class FaceConnectivityStructTest : public ::testing::Test
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
TYPED_TEST_SUITE(FaceConnectivityStructTest, TestTypes);

TYPED_TEST(FaceConnectivityStructTest, BuildFromSingleCube)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testSingleCubeFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, BuildFromTwoAdjacentCubes)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testTwoAdjacentCubesFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, GetGlobalFaceReturnsValidIds)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testGlobalFaceIds<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, GetGlobalFaceUniqueness)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testGlobalFaceUniqueness<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, SingleCubeAllFacesAreBoundary)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testSingleCubeAllBoundary<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, TwoCubesSharedFaceNotBoundary)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testSharedFaceNotBoundary<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, TwoCubesCountBoundaryFaces)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testBoundaryFaceCount<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, FaceNodesValid)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testFaceNodes<F, S>(fc, 8);
}

TYPED_TEST(FaceConnectivityStructTest, FaceNodesUnique)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testFaceNodesUnique<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, InternalFaceHasTwoOwners)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(2, 1, 1, 1);
  testInternalFaceOwners<F, S>(fc);
}

TYPED_TEST(FaceConnectivityStructTest, BoundaryFaceHasNoNeighbor)
{
  using F = typename TestFixture::FloatType;
  using S = typename TestFixture::ScalarType;
  FaceConnectivityStruct<F, S> fc(1, 1, 1, 1);
  testBoundaryFaceNoNeighbor<F, S>(fc);
}

}  // namespace
}  // namespace model
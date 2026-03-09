#include <gtest/gtest.h>

#include "cartesian_unstruct_boundary_classifier.h"

namespace model
{
namespace test
{

// ---------------------------------------------------------------------------
// Helper: fill VECTOR_REAL_VIEW from a std::initializer_list of floats
// ---------------------------------------------------------------------------
static VECTOR_REAL_VIEW makeCoords(std::initializer_list<float> vals)
{
  auto v = allocateVector<VECTOR_REAL_VIEW>(static_cast<int>(vals.size()),
                                             "test_coords");
  int i = 0;
  for (float val : vals) v(i++) = val;
  return v;
}

// domain [0,2]^3, tol = 1e-5
static constexpr float TOL = 1e-5f;

class CartesianUnstructBoundaryClassifierTest : public ::testing::Test
{
 protected:
  // Full-domain classifier: global == local [0,2]^3
  CartesianUnstructBoundaryClassifier<float, int> classifier{
      0.f, 2.f, 0.f, 2.f, 0.f, 2.f, TOL, true};
};

TEST_F(CartesianUnstructBoundaryClassifierTest, InteriorNodeIsInterior)
{
  auto x = makeCoords({1.f});
  auto y = makeCoords({1.f});
  auto z = makeCoords({1.f});

  auto flags = classifier.classify(1, x, y, z);
  EXPECT_EQ(flags(0), static_cast<int>(BoundaryFlag::InteriorNode));
}

TEST_F(CartesianUnstructBoundaryClassifierTest,
       ZmaxNodeIsSurfaceWhenFreeSurfaceOnTop)
{
  auto x = makeCoords({1.f});
  auto y = makeCoords({1.f});
  auto z = makeCoords({2.f});  // at z_max

  auto flags = classifier.classify(1, x, y, z);
  EXPECT_EQ(flags(0), static_cast<int>(BoundaryFlag::Surface));
}

TEST_F(CartesianUnstructBoundaryClassifierTest,
       ZmaxNodeIsDampingWhenFreeSurfaceOff)
{
  CartesianUnstructBoundaryClassifier<float, int> noSurface{
      0.f, 2.f, 0.f, 2.f, 0.f, 2.f, TOL, false};

  auto x = makeCoords({1.f});
  auto y = makeCoords({1.f});
  auto z = makeCoords({2.f});

  auto flags = noSurface.classify(1, x, y, z);
  EXPECT_EQ(flags(0), static_cast<int>(BoundaryFlag::Damping));
}

TEST_F(CartesianUnstructBoundaryClassifierTest, AllSixFacesAreDamping)
{
  // One representative node on each non-zmax face
  auto x = makeCoords({0.f, 2.f, 1.f, 1.f, 1.f});
  auto y = makeCoords({1.f, 1.f, 0.f, 2.f, 1.f});
  auto z = makeCoords({1.f, 1.f, 1.f, 1.f, 0.f});  // last node on z_min

  auto flags = classifier.classify(5, x, y, z);
  for (int n = 0; n < 5; ++n)
    EXPECT_EQ(flags(n), static_cast<int>(BoundaryFlag::Damping)) << "node " << n;
}

TEST_F(CartesianUnstructBoundaryClassifierTest, MixedNodeTypes)
{
  // Layout: interior, zmax (surface), xmin (damping), zmin (damping)
  auto x = makeCoords({1.f, 1.f, 0.f, 1.f});
  auto y = makeCoords({1.f, 1.f, 1.f, 1.f});
  auto z = makeCoords({1.f, 2.f, 1.f, 0.f});

  auto flags = classifier.classify(4, x, y, z);
  EXPECT_EQ(flags(0), static_cast<int>(BoundaryFlag::InteriorNode));
  EXPECT_EQ(flags(1), static_cast<int>(BoundaryFlag::Surface));
  EXPECT_EQ(flags(2), static_cast<int>(BoundaryFlag::Damping));
  EXPECT_EQ(flags(3), static_cast<int>(BoundaryFlag::Damping));
}

// ---------------------------------------------------------------------------
// MPI subdomain test: local mesh in [0,2] x [0,2] x [0,2] with global x in
// [0,4] → node at x=2 is an internal MPI edge, not a physical boundary.
// ---------------------------------------------------------------------------
TEST(CartesianUnstructBoundaryClassifierMpiTest,
     InternalPartitionFaceIsInterior)
{
  CartesianUnstructBoundaryClassifier<float, int> classifier{
      0.f, 4.f,  // global x: local covers only half
      0.f, 2.f,  // global y
      0.f, 2.f,  // global z
      TOL, true};

  // Node at x=2 (local xmax), y=1, z=1 — NOT on any global face
  auto x = makeCoords({2.f});
  auto y = makeCoords({1.f});
  auto z = makeCoords({1.f});

  auto flags = classifier.classify(1, x, y, z);
  EXPECT_EQ(flags(0), static_cast<int>(BoundaryFlag::InteriorNode));
}

TEST(CartesianUnstructBoundaryClassifierMpiTest,
     GlobalXminIsStillDampingInSubdomain)
{
  CartesianUnstructBoundaryClassifier<float, int> classifier{
      0.f, 4.f, 0.f, 2.f, 0.f, 2.f, TOL, true};

  auto x = makeCoords({0.f});
  auto y = makeCoords({1.f});
  auto z = makeCoords({1.f});

  auto flags = classifier.classify(1, x, y, z);
  EXPECT_EQ(flags(0), static_cast<int>(BoundaryFlag::Damping));
}

}  // namespace test
}  // namespace model

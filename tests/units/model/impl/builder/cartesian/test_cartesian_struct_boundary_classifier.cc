#include <gtest/gtest.h>

#include "cartesian_struct_boundary_classifier.h"

namespace model
{
namespace test
{

// ---------------------------------------------------------------------------
// Node index helper: n = k*(nx*ny) + j*nx + i
// ---------------------------------------------------------------------------
static int nodeIndex3(int nx, int ny, int i, int j, int k)
{
  return k * (nx * ny) + j * nx + i;
}

// ---------------------------------------------------------------------------
// Tests for full-domain classification (all six faces are global boundaries)
// ---------------------------------------------------------------------------
// 2×2×2 element mesh, order 1 → 3×3×3 nodes, domain [0, 2]^3
static constexpr int NX = 3, NY = 3, NZ = 3;
static constexpr int N_NODE = NX * NY * NZ;  // 27
static constexpr float L = 2.f;
static constexpr float TOL = L * 1e-4f;

class CartesianStructBoundaryClassifierTest : public ::testing::Test
{
 protected:
  // All six faces are global boundaries; free surface on top
  CartesianStructBoundaryClassifier<float, int> classifier{
      0.f, L, 0.f, L, 0.f, L, TOL, true};

  VECTOR_INT_VIEW classify() const
  {
    return classifier.classify(N_NODE, NX, NY, NZ, 0.f, 0.f, 0.f, L, L, L);
  }
};

TEST_F(CartesianStructBoundaryClassifierTest, InteriorNodeIsInterior)
{
  auto flags = classify();
  // Only (1,1,1) is interior in a 3×3×3 grid
  EXPECT_EQ(flags(nodeIndex3(NX, NY, 1, 1, 1)),
            static_cast<int>(BoundaryFlag::InteriorNode));
}

TEST_F(CartesianStructBoundaryClassifierTest,
       ZmaxNodesAreSurfaceWhenFreeSurfaceOnTop)
{
  auto flags = classify();
  for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i)
      EXPECT_EQ(flags(nodeIndex3(NX, NY, i, j, NZ - 1)),
                static_cast<int>(BoundaryFlag::Surface))
          << "Node (" << i << "," << j << "," << (NZ - 1) << ")";
}

TEST_F(CartesianStructBoundaryClassifierTest,
       ZmaxNodesAreDampingWhenFreeSurfaceOff)
{
  CartesianStructBoundaryClassifier<float, int> noSurface{
      0.f, L, 0.f, L, 0.f, L, TOL, false};
  auto flags =
      noSurface.classify(N_NODE, NX, NY, NZ, 0.f, 0.f, 0.f, L, L, L);
  for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i)
      EXPECT_EQ(flags(nodeIndex3(NX, NY, i, j, NZ - 1)),
                static_cast<int>(BoundaryFlag::Damping))
          << "Node (" << i << "," << j << "," << (NZ - 1) << ")";
}

TEST_F(CartesianStructBoundaryClassifierTest, NonZmaxBoundaryNodesAreDamping)
{
  auto flags = classify();
  // xmin face (i=0), k < NZ-1
  for (int k = 0; k < NZ - 1; ++k)
    for (int j = 0; j < NY; ++j)
      EXPECT_EQ(flags(nodeIndex3(NX, NY, 0, j, k)),
                static_cast<int>(BoundaryFlag::Damping));
  // xmax face (i=NX-1), k < NZ-1
  for (int k = 0; k < NZ - 1; ++k)
    for (int j = 0; j < NY; ++j)
      EXPECT_EQ(flags(nodeIndex3(NX, NY, NX - 1, j, k)),
                static_cast<int>(BoundaryFlag::Damping));
}

TEST_F(CartesianStructBoundaryClassifierTest, NodeCountsByFlag)
{
  auto flags = classify();
  int n_interior = 0, n_surface = 0, n_damping = 0;
  for (int n = 0; n < N_NODE; ++n)
  {
    if (flags(n) == static_cast<int>(BoundaryFlag::InteriorNode))
      ++n_interior;
    else if (flags(n) == static_cast<int>(BoundaryFlag::Surface))
      ++n_surface;
    else if (flags(n) == static_cast<int>(BoundaryFlag::Damping))
      ++n_damping;
  }
  // 3×3×3 = 27 nodes: 1 interior, 9 surface (z=max), 17 damping
  EXPECT_EQ(n_interior, 1);
  EXPECT_EQ(n_surface, 9);
  EXPECT_EQ(n_damping, 17);
}

// ---------------------------------------------------------------------------
// MPI subdomain: local mesh occupies [0, L]^3 but the global x domain is
// [0, 2*L], so the local x_max face is an internal MPI edge.
// ---------------------------------------------------------------------------
TEST(CartesianStructBoundaryClassifierMpiTest, InternalPartitionFaceIsInterior)
{
  // Global x extends to 2*L; local x_max (= L) is NOT a global boundary
  CartesianStructBoundaryClassifier<float, int> classifier{
      0.f, 2 * L, 0.f, L, 0.f, L, TOL, true};
  auto flags =
      classifier.classify(N_NODE, NX, NY, NZ, 0.f, 0.f, 0.f, L, L, L);

  // Node at i=NX-1, j=1, k=1: x_max is internal → interior
  EXPECT_EQ(flags(nodeIndex3(NX, NY, NX - 1, 1, 1)),
            static_cast<int>(BoundaryFlag::InteriorNode));

  // Node at i=0, j=1, k=1: x_min is global → damping
  EXPECT_EQ(flags(nodeIndex3(NX, NY, 0, 1, 1)),
            static_cast<int>(BoundaryFlag::Damping));

  // Node at i=NX-1, j=1, k=NZ-1: z_max is global → surface
  EXPECT_EQ(flags(nodeIndex3(NX, NY, NX - 1, 1, NZ - 1)),
            static_cast<int>(BoundaryFlag::Surface));
}

}  // namespace test
}  // namespace model

#include <gtest/gtest.h>

#include "cartesian_struct_boundary_classifier.h"
#include "model_struct.h"

namespace model
{
namespace test
{

// ---------------------------------------------------------------------------
// Helper: build a minimal ModelStructData for a domain of 'ne x ne x ne'
// order-1 elements spanning [0, L]^3 with origin at (ox, oy, oz).
// ---------------------------------------------------------------------------
static model::ModelStructData<float, int> makeStructData(int ne, float L,
                                                          float ox = 0.f,
                                                          float oy = 0.f,
                                                          float oz = 0.f)
{
  model::ModelStructData<float, int> data;
  data.ex_ = ne;
  data.ey_ = ne;
  data.ez_ = ne;
  data.dx_ = L;
  data.dy_ = L;
  data.dz_ = L;
  data.ox_ = ox;
  data.oy_ = oy;
  data.oz_ = oz;
  data.isModelOnNodes_ = true;
  data.isElastic_ = false;
  return data;
}

// ---------------------------------------------------------------------------
// Node index helpers for Order=1, ne elements per side → (ne+1)^3 nodes
// Ordering: n = k*(nx*ny) + j*nx + i
// ---------------------------------------------------------------------------
static int nodeIndex3(int nx, int ny, int i, int j, int k)
{
  return k * (nx * ny) + j * nx + i;
}

// ---------------------------------------------------------------------------
// Tests for full-domain classification (global == local bounds)
// ---------------------------------------------------------------------------
class CartesianStructBoundaryClassifierTest : public ::testing::Test
{
 protected:
  // 2x2x2 element mesh, each element has unit length → domain [0,2]^3
  static constexpr int ne = 2;
  static constexpr float L = 2.f;
  static constexpr int nx = ne + 1;  // 3 nodes per side for Order=1
  static constexpr int ny = ne + 1;
  static constexpr int nz = ne + 1;

  model::ModelStruct<float, int, 1> mesh{makeStructData(ne, L)};

  // Classifier with global == local bounds
  CartesianStructBoundaryClassifier<float, int, 1> classifier{0.f, L, 0.f, L,
                                                               0.f, L, true};
};

TEST_F(CartesianStructBoundaryClassifierTest, InteriorNodeIsInterior)
{
  // Center node (1,1,1) is the only fully interior node
  const int n_interior = nodeIndex3(nx, ny, 1, 1, 1);
  auto flags = classifier.classify(mesh);
  EXPECT_EQ(flags(n_interior), static_cast<int>(BoundaryFlag::InteriorNode));
}

TEST_F(CartesianStructBoundaryClassifierTest,
       ZmaxNodesAreSurfaceWhenFreeSurfaceOnTop)
{
  auto flags = classifier.classify(mesh);
  // All nodes with k=nz-1 (z=2) should be Surface
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      EXPECT_EQ(flags(nodeIndex3(nx, ny, i, j, nz - 1)),
                static_cast<int>(BoundaryFlag::Surface))
          << "Node (" << i << "," << j << "," << (nz - 1) << ")";
}

TEST_F(CartesianStructBoundaryClassifierTest,
       ZmaxNodesAreDampingWhenFreeSurfaceOff)
{
  CartesianStructBoundaryClassifier<float, int, 1> noSurface{0.f, L, 0.f, L,
                                                              0.f, L, false};
  auto flags = noSurface.classify(mesh);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      EXPECT_EQ(flags(nodeIndex3(nx, ny, i, j, nz - 1)),
                static_cast<int>(BoundaryFlag::Damping))
          << "Node (" << i << "," << j << "," << (nz - 1) << ")";
}

TEST_F(CartesianStructBoundaryClassifierTest, NonZmaxBoundaryNodesAreDamping)
{
  auto flags = classifier.classify(mesh);
  // xmin face (i=0) with k < nz-1
  for (int k = 0; k < nz - 1; ++k)
    for (int j = 0; j < ny; ++j)
      EXPECT_EQ(flags(nodeIndex3(nx, ny, 0, j, k)),
                static_cast<int>(BoundaryFlag::Damping));
  // xmax face (i=nx-1) with k < nz-1
  for (int k = 0; k < nz - 1; ++k)
    for (int j = 0; j < ny; ++j)
      EXPECT_EQ(flags(nodeIndex3(nx, ny, nx - 1, j, k)),
                static_cast<int>(BoundaryFlag::Damping));
}

TEST_F(CartesianStructBoundaryClassifierTest, NodeCountsByFlag)
{
  auto flags = classifier.classify(mesh);
  int n_interior = 0, n_surface = 0, n_damping = 0;
  for (int n = 0; n < mesh.getNumberOfNodes(); ++n)
  {
    if (flags(n) == static_cast<int>(BoundaryFlag::InteriorNode)) ++n_interior;
    else if (flags(n) == static_cast<int>(BoundaryFlag::Surface)) ++n_surface;
    else if (flags(n) == static_cast<int>(BoundaryFlag::Damping)) ++n_damping;
  }
  // 3x3x3 = 27 nodes total
  // Interior: only (1,1,1) → 1
  // Surface (z=2): 3x3 = 9
  // Damping: all others = 17
  EXPECT_EQ(n_interior, 1);
  EXPECT_EQ(n_surface, 9);
  EXPECT_EQ(n_damping, 17);
}

// ---------------------------------------------------------------------------
// MPI subdomain test: local mesh occupies [0, 2]^3 but global domain is
// [0, 4] x [0, 2] x [0, 2].  The face at x=2 is an internal MPI edge and
// must NOT be classified as a physical boundary.
// ---------------------------------------------------------------------------
TEST(CartesianStructBoundaryClassifierMpiTest, InternalPartitionFaceIsInterior)
{
  constexpr int ne = 2;
  constexpr float L = 2.f;
  constexpr int nx = ne + 1;
  constexpr int ny = ne + 1;
  constexpr int nz = ne + 1;

  model::ModelStruct<float, int, 1> mesh{makeStructData(ne, L)};

  // Global domain is twice as wide in X
  CartesianStructBoundaryClassifier<float, int, 1> classifier{
      0.f, 4.f,   // global x bounds (local only covers [0,2])
      0.f, L,     // global y bounds (same as local)
      0.f, L,     // global z bounds
      true};

  auto flags = classifier.classify(mesh);

  // Nodes on the local x_max face (i = nx-1 = 2) that are NOT on any other
  // global boundary should be InteriorNode (they sit on an internal MPI edge)
  // — specifically, node (2, 1, 1) is not on any global face.
  const int n = nodeIndex3(nx, ny, nx - 1, 1, 1);
  EXPECT_EQ(flags(n), static_cast<int>(BoundaryFlag::InteriorNode));

  // But node (0, 1, 1) is still on the global x_min → Damping
  const int n_xmin = nodeIndex3(nx, ny, 0, 1, 1);
  EXPECT_EQ(flags(n_xmin), static_cast<int>(BoundaryFlag::Damping));

  // Node (2, 1, 2) is on the global z_max → Surface
  const int n_surf = nodeIndex3(nx, ny, nx - 1, 1, nz - 1);
  EXPECT_EQ(flags(n_surf), static_cast<int>(BoundaryFlag::Surface));
}

}  // namespace test
}  // namespace model

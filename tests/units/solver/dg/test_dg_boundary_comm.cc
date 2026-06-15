/**
 * @file test_dg_boundary_comm.cc
 * @brief Unit tests for DGBoundaryComm (single-rank path; no MPI runtime required).
 */
#include <gtest/gtest.h>

#include "data_type.h"
#include "dg_boundary_comm.h"

namespace solver {
namespace fe {
namespace test {

TEST(DGBoundaryComm, SingleRank_NotActive) {
  DGBoundaryComm comm;
  comm.setup(0, 1, 3, 2, 2, 8, 12);
  EXPECT_FALSE(comm.isActive());
}

TEST(DGBoundaryComm, MultiRank_IsActive) {
  DGBoundaryComm comm;
  comm.setup(0, 2, 3, 2, 2, 8, 12);
  EXPECT_TRUE(comm.isActive());
}

TEST(DGBoundaryComm, Setup_IndicesCorrect) {
  // rank=0, size=2, ex_local=2, ey=2, ez=2, n_dof=8, n_local=8, n_ghost=4
  // left_elems[m]   = j*ex_local + k*ex_local*ey
  // left_ghosts[m]  = n_local + k*ey + j
  // right_elems[m]  = (ex_local-1) + j*ex_local + k*ex_local*ey
  // right_ghosts[m] = n_local + n_ghost + k*ey + j
  DGBoundaryComm comm;
  comm.setup(0, 2, 2, 2, 2, 8, 8);

  const auto& le = comm.leftElems();
  const auto& lg = comm.leftGhosts();
  const auto& re = comm.rightElems();
  const auto& rg = comm.rightGhosts();

  ASSERT_EQ(static_cast<int>(le.size()), 4);

  // m=0: k=0, j=0
  EXPECT_EQ(le[0], 0);
  EXPECT_EQ(lg[0], 8);
  EXPECT_EQ(re[0], 1);
  EXPECT_EQ(rg[0], 12);

  // m=1: k=0, j=1
  EXPECT_EQ(le[1], 2);
  EXPECT_EQ(lg[1], 9);
  EXPECT_EQ(re[1], 3);
  EXPECT_EQ(rg[1], 13);

  // m=2: k=1, j=0
  EXPECT_EQ(le[2], 4);
  EXPECT_EQ(lg[2], 10);
  EXPECT_EQ(re[2], 5);
  EXPECT_EQ(rg[2], 14);

  // m=3: k=1, j=1
  EXPECT_EQ(le[3], 6);
  EXPECT_EQ(lg[3], 11);
  EXPECT_EQ(re[3], 7);
  EXPECT_EQ(rg[3], 15);
}

TEST(DGBoundaryComm, Exchange_SingleRank_NoOp) {
  DGBoundaryComm comm;
  comm.setup(0, 1, 2, 2, 2, 8, 8);
  arrayReal field("field", 8, 8);
  Kokkos::deep_copy(field, 42.f);
  comm.exchange(field);
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field);
  EXPECT_FLOAT_EQ(h(0, 0), 42.f);
}

}  // namespace test
}  // namespace fe
}  // namespace solver

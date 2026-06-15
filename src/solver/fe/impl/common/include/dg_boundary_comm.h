#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_BOUNDARY_COMM_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_BOUNDARY_COMM_H_

#include <vector>

#include "data_type.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace solver {
namespace fe {

/**
 * @brief Handles ghost element row exchange for DG 1D X-decomposition.
 *
 * Ghost layout (field has n_local + 2*n_ghost rows):
 *   rows [0,          n_local)            local elements
 *   rows [n_local,    n_local+n_ghost)     left ghost  (rank-1 right slice)
 *   rows [n_local+n_ghost, n_local+2*n_ghost) right ghost (rank+1 left slice)
 *
 * Exchange is a two-phase paired Sendrecv: no-op when size==1.
 */
class DGBoundaryComm {
 public:
  /**
   * @brief Configure the communicator for the local partition.
   *
   * @param rank          MPI rank of this process.
   * @param size          Total number of MPI ranks.
   * @param ex_local      Local element count in X.
   * @param ey            Element count in Y (global, same on all ranks).
   * @param ez            Element count in Z (global, same on all ranks).
   * @param n_dof         DOF per element.
   * @param n_local       Total local element count (ex_local * ey * ez).
   */
  void setup(int rank, int size, int ex_local, int ey, int ez, int n_dof, int n_local) {
    rank_ = rank;
    size_ = size;
    n_local_ = n_local;
    n_ghost_ = ey * ez;
    n_dof_ = n_dof;

    left_elems_.resize(n_ghost_);
    left_ghosts_.resize(n_ghost_);
    right_elems_.resize(n_ghost_);
    right_ghosts_.resize(n_ghost_);

    int m = 0;
    for (int k = 0; k < ez; ++k) {
      for (int j = 0; j < ey; ++j) {
        left_elems_[m] = j * ex_local + k * ex_local * ey;
        left_ghosts_[m] = n_local + k * ey + j;
        right_elems_[m] = (ex_local - 1) + j * ex_local + k * ex_local * ey;
        right_ghosts_[m] = n_local + n_ghost_ + k * ey + j;
        ++m;
      }
    }
  }

  bool isActive() const { return size_ > 1; }

  const std::vector<int>& leftElems() const { return left_elems_; }
  const std::vector<int>& leftGhosts() const { return left_ghosts_; }
  const std::vector<int>& rightElems() const { return right_elems_; }
  const std::vector<int>& rightGhosts() const { return right_ghosts_; }

  /**
   * @brief Exchange boundary DOF slices with MPI neighbours.
   *
   * Fills ghost rows in @p field from neighbouring ranks.
   * No-op when size==1.
   *
   * @param field  2D Kokkos array of shape [n_local + 2*n_ghost, n_dof].
   */
  void exchange(arrayReal& field) const {
    if (size_ == 1) return;
#ifdef USE_MPI
    auto h_field = Kokkos::create_mirror_view(field);
    Kokkos::deep_copy(h_field, field);

    const int left_rank = (rank_ > 0) ? rank_ - 1 : MPI_PROC_NULL;
    const int right_rank = (rank_ < size_ - 1) ? rank_ + 1 : MPI_PROC_NULL;
    const int buf_size = n_ghost_ * n_dof_;

    std::vector<float> send_right(buf_size), recv_left(buf_size, 0.f);
    std::vector<float> send_left(buf_size), recv_right(buf_size, 0.f);

    for (int m = 0; m < n_ghost_; ++m)
      for (int d = 0; d < n_dof_; ++d) send_right[m * n_dof_ + d] = h_field(right_elems_[m], d);

    for (int m = 0; m < n_ghost_; ++m)
      for (int d = 0; d < n_dof_; ++d) send_left[m * n_dof_ + d] = h_field(left_elems_[m], d);

    // Phase 1: send right slice → right_rank, recv from left_rank → left ghost
    MPI_Sendrecv(send_right.data(), buf_size, MPI_FLOAT, right_rank, 0, recv_left.data(), buf_size, MPI_FLOAT,
                 left_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Phase 2: send left slice → left_rank, recv from right_rank → right ghost
    MPI_Sendrecv(send_left.data(), buf_size, MPI_FLOAT, left_rank, 1, recv_right.data(), buf_size, MPI_FLOAT,
                 right_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (rank_ > 0)
      for (int m = 0; m < n_ghost_; ++m)
        for (int d = 0; d < n_dof_; ++d) h_field(n_local_ + m, d) = recv_left[m * n_dof_ + d];

    if (rank_ < size_ - 1)
      for (int m = 0; m < n_ghost_; ++m)
        for (int d = 0; d < n_dof_; ++d) h_field(n_local_ + n_ghost_ + m, d) = recv_right[m * n_dof_ + d];

    Kokkos::deep_copy(field, h_field);
#endif
  }

 private:
  int rank_ = 0, size_ = 1;
  int n_local_ = 0, n_ghost_ = 0, n_dof_ = 0;
  std::vector<int> left_elems_, left_ghosts_, right_elems_, right_ghosts_;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_DG_BOUNDARY_COMM_H_

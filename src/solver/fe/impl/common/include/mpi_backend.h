#ifndef SRC_SOLVER_FE_IMPL_INCLUDE_MPI_BACKEND_H_
#define SRC_SOLVER_FE_IMPL_INCLUDE_MPI_BACKEND_H_

#include <mpi.h>

#include <map>
#include <stdexcept>
#include <vector>

#include "boundary_synchronizer.h"

namespace solver {
namespace fe {

/**
 * @brief Boundary exchange backend based on non-blocking MPI point-to-point calls.
 *
 * Every neighbor rank receives the buffer stored under its rank in the send map and
 * is expected to send back a buffer of the same size. Uses MPI_COMM_WORLD.
 */
class MPIBackend : public BoundarySynchronizer::Backend {
 public:
  MPIBackend() = default;
  ~MPIBackend() = default;

  /**
   * @brief Exchanges one buffer with each neighbor rank and blocks until all transfers complete.
   * @param[in] sendBuffers Data to send, keyed by neighbor rank.
   * @param[out] recvBuffers Cleared, then filled with one buffer per neighbor rank, each
   *             of the same size as the corresponding send buffer.
   * @throws std::runtime_error If MPI_Waitall fails.
   */
  void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                std::map<int, std::vector<float>>& recvBuffers) override {
    recvBuffers.clear();
    std::vector<MPI_Request> requests;
    requests.reserve(sendBuffers.size() * 2);

    // Receives are posted first; the received size is assumed equal to the sent size
    // (symmetric partition).
    for (const auto& [neighborRank, sendData] : sendBuffers) {
      recvBuffers[neighborRank].resize(sendData.size());

      MPI_Request req;
      int tag = 0;
      MPI_Irecv(recvBuffers[neighborRank].data(), static_cast<int>(sendData.size()), MPI_FLOAT, neighborRank, tag,
                MPI_COMM_WORLD, &req);
      requests.push_back(req);
    }

    for (const auto& [neighborRank, sendData] : sendBuffers) {
      MPI_Request req;
      int tag = 0;
      // MPI-3 signatures take a non-const send buffer, but it is only read.
      MPI_Isend(const_cast<float*>(sendData.data()), static_cast<int>(sendData.size()), MPI_FLOAT, neighborRank, tag,
                MPI_COMM_WORLD, &req);
      requests.push_back(req);
    }

    if (!requests.empty()) {
      std::vector<MPI_Status> statuses(requests.size());
      int err = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), statuses.data());
      if (err != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Waitall failed during boundary exchange");
      }
    }
  }
};

}  // namespace fe
}  // namespace solver

#endif  // SRC_SOLVER_FE_IMPL_INCLUDE_MPI_BACKEND_H_

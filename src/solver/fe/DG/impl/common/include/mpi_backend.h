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
 * @brief MPI implementation for boundary synchronization.
 * Uses Non-blocking communication (Isend/Irecv) to exchange ghost layers.
 */
class MPIBackend : public BoundarySynchronizer::Backend {
 public:
  MPIBackend() = default;
  ~MPIBackend() = default;

  void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                std::map<int, std::vector<float>>& recvBuffers) override {
    recvBuffers.clear();
    std::vector<MPI_Request> requests;
    // Reserve space for 2 requests per neighbor (send + recv)
    requests.reserve(sendBuffers.size() * 2);

    // Post Non-Blocking Receives
    // /!\ We assume we receive the same amount we send (symmetric mesh
    // partition)
    for (const auto& [neighborRank, sendData] : sendBuffers) {
      // Allocate receive buffer
      recvBuffers[neighborRank].resize(sendData.size());

      MPI_Request req;
      int tag = 0;  // Simple tag
      MPI_Irecv(recvBuffers[neighborRank].data(), static_cast<int>(sendData.size()), MPI_FLOAT, neighborRank, tag,
                MPI_COMM_WORLD, &req);
      requests.push_back(req);
    }

    // Post Non-Blocking Sends
    for (const auto& [neighborRank, sendData] : sendBuffers) {
      MPI_Request req;
      int tag = 0;
      // const_cast is safe here because MPI_Isend treats buffer as const
      MPI_Isend(const_cast<float*>(sendData.data()), static_cast<int>(sendData.size()), MPI_FLOAT, neighborRank, tag,
                MPI_COMM_WORLD, &req);
      requests.push_back(req);
    }

    // Wait for all communications to finish
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

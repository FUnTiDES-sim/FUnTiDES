#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_BOUNDARY_SYNCHRONIZER_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_BOUNDARY_SYNCHRONIZER_H_
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include "parallel_topology.h"

using namespace utils;

namespace solver {
namespace fe {
/**
 * @brief Sums the values of a nodal field at the nodes shared between partitions.
 *
 * Each rank packs its values at the shared nodes, hands them to a pluggable
 * Backend that carries them to the neighbor ranks, and adds the values
 * received from the neighbors to its own. Used to complete the mass matrix
 * (once after assembly) and the force vector (every time step, after
 * computeForces() and before the solution update).
 *
 * The synchronizer is driven by the caller, not by the solver. Values travel
 * as float, whatever the type of the field. When the topology is not
 * distributed, synchronize() does nothing.
 */
class BoundarySynchronizer {
 public:
  /**
   * @brief Communication strategy used to exchange shared-node values between ranks.
   */
  struct Backend {
    virtual ~Backend() = default;

    /**
     * @brief Sends one buffer to each neighbor rank and receives one from each.
     *
     * The receive buffers must hold, for each neighbor, as many values as the
     * send buffer for that neighbor, in the same node order.
     *
     * @param[in] sendBuffers Values to send, keyed by neighbor rank.
     * @param[out] recvBuffers Values received, keyed by neighbor rank. The
     *             implementation replaces any previous content.
     *
     * @throws std::runtime_error if the communication fails.
     */
    virtual void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                          std::map<int, std::vector<float>>& recvBuffers) = 0;
  };

  /**
   * @brief Builds a synchronizer that takes ownership of a backend.
   *
   * @param[in] backend Communication strategy, must not be null.
   *
   * @throws std::invalid_argument if backend is null.
   */
  explicit BoundarySynchronizer(std::unique_ptr<Backend> backend) : m_backend(std::move(backend)) {
    if (!m_backend) {
      throw std::invalid_argument("Backend cannot be null");
    }
  }

  /**
   * @brief Adds to a nodal field the values held by the neighbor ranks at the shared nodes.
   *
   * After the call, each shared node holds the sum of the local value and the
   * values of all ranks sharing it. Does nothing if the topology is not
   * distributed.
   *
   * @tparam ViewType Type whose operator()(int) gives read and write access to
   *         the value at a node index (for example a Kokkos::View). The
   *         storage must be accessible from the host.
   *
   * @param[in,out] field Nodal field, indexed by local node index.
   * @param[in] topo Topology giving, for each neighbor rank, the local indices
   *            of the shared nodes.
   *
   * @throws std::runtime_error if the exchange or the accumulation fails
   *         (including a buffer size mismatch or missing neighbor data).
   */
  template <typename ViewType>
  void synchronize(ViewType& field, const ParallelTopology& topo) {
    if (!topo.isDistributed()) {
      return;
    }

    try {
      auto sendBufs = pack(field, topo);

      std::map<int, std::vector<float>> recvBufs;
      m_backend->exchange(sendBufs, recvBufs);

      accumulate(field, recvBufs, topo);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Boundary synchronization failed: ") + e.what());
    }
  }

 private:
  std::unique_ptr<Backend> m_backend;  ///< Communication strategy, never null.

  /**
   * @brief Copies the field values at the shared nodes into one buffer per neighbor rank.
   *
   * @tparam ViewType Type whose operator()(int) gives the value at a node index.
   * @param[in] field Nodal field.
   * @param[in] topo Topology giving the shared nodes.
   *
   * @return Map from neighbor rank to values, in the order of
   *         topo.sharedNodes for that rank.
   */
  template <typename ViewType>
  static std::map<int, std::vector<float>> pack(const ViewType& field, const ParallelTopology& topo) {
    std::map<int, std::vector<float>> buffers;

    for (const auto& [neighborRank, nodeIndices] : topo.sharedNodes) {
      auto& buf = buffers[neighborRank];
      buf.reserve(nodeIndices.size());

      for (int nodeIdx : nodeIndices) {
        buf.push_back(static_cast<float>(field(nodeIdx)));
      }
    }

    return buffers;
  }

  /**
   * @brief Adds the values received from each neighbor rank to the field at the shared nodes.
   *
   * @tparam ViewType Type whose operator()(int) gives read and write access to
   *         the value at a node index.
   * @param[in,out] field Nodal field.
   * @param[in] recvBufs Received values, keyed by neighbor rank, in the order
   *            of topo.sharedNodes for that rank.
   * @param[in] topo Topology giving the shared nodes.
   *
   * @throws std::length_error if a received buffer size differs from the
   *         number of shared nodes with that rank.
   * @throws std::runtime_error if the topology is distributed and no buffer
   *         was received from an expected neighbor.
   */
  template <typename ViewType>
  static void accumulate(ViewType& field, const std::map<int, std::vector<float>>& recvBufs,
                         const ParallelTopology& topo) {
    for (const auto& [neighborRank, nodeIndices] : topo.sharedNodes) {
      auto it = recvBufs.find(neighborRank);

      if (it == recvBufs.end()) {
        if (topo.isDistributed()) {
          throw std::runtime_error("Expected data from rank " + std::to_string(neighborRank) +
                                   " but received nothing. Exchange failed or topology mismatch.");
        }
        continue;
      }

      const auto& buf = it->second;

      if (buf.size() != nodeIndices.size()) {
        throw std::length_error("Buffer size mismatch from rank " + std::to_string(neighborRank) + ": expected " +
                                std::to_string(nodeIndices.size()) + " values, got " + std::to_string(buf.size()));
      }

      for (size_t i = 0; i < nodeIndices.size(); ++i) {
        int nodeIdx = nodeIndices[i];
        field(nodeIdx) += buf[i];
      }
    }
  }
};

/**
 * @brief Backend for a single rank: performs no communication.
 */
class SerialBackend : public BoundarySynchronizer::Backend {
 public:
  /**
   * @brief Empties recvBuffers, since there is no neighbor to receive from.
   */
  void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                std::map<int, std::vector<float>>& recvBuffers) override {
    recvBuffers.clear();
  }
};

/**
 * @brief Backend that logs the exchanges to stdout without communicating.
 *
 * For each neighbor rank it prints the number of values sent and the number
 * "received", and fills the receive buffer with zeros, so the accumulation
 * leaves the field unchanged. Meant for testing the distributed call sequence
 * without MPI.
 */
class DebugBackend : public BoundarySynchronizer::Backend {
 private:
  int m_rank;  ///< Rank printed in the log lines.

 public:
  /**
   * @param[in] rank Rank printed in the log lines.
   */
  explicit DebugBackend(int rank = 0) : m_rank(rank) {}

  /**
   * @brief Logs the send sizes and sets each receive buffer to zeros of the same size.
   */
  void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                std::map<int, std::vector<float>>& recvBuffers) override {
    std::cout << "[Rank " << m_rank << "] Boundary Synchronization:\n";

    for (const auto& [neighbor, data] : sendBuffers) {
      std::cout << "  → Send " << data.size() << " values to rank " << neighbor << "\n";
    }

    recvBuffers.clear();
    for (const auto& [neighbor, data] : sendBuffers) {
      recvBuffers[neighbor].resize(data.size(), 0.0f);
      std::cout << "  ← Recv " << data.size() << " values from rank " << neighbor << " (zeroed for debug)\n";
    }
  }
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_BOUNDARY_SYNCHRONIZER_H_

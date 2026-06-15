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
// @brief Handles data exchange at partition boundaries.
//
// Manages synchronization of boundary node values across distributed ranks.
// Uses a pluggable Backend pattern to support different communication
// strategies (serial, debug, MPI, GPU direct, etc.).
//
// @details
// The synchronization process has three phases:
// 1. **Pack**: Extract values from boundary nodes into buffers
// 2. **Exchange**: Send buffers to neighbors, receive from neighbors (via
// Backend)
// 3. **Accumulate**: Sum received values with local values (assembly semantics)
//
// This is **decoupled from the solver**: the external orchestrator decides
// when and how to synchronize. This enables Python integration, adaptive
// refinement, and custom communication strategies.
//
// @note
// For single-rank execution, synchronization is a no-op (no neighbors exist).
// For distributed execution, must be called after computeForces() and before
// updateSolutionForward() or updateSolutionBackward() to ensure boundary values are complete.
class BoundarySynchronizer {
 public:
  // @brief Controls how received values are applied to the local field.
  //
  // SUM:  field(idx) += recv  — SEM force/mass assembly (partial contributions)
  // COPY: field(idx) =  recv  — DG wavefield exchange (overwrite with neighbor values)
  enum class SyncMode { SUM, COPY };

  // @brief Backend interface for communication.
  //
  // Implementations can provide different communication strategies:
  // - SerialBackend: no-op (single rank)
  // - DebugBackend: logs all exchanges for diagnostics
  // - MPIBackend: real MPI communication (future)
  // - GPUBackend: GPU direct peer-to-peer (future)
  struct Backend {
    virtual ~Backend() = default;

    // @brief Exchange boundary data with neighboring ranks.
    //
    // Packs send buffers from local rank, sends to neighbors, and receives
    // from neighbors. The actual mechanism (MPI, loopback, GPU DMA, etc.)
    // is implementation-specific.
    //
    // @param[in] sendBuffers Data to send (neighbor rank -> values)
    // @param[out] recvBuffers Received data from neighbors
    //
    // @note
    // Implementation is responsible for:
    // - Clearing or resizing recvBuffers as appropriate
    // - Handling communication failures (throwing std::runtime_error)
    // - Ensuring buffer order matches sendBuffers (same neighbor ranks)
    virtual void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                          std::map<int, std::vector<float>>& recvBuffers) = 0;
  };

  // @brief Constructor.
  //
  // @param[in] backend Communication strategy (caller manages lifetime)
  //
  // @throws std::invalid_argument if backend is null
  explicit BoundarySynchronizer(std::unique_ptr<Backend> backend) : m_backend(std::move(backend)) {
    if (!m_backend) {
      throw std::invalid_argument("Backend cannot be null");
    }
  }

  // @brief Synchronize a field across partition boundaries.
  //
  // @tparam ViewType Type providing operator()(int) access to field values
  //                  (e.g., std::vector, Kokkos::View)
  //
  // @param[in,out] field Field to synchronize (modified in-place)
  // @param[in] topo Topology describing shared nodes
  // @param[in] mode SUM (default): accumulate partial contributions (SEM forces/mass);
  //                 COPY: overwrite with neighbor values (DG wavefield exchange)
  //
  // @throws std::runtime_error if exchange or accumulation fails
  template <typename ViewType>
  void synchronize(ViewType& field, const ParallelTopology& topo,
                   SyncMode mode = SyncMode::SUM) {
    if (!topo.isDistributed())
      return;

    try {
      auto sendBufs = pack(field, topo);
      std::map<int, std::vector<float>> recvBufs;
      m_backend->exchange(sendBufs, recvBufs);
      accumulate(field, recvBufs, topo, mode);
    } catch (const std::exception& e) {
      throw std::runtime_error(std::string("Boundary synchronization failed: ") + e.what());
    }
  }

 private:
  std::unique_ptr<Backend> m_backend;  //< Communication strategy

  // @brief Extract boundary node values into buffers keyed by neighbor rank.
  //
  // @tparam ViewType Type providing operator()(int) access
  // @param[in] field Field to extract from
  // @param[in] topo Topology describing which nodes are shared
  //
  // @return Map from neighbor rank to vector of extracted values
  //         (in same order as topo.sharedNodes[rank])
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

  // @brief Apply received values into field at boundary nodes.
  //
  // @tparam ViewType Type providing operator()(int) access
  // @param[in,out] field Field to update (modified in-place)
  // @param[in] recvBufs Received data from neighbors
  // @param[in] topo Topology describing shared nodes
  // @param[in] mode SUM: field(idx) += recv; COPY: field(idx) = recv
  //
  // @throws std::length_error if received buffer size != expected node count
  // @throws std::runtime_error if expected neighbor data is missing
  template <typename ViewType>
  static void accumulate(ViewType& field, const std::map<int, std::vector<float>>& recvBufs,
                         const ParallelTopology& topo, SyncMode mode) {
    for (const auto& [neighborRank, nodeIndices] : topo.sharedNodes) {
      auto it = recvBufs.find(neighborRank);

      if (it == recvBufs.end()) {
        if (topo.isDistributed())
          throw std::runtime_error("Expected data from rank " + std::to_string(neighborRank) +
                                   " but received nothing. Exchange failed or topology mismatch.");
        continue;
      }

      const auto& buf = it->second;

      if (buf.size() != nodeIndices.size())
        throw std::length_error("Buffer size mismatch from rank " + std::to_string(neighborRank) + ": expected " +
                                std::to_string(nodeIndices.size()) + " values, got " + std::to_string(buf.size()));

      for (size_t i = 0; i < nodeIndices.size(); ++i) {
        int nodeIdx = nodeIndices[i];
        if (mode == SyncMode::COPY)
          field(nodeIdx) = buf[i];
        else
          field(nodeIdx) += buf[i];
      }
    }
  }
};

// @brief Serial backend: No-op communication.
//
// Used when running on a single rank or in serial testing mode.
// Simply clears receive buffers to prevent spurious data accumulation.
class SerialBackend : public BoundarySynchronizer::Backend {
 public:
  // @brief Exchange (no-op for serial execution).
  //
  // Clears recvBuffers since there are no neighbors in serial mode.
  void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                std::map<int, std::vector<float>>& recvBuffers) override {
    // Single rank: clear recv buffers (no data from neighbors)
    recvBuffers.clear();
  }
};

// @brief Debug backend: Logs synchronization for diagnostics.
//
// Useful for testing distributed logic without MPI. Prints all exchange
// operations and populates receive buffers with zeros (safe default for
// testing assembly logic).
//
// Output format:
// ```
// [Rank N] Boundary Synchronization:
//   → Send M values to rank X
//   ← Recv M values from rank X (zeroed for debug)
// ```
class DebugBackend : public BoundarySynchronizer::Backend {
 private:
  int m_rank;  //< Current rank (for logging)

 public:
  // @brief Constructor.
  // @param[in] rank Current MPI rank (default 0)
  explicit DebugBackend(int rank = 0) : m_rank(rank) {}

  // @brief Exchange (with logging).
  //
  // Logs all send/receive operations and populates recvBuffers with zeros.
  // Useful for verifying synchronization is called at correct times.
  void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                std::map<int, std::vector<float>>& recvBuffers) override {
    std::cout << "[Rank " << m_rank << "] Boundary Synchronization:\n";

    for (const auto& [neighbor, data] : sendBuffers) {
      std::cout << "  → Send " << data.size() << " values to rank " << neighbor << "\n";
    }

    // In debug mode, populate recv with zeros (safe default)
    // In real MPI, this would receive actual data from neighbors
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

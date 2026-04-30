#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_BOUNDARY_SYNCHRONIZER_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_BOUNDARY_SYNCHRONIZER_H_
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include "parallel_topology.h"

// This do not need to be duplicated but I need ti check inside the solver to not break mpi

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
// updateSolution() to ensure boundary values are complete.
class BoundarySynchronizer {
 public:
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

  // @brief Synchronize a field (mass matrix or forces).
  //
  // For distributed execution, extracts boundary node values, exchanges with
  // neighbors, and accumulates received values (summing). This is used for:
  // - Mass matrix assembly: sums partial mass contributions at boundaries
  // - Force assembly: sums partial stiffness contributions at boundaries
  //
  // @tparam ViewType Type providing operator()(int) access to field values
  //                  (e.g., std::vector, Kokkos::View)
  //
  // @param[in,out] field Field to synchronize (modified in-place)
  // @param[in] topo Topology describing shared nodes
  //
  // @throws std::runtime_error if exchange or accumulation fails
  //
  // @details
  // Process:
  // 1. If not distributed, return immediately (no-op)
  // 2. Pack boundary values from field
  // 3. Exchange with neighbors via backend
  // 4. Accumulate (sum) received values back into field
  //
  // After synchronization, boundary node values contain contributions from
  // all neighboring elements, not just local elements.
  //
  // @note
  // For single-rank execution (isDistributed() == false), this is a no-op.
  // For distributed execution, caller must synchronize:
  // - Mass matrix: once after computeFEInit()
  // - Force vectors: every time step after computeForces()
  template <typename ViewType>
  void synchronize(ViewType& field, const ParallelTopology& topo) {
    if (!topo.isDistributed()) {
      // Single rank: no synchronization needed
      return;
    }

    try {
      // Pack boundary values
      auto sendBufs = pack(field, topo);

      // Exchange with neighbors
      std::map<int, std::vector<float>> recvBufs;
      m_backend->exchange(sendBufs, recvBufs);

      // Accumulate (sum) received values
      accumulate(field, recvBufs, topo);
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

  // @brief Accumulate (sum) received values back into field at boundary nodes.
  //
  // Sums contributions from neighbor ranks at shared nodes. This assembly
  // operation is critical for correct boundary behavior in FEM:
  // - Boundary nodes see partial contributions from local elements
  // - Boundary nodes see partial contributions from neighbor elements
  // - Summing produces the complete value
  //
  // @tparam ViewType Type providing operator()(int) access
  // @param[in,out] field Field to accumulate into (modified in-place)
  // @param[in] recvBufs Received data from neighbors
  // @param[in] topo Topology describing shared nodes
  //
  // @throws std::length_error if received buffer size != expected node count
  // @throws std::runtime_error if expected neighbor data is missing
  //
  // @details
  // For each neighbor rank:
  // 1. Verify received buffer exists (error in distributed mode)
  // 2. Verify buffer size matches node count
  // 3. Sum received values into field at boundary nodes
  template <typename ViewType>
  static void accumulate(ViewType& field, const std::map<int, std::vector<float>>& recvBufs,
                         const ParallelTopology& topo) {
    for (const auto& [neighborRank, nodeIndices] : topo.sharedNodes) {
      auto it = recvBufs.find(neighborRank);

      if (it == recvBufs.end()) {
        // In distributed mode, missing data is an error
        // In serial mode (no actual neighbors), it's OK to skip
        if (topo.isDistributed()) {
          throw std::runtime_error("Expected data from rank " + std::to_string(neighborRank) +
                                   " but received nothing. Exchange failed or topology mismatch.");
        }
        continue;
      }

      const auto& buf = it->second;

      // Validate buffer size matches node count
      if (buf.size() != nodeIndices.size()) {
        throw std::length_error("Buffer size mismatch from rank " + std::to_string(neighborRank) + ": expected " +
                                std::to_string(nodeIndices.size()) + " values, got " + std::to_string(buf.size()));
      }

      // Accumulate (sum) received values
      for (size_t i = 0; i < nodeIndices.size(); ++i) {
        int nodeIdx = nodeIndices[i];
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

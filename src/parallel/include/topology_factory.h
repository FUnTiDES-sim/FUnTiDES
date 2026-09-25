#ifndef FUNTIDES_PARALLEL_INCLUDE_TOPOLOGY_FACTORY_H_
#define FUNTIDES_PARALLEL_INCLUDE_TOPOLOGY_FACTORY_H_
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "model.h"
#include "parallel_topology.h"

using namespace utils;

/**
 * @brief Tolerance used to decide whether a node lies on a partition boundary.
 */
struct TopologyTolerance {
  double absolute = 1e-6;    ///< Absolute tolerance on the x coordinate, same unit as the mesh coordinates
  bool auto_compute = true;  ///< If true, createFromMesh() overrides `absolute` with a fraction of the minimum spacing
};

/**
 * @brief Discovers the distributed topology of a mesh from its node coordinates.
 *
 * Nodes whose x coordinate matches the left or right edge of the local subdomain
 * (within a tolerance) are recorded as shared with the neighbor rank. Only a
 * decomposition along x, with rank r adjacent to ranks r-1 and r+1, is supported.
 *
 * @warning The tolerance is critical: too large and nodes near a boundary are
 * wrongly shared, too small and true boundary nodes are missed.
 */
class TopologyFactory {
 public:
  /**
   * @brief Builds the ParallelTopology of one rank by comparing node coordinates to the subdomain edges.
   *
   * For a serial run (size <= 1) the returned topology has no shared node. Otherwise
   * every expected neighbor must own at least one shared node.
   *
   * @tparam FloatType Floating point type of the mesh coordinates.
   * @tparam ScalarType Integer type of the mesh node indices.
   *
   * @param[in] mesh Mesh to inspect; its node indices are local to this rank.
   * @param[in] rank MPI rank of the caller, in [0, size).
   * @param[in] size Total number of ranks.
   * @param[in] origin_x X coordinate of the left edge of the local subdomain.
   * @param[in] domain_width_x Width of the local subdomain along x, must be positive.
   * @param[in] tol Boundary detection tolerance.
   *
   * @return Topology whose sharedNodes maps each neighbor rank to the local indices of the shared nodes.
   *
   * @throws std::invalid_argument if rank is outside [0, size) or domain_width_x <= 0.
   * @throws std::logic_error if a node lies on both edges, or if an expected neighbor has no shared node.
   *
   * @todo VERIFY: is tol.absolute expressed in mesh coordinate units (meters)?
   */
  template <typename FloatType, typename ScalarType>
  static ParallelTopology createFromMesh(const model::ModelApi<FloatType, ScalarType>& mesh, int rank, int size,
                                         FloatType origin_x, FloatType domain_width_x, TopologyTolerance tol = {}) {
    if (rank < 0 || rank >= size) {
      throw std::invalid_argument("Invalid rank " + std::to_string(rank) + " for numRanks " + std::to_string(size));
    }
    if (domain_width_x <= 0) {
      throw std::invalid_argument("Invalid domain_width_x: " + std::to_string(domain_width_x) + " (must be > 0)");
    }

    ParallelTopology topo;
    topo.myRank = rank;
    topo.numRanks = size;

    if (size <= 1) {
      return topo;
    }

    if (tol.auto_compute) {
      try {
        FloatType minDx = mesh.getMinSpacing();
        if (minDx > 0) {
          tol.absolute = minDx * 1e-4;
        }
      } catch (...) {
        tol.auto_compute = false;
      }
    }

    bool hasLeft = (rank > 0);
    bool hasRight = (rank < size - 1);

    FloatType left_x = origin_x;
    FloatType right_x = origin_x + domain_width_x;

    ScalarType numNodes = mesh.getNumberOfNodes();

    for (ScalarType i = 0; i < numNodes; ++i) {
      FloatType x = mesh.nodeCoord(i, 0);

      bool onLeft = hasLeft && (std::abs(x - left_x) < tol.absolute);
      bool onRight = hasRight && (std::abs(x - right_x) < tol.absolute);

      if (onLeft && onRight) {
        throw std::logic_error("Topology Error: Node " + std::to_string(i) +
                               " detected on both left and right boundaries.");
      }

      if (onLeft) {
        topo.sharedNodes[rank - 1].push_back(static_cast<int>(i));
      }
      if (onRight) {
        topo.sharedNodes[rank + 1].push_back(static_cast<int>(i));
      }
    }

    if (hasLeft && topo.sharedNodes[rank - 1].empty()) {
      throw std::logic_error("Topology Error: Rank " + std::to_string(rank) + " missing left boundary nodes.");
    }
    if (hasRight && topo.sharedNodes[rank + 1].empty()) {
      throw std::logic_error("Topology Error: Rank " + std::to_string(rank) + " missing right boundary nodes.");
    }

    return topo;
  }
};
#endif  // FUNTIDES_PARALLEL_INCLUDE_TOPOLOGY_FACTORY_H_

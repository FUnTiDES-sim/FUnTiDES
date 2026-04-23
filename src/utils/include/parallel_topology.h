#ifndef FUNTIDES_UTILS_INCLUDE_PARALLEL_TOPOLOGY_H_
#define FUNTIDES_UTILS_INCLUDE_PARALLEL_TOPOLOGY_H_
#include <map>
#include <vector>

namespace utils {

/**
 * @brief Describes the distributed connectivity of the mesh.
 *
 * Contains lists of local node indices shared with neighbor ranks.
 * Used to identify boundary nodes that require synchronization during
 * distributed execution.
 *
 * @details
 * For each neighbor rank, stores a vector of local node indices that
 * lie on the partition boundary. This information is discovered by
 * TopologyFactory and used by BoundarySynchronizer to coordinate
 * data exchange.
 */
struct ParallelTopology {
  int myRank = 0;    //< Current MPI rank (0-based)
  int numRanks = 1;  //< Total number of ranks

  // Map: Neighbor Rank ID -> List of LOCAL node indices shared with that
  // neighbor.
  std::map<int, std::vector<int>> sharedNodes;

  // @brief Returns true if mesh is distributed across multiple ranks.
  // @return true if numRanks > 1, false otherwise
  bool isDistributed() const { return numRanks > 1; }

  // @brief Returns the total number of boundary nodes across all neighbors.
  // @return Sum of all shared node counts
  size_t getTotalBoundaryNodes() const {
    size_t total = 0;
    for (const auto& [rank, nodes] : sharedNodes) {
      total += nodes.size();
    }
    return total;
  }
};

}  // namespace utils
#endif  // FUNTIDES_UTILS_INCLUDE_PARALLEL_TOPOLOGY_H_

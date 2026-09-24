#ifndef FUNTIDES_PARALLEL_INCLUDE_PARALLEL_TOPOLOGY_H_
#define FUNTIDES_PARALLEL_INCLUDE_PARALLEL_TOPOLOGY_H_
#include <map>
#include <vector>

namespace utils {

/**
 * @brief Distributed connectivity of the local mesh subdomain.
 *
 * For each neighbor rank, lists the local node indices lying on the partition
 * boundary shared with that rank. Used to decide which nodes need
 * synchronization during distributed execution.
 */
struct ParallelTopology {
  int myRank = 0;    ///< Rank of the current process (0-based).
  int numRanks = 1;  ///< Total number of ranks.

  /// Neighbor rank id -> local node indices shared with that neighbor.
  std::map<int, std::vector<int>> sharedNodes;

  /**
   * @brief Tells whether the mesh is split across several ranks.
   * @return true if numRanks > 1, false otherwise.
   */
  bool isDistributed() const { return numRanks > 1; }

  /**
   * @brief Counts the shared nodes over all neighbors.
   * @return Sum of the sizes of all lists in sharedNodes. A node shared with
   *         several neighbors is counted once per neighbor.
   */
  size_t getTotalBoundaryNodes() const {
    size_t total = 0;
    for (const auto& [rank, nodes] : sharedNodes) {
      total += nodes.size();
    }
    return total;
  }
};

}  // namespace utils
#endif  // FUNTIDES_PARALLEL_INCLUDE_PARALLEL_TOPOLOGY_H_

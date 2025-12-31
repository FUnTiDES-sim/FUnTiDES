#ifndef SRC_MODEL_MESH_API_INCLUDE_TOPOLOGY_FACTORY_H_
#define SRC_MODEL_MESH_API_INCLUDE_TOPOLOGY_FACTORY_H_

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "model.h"
#include "parallel_topology.h"

using namespace utils;

// @brief Configuration for boundary detection tolerance.
//
// Controls how TopologyFactory identifies boundary nodes by comparing
// node coordinates to expected boundary locations.
struct TopologyTolerance
{
  double absolute = 1e-6;    //< Absolute tolerance for coordinate comparison
  bool auto_compute = true;  //< Auto-compute from element spacing if true
};

// @brief Factory for discovering distributed topology from mesh coordinates.
//
// This class inspects a mesh and identifies which nodes lie on partition
// boundaries. It does this by comparing node coordinates to the expected
// boundary locations, rather than assuming index formulas. This approach
// is **correct by construction** - if the mesh structure changes internally,
// the topology discovery still works correctly.
//
// @details
// The topology discovery process:
// 1. Determines boundary locations from origin and domain width
// 2. Iterates all mesh nodes, comparing coordinates to boundaries
// 3. Validates that:
//    - No node appears on both left and right boundaries
//    - Expected boundaries have at least one node
// 4. Returns topology with shared node lists
//
// @warning
// Tolerance is critical. If too large, nodes far from boundaries may be
// included; if too small, actual boundary nodes may be missed. Use
// auto_compute=true for most cases, or provide a tolerance based on
// element size.
class TopologyFactory
{
 public:
  // @brief Generates topology by inspecting mesh coordinates.
  //
  // Scans all nodes in the mesh, comparing their coordinates to the
  // expected partition boundaries. Nodes within tolerance of a boundary
  // are marked as shared with the neighbor rank.
  //
  // @tparam FloatType Floating point type (float, double)
  // @tparam ScalarType Integer type for node counts (int, long)
  //
  // @param[in] mesh Mesh to inspect (must provide nodeCoord, getNumberOfNodes)
  // @param[in] rank Current MPI rank (0-based)
  // @param[in] size Total number of ranks
  // @param[in] origin_x Local subdomain origin in X direction
  // @param[in] domain_width_x Local subdomain width in X direction
  // @param[in] tol Tolerance configuration (optional, defaults provided)
  //
  // @return ParallelTopology with shared nodes identified
  //
  // @throws std::invalid_argument
  //   - if rank < 0 or rank >= size
  //   - if domain_width_x <= 0
  //
  // @throws std::logic_error
  //   - if node appears on both left and right boundaries
  //   - if expected left boundary found no nodes (rank > 0)
  //   - if expected right boundary found no nodes (rank < size-1)
  //
  // @note
  // For serial execution (size == 1), returns empty topology immediately.
  // For distributed execution, validates that all expected boundaries exist.
  template <typename FloatType, typename ScalarType>
  static ParallelTopology createFromMesh(
      const model::ModelApi<FloatType, ScalarType>& mesh, int rank, int size,
      FloatType origin_x, FloatType domain_width_x, TopologyTolerance tol = {})
  {
    // Input validation
    if (rank < 0 || rank >= size)
    {
      throw std::invalid_argument("Invalid rank " + std::to_string(rank) +
                                  " for numRanks " + std::to_string(size));
    }
    if (domain_width_x <= 0)
    {
      throw std::invalid_argument(
          "Invalid domain_width_x: " + std::to_string(domain_width_x) +
          " (must be > 0)");
    }

    ParallelTopology topo;
    topo.myRank = rank;
    topo.numRanks = size;

    if (size <= 1)
    {
      return topo;
    }

    // Auto-compute tolerance
    if (tol.auto_compute)
    {
      try
      {
        FloatType minDx = mesh.getMinSpacing();
        if (minDx > 0)
        {
          tol.absolute = minDx * 1e-4;
        }
      }
      catch (...)
      {
        tol.auto_compute = false;
      }
    }

    bool hasLeft = (rank > 0);
    bool hasRight = (rank < size - 1);

    FloatType left_x = origin_x;
    FloatType right_x = origin_x + domain_width_x;

    ScalarType numNodes = mesh.getNumberOfNodes();

    for (ScalarType i = 0; i < numNodes; ++i)
    {
      FloatType x = mesh.nodeCoord(i, 0);

      bool onLeft = hasLeft && (std::abs(x - left_x) < tol.absolute);
      bool onRight = hasRight && (std::abs(x - right_x) < tol.absolute);

      if (onLeft && onRight)
      {
        throw std::logic_error("Topology Error: Node " + std::to_string(i) +
                               " detected on both left and right boundaries.");
      }

      if (onLeft)
      {
        topo.sharedNodes[rank - 1].push_back(static_cast<int>(i));
      }
      if (onRight)
      {
        topo.sharedNodes[rank + 1].push_back(static_cast<int>(i));
      }
    }

    // Validation
    if (hasLeft && topo.sharedNodes[rank - 1].empty())
    {
      throw std::logic_error("Topology Error: Rank " + std::to_string(rank) +
                             " missing left boundary nodes.");
    }
    if (hasRight && topo.sharedNodes[rank + 1].empty())
    {
      throw std::logic_error("Topology Error: Rank " + std::to_string(rank) +
                             " missing right boundary nodes.");
    }

    return topo;
  }
};

#endif  // SRC_MODEL_MESH_API_INCLUDE_TOPOLOGY_FACTORY_H_

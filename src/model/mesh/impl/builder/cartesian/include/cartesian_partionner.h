#ifndef SRC_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_
#define SRC_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_

#include <algorithm>
#include <stdexcept>
#include <string>

#include "partitioning.h"

namespace model
{

// @brief 1D domain decomposition along X-axis.
//
// Distributes the global X-range evenly across ranks, with remainder
// elements distributed to the first few ranks. This ensures load balancing
// while maintaining consecutive element assignment.
//
// Example with 10 elements across 3 ranks:
// - Rank 0: elements 0-3 (4 elements)
// - Rank 1: elements 4-7 (4 elements)
// - Rank 2: elements 8-9 (2 elements)
//
// @tparam FloatType Floating point type (float, double)
// @tparam ScalarType Integer type for indices (int, long)
//
// @details
// For each rank, computes:
// - Local element count: base_ex + (1 if rank < remainder else 0)
// - Global element offset: rank * base_ex + min(rank, remainder)
// - Local physical size: local_ex * dx
// - Global origin: global_origin + offset * dx
//
// The global origin is critical for distributed execution, as it ensures
// TopologyFactory can identify boundary nodes by comparing coordinates.
template <typename FloatType, typename ScalarType>
class CartesianXPartitioner : public PartitioningStrategy<FloatType, ScalarType>
{
 public:
  // @brief Partition domain along X-axis.
  //
  // @param[in] global Global mesh parameters with global origin and dimensions
  // @param[in] rank Current MPI rank (0-based)
  // @param[in] size Total number of ranks
  //
  // @return Local parameters with:
  //   - ex: local element count
  //   - lx: local physical size
  //   - origin_x: global origin of this subdomain
  //   - y, z parameters: copied from global (unchanged)
  //
  // @throws std::invalid_argument if rank < 0 or rank >= size
  CartesianParams<FloatType, ScalarType> partition(
      const CartesianParams<FloatType, ScalarType>& global, int rank,
      int size) const override
  {
    if (rank < 0 || rank >= size)
    {
      throw std::invalid_argument("Invalid rank: " + std::to_string(rank) +
                                  " for size: " + std::to_string(size));
    }

    auto local = global;  // Copy global settings

    // Compute local element count
    ScalarType base_ex = global.ex / size;
    ScalarType remainder = global.ex % size;

    local.ex = base_ex + (rank < remainder ? 1 : 0);

    // Compute global element offset
    ScalarType element_offset_x =
        rank * base_ex + std::min(static_cast<ScalarType>(rank), remainder);

    // Element size
    FloatType dx = global.lx / global.ex;

    // Set local physical parameters
    local.lx = local.ex * dx;

    // CRITICAL: Set global origin for this subdomain
    // This ensures nodeCoord() returns absolute global coordinates
    local.origin_x = global.origin_x + element_offset_x * dx;
    local.origin_y = global.origin_y;
    local.origin_z = global.origin_z;

    return local;
  }
};

}  // namespace model

#endif  // SRC_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_

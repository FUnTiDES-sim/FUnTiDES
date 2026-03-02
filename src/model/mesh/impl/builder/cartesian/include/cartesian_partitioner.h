#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cartesian_params.h"
#include "partitioning.h"

namespace model
{

/**
 * @brief 1D domain decomposition along X-axis.
 *
 * Distributes the global X-range evenly across ranks, with remainder
 * elements distributed to the first few ranks. This ensures load balancing
 * while maintaining consecutive element assignment.
 *
 * Example with 10 elements across 3 ranks:
 * - Rank 0: elements 0-3 (4 elements)
 * - Rank 1: elements 4-7 (4 elements)
 * - Rank 2: elements 8-9 (2 elements)
 *
 * @tparam FloatType Floating point type (float, double)
 * @tparam ScalarType Integer type for indices (int, long)
 *
 * @details
 * For each rank, computes:
 * - Local element count: base_ex + (1 if rank < remainder else 0)
 * - Global element offset: rank * base_ex + min(rank, remainder)
 * - Local physical size: local_ex * dx
 * - Global origin: global_origin + offset * dx
 *
 * The global origin is critical for distributed execution, as it ensures
 * TopologyFactory can identify boundary nodes by comparing coordinates.
 */
template <typename FloatType, typename ScalarType>
class CartesianXPartitioner
    : public PartitioningStrategy<CartesianParams<FloatType, ScalarType>>
{
 public:
  using Params = CartesianParams<FloatType, ScalarType>;

  /**
   * @brief Implements the partition method for Cartesian parameters.
   * Matches the signature: LocalParams partition(const GlobalParams&, int, int)
   * const
   */
  Params partition(const Params& global, int rank, int size) const override
  {
    // Validation
    if (size <= 0)
    {
      throw std::invalid_argument("CartesianPartitioner: size must be > 0");
    }
    if (rank < 0 || rank >= size)
    {
      throw std::invalid_argument(
          "CartesianPartitioner: rank must be between 0 and size-1");
    }

    auto local = global;

    // 1D Decomposition along X axis
    ScalarType base_ex = global.ex / size;
    ScalarType remainder = global.ex % size;

    // Determine local element count
    local.ex = base_ex + (rank < remainder ? 1 : 0);

    // Calculate Global Offset (in elements)
    ScalarType element_offset_x =
        rank * base_ex + std::min((ScalarType)rank, remainder);

    // Calculate Element Size
    FloatType dx = global.lx / global.ex;

    // Set Local Physical Size
    local.lx = local.ex * dx;

    // Critical: Set Global Origin for this subdomain
    local.origin_x = global.origin_x + element_offset_x * dx;
    local.origin_y = global.origin_y;
    local.origin_z = global.origin_z;

    // Store global dimensions for reference
    local.global_lx = global.lx;
    local.global_ly = global.ly;
    local.global_lz = global.lz;

    // Store global origins
    local.global_origin_x = global.origin_x;
    local.global_origin_y = global.origin_y;
    local.global_origin_z = global.origin_z;

    return local;
  }
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_

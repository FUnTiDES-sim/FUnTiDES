#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "cartesian_params.h"
#include "partitioning.h"

namespace model {

/**
 * @brief 1D domain decomposition of a Cartesian box along the X axis.
 *
 * The global element count along X is split evenly across ranks; the first
 * (ex % size) ranks receive one extra element, so each rank owns a contiguous
 * block of elements. Y and Z are not split.
 *
 * Example with 10 elements across 3 ranks: rank 0 owns elements 0-3, rank 1
 * owns 4-7, rank 2 owns 8-9.
 *
 * Each local Params carries the local element count and size along X, the
 * origin of the subdomain, and the global sizes and origin. The global origin
 * lets the topology factory identify boundary nodes by comparing coordinates
 * in distributed runs.
 *
 * @tparam FloatType Floating point type of the coordinates and sizes.
 * @tparam ScalarType Integer type of the element counts.
 */
template <typename FloatType, typename ScalarType>
class CartesianXPartitioner : public PartitioningStrategy<CartesianParams<FloatType, ScalarType>> {
 public:
  using Params = CartesianParams<FloatType, ScalarType>;

  /**
   * @brief Builds the parameters of the subdomain owned by one rank.
   *
   * All fields of @p global not listed in the class description are copied unchanged.
   *
   * @param[in] global Parameters of the whole domain; global.ex must be at least 1.
   * @param[in] rank Rank of the subdomain, in [0, size).
   * @param[in] size Number of ranks, strictly positive.
   * @return Parameters of the subdomain of @p rank.
   * @throws std::invalid_argument If size <= 0 or rank is outside [0, size).
   */
  Params partition(const Params& global, int rank, int size) const override {
    if (size <= 0) {
      throw std::invalid_argument("CartesianPartitioner: size must be > 0");
    }
    if (rank < 0 || rank >= size) {
      throw std::invalid_argument("CartesianPartitioner: rank must be between 0 and size-1");
    }

    auto local = global;

    ScalarType base_ex = global.ex / size;
    ScalarType remainder = global.ex % size;

    local.ex = base_ex + (rank < remainder ? 1 : 0);

    // Offset of the first local element, in elements.
    ScalarType element_offset_x = rank * base_ex + std::min((ScalarType)rank, remainder);

    FloatType dx = global.lx / global.ex;

    local.lx = local.ex * dx;

    local.origin_x = global.origin_x + element_offset_x * dx;
    local.origin_y = global.origin_y;
    local.origin_z = global.origin_z;

    local.global_lx = global.lx;
    local.global_ly = global.ly;
    local.global_lz = global.lz;

    local.global_origin_x = global.origin_x;
    local.global_origin_y = global.origin_y;
    local.global_origin_z = global.origin_z;

    return local;
  }
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARTITIONER_H_

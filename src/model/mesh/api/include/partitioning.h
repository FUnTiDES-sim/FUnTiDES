#ifndef SRC_MODEL_MESH_API_INCLUDE_PARTITIONING_H_
#define SRC_MODEL_MESH_API_INCLUDE_PARTITIONING_H_

#include "cartesian_params.h"

namespace model
{

// @brief Interface for domain decomposition strategy.
//
// Responsible for determining local mesh parameters from global ones.
// Different partitioning strategies (1D X-decomposition, 2D XY-decomposition,
// METIS-based, etc.) can be implemented by subclassing this interface.
//
// @tparam FloatType Floating point type (float, double)
// @tparam ScalarType Integer type for indices (int, long)
template <typename FloatType, typename ScalarType>
class PartitioningStrategy
{
 public:
  virtual ~PartitioningStrategy() = default;

  // @brief Partition global domain into local subdomain for given rank.
  //
  // Takes global mesh parameters and computes the local parameters for
  // the given rank. This includes:
  // - Number of elements in local domain
  // - Physical size of local domain
  // - Global origin of local domain (for coordinate mapping)
  //
  // @param[in] globalParams Global mesh parameters
  // @param[in] rank Current MPI rank (0-based)
  // @param[in] numRanks Total number of ranks
  //
  // @return CartesianParams containing local mesh parameters with correct
  //         origin, element count, and dimensions
  //
  // @throws std::invalid_argument if rank or numRanks are invalid
  virtual CartesianParams<FloatType, ScalarType> partition(
      const CartesianParams<FloatType, ScalarType>& globalParams, int rank,
      int numRanks) const = 0;
};

}  // namespace model

#endif  // SRC_MODEL_MESH_API_INCLUDE_PARTITIONING_H_

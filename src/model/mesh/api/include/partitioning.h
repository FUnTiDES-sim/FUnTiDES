#ifndef FUNTIDES_MODEL_MESH_API_INCLUDE_PARTITIONING_H_
#define FUNTIDES_MODEL_MESH_API_INCLUDE_PARTITIONING_H_
#include "cartesian_params.h"

namespace model {

/** @brief Interface for domain decomposition strategy.
 * Responsible for determining local mesh parameters from global ones.
 * Different partitioning strategies (1D X-decomposition, 2D XY-decomposition,
 * METIS-based, etc.) can be implemented by subclassing this interface.
 *
 * @tparam GlobalParams Type describing the global problem (e.g.,
 * CartesianParams, or a filename string)
 * @tparam LocalParams Type describing the local partition (e.g.,
 * CartesianParams, or a list of element IDs)
 */
template <typename GlobalParams, typename LocalParams = GlobalParams>
class PartitioningStrategy {
 public:
  virtual ~PartitioningStrategy() = default;

  /** @brief Partition global domain into local subdomain for given rank.
   *
   * Takes global mesh parameters and computes the local parameters for
   * the given rank. This includes:
   * - Number of elements in local domain
   * - Physical size of local domain
   * - Global origin of local domain (for coordinate mapping)
   *
   * @param[in] globalParams Description of the global domain/mesh
   * @param[in] rank Current MPI rank (0-based)
   * @param[in] numRanks Total number of ranks
   *
   * @return CartesianParams containing local mesh parameters with correct
   *         origin, element count, and dimensions
   *
   * @throws std::invalid_argument if rank or
   * numRanks are invalid
   */
  virtual LocalParams partition(const GlobalParams& globalParams, int rank, int numRanks) const = 0;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_PARTITIONING_H_

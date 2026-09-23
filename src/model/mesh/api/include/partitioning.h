#ifndef FUNTIDES_MODEL_MESH_API_INCLUDE_PARTITIONING_H_
#define FUNTIDES_MODEL_MESH_API_INCLUDE_PARTITIONING_H_
#include "cartesian_params.h"

namespace model {

/**
 * @brief Domain decomposition strategy: derives the parameters of one rank's subdomain from
 * the parameters of the global problem.
 *
 * @tparam GlobalParams Description of the global problem.
 * @tparam LocalParams Description of one subdomain.
 */
template <typename GlobalParams, typename LocalParams = GlobalParams>
class PartitioningStrategy {
 public:
  virtual ~PartitioningStrategy() = default;

  /**
   * @brief Compute the subdomain owned by a rank.
   *
   * The result gives the local element count, the local physical size and the origin of the
   * subdomain in global coordinates.
   *
   * @param[in] globalParams Description of the global domain.
   * @param[in] rank Rank whose subdomain is requested, in [0, numRanks).
   * @param[in] numRanks Total number of ranks.
   * @return Parameters of the subdomain of rank.
   * @throws std::invalid_argument if numRanks <= 0 or rank is outside [0, numRanks).
   */
  virtual LocalParams partition(const GlobalParams& globalParams, int rank, int numRanks) const = 0;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_PARTITIONING_H_

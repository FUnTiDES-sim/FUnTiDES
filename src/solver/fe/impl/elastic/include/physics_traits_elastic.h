#ifndef FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_
#include "physics_traits.h"
#include "rhs_elastic.h"
#include "wavefield_elastic.h"

namespace solver
{
namespace fe
{

/**
 * @brief Specialization for elastic physics.
 *
 * Elastic wave propagation uses three displacement components (ux, uy, uz).
 */
template <>
struct PhysicsTraits<utils::enums::physicType::kElastic>
{
  /// Human-readable name for logging
  static constexpr const char* kName = "Elastic";

  /// Concrete types for device access
  using WavefieldType = WavefieldElastic;
  using RhsType = RhsElastic;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_

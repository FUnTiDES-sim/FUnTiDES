#ifndef FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_
#include "physics_traits.h"
#include "rhs_elastic.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {

/**
 * @brief Elastic specialization of PhysicsTraits.
 *
 * Associates the elastic wavefield and source types (three displacement
 * components: ux, uy, uz) with the elastic physics tag.
 */
template <>
struct PhysicsTraits<utils::enums::physicType::kElastic> {
  /// Human-readable name for logging.
  static constexpr const char* kName = "Elastic";

  /// Wavefield type holding the elastic fields.
  using WavefieldType = WavefieldElastic;
  /// Source term type for elastic physics.
  using RhsType = RhsElastic;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_

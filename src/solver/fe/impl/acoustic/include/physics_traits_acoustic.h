#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_

#include "physics_traits.h"
#include "rhs_acoustic.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Traits of the acoustic physics (single scalar pressure field).
 *
 * Maps the acoustic physic type to its wavefield and source (right-hand side)
 * types.
 */
template <>
struct PhysicsTraits<utils::enums::physicType::kAcoustic> {
  /// Human-readable name, for logging.
  static constexpr const char* kName = "Acoustic";

  using WavefieldType = WavefieldAcoustic;  ///< Wavefield storage type.
  using RhsType = RhsAcoustic;              ///< Source term type.
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_

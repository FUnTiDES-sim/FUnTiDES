#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_

#include "sem_enums.h"

namespace solver {
namespace fe {

/**
 * @brief Maps a physics type to the concrete wavefield and RHS types that
 * device kernels use instead of the Wavefield and Rhs base classes.
 *
 * The primary template is empty; each supported physics provides a
 * specialization defining:
 * - kName: name of the physics, for logging;
 * - WavefieldType: concrete Wavefield type;
 * - RhsType: concrete Rhs type.
 *
 * @tparam PHYSICS Physics type.
 * @see docs/design.md, "Device calls on mesh objects".
 */
template <utils::enums::physicType PHYSICS>
struct PhysicsTraits {
  static constexpr const char* kName = "";  ///< Physics name, for logging.
  using WavefieldType = void;  ///< Concrete Wavefield type.
  using RhsType = void;  ///< Concrete Rhs type.
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_
#ifndef FUNTIDES_GRADIENT_API_INCLUDE_PHYSICS_TRAITS_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_PHYSICS_TRAITS_H_

#include "sem_enums.h"

namespace gradient {

/**
 * @brief Maps a physics type to the concrete wavefield view and gradient types
 * used by the differentiator data of that physics.
 *
 * The primary template is empty; each supported physics provides a
 * specialization defining:
 * - kName: name of the physics, for logging;
 * - WavefieldViewForwardType: concrete forward WavefieldView type;
 * - WavefieldViewBackwardType: concrete adjoint WavefieldView type;
 * - GradientType: concrete Gradient type.
 *
 * @tparam PHYSICS Physics type.
 * @see docs/design.md, "Device calls on mesh objects".
 */
template <utils::enums::physicType PHYSICS>
struct PhysicsTraits {
  static constexpr const char* kName = "";  ///< Physics name, for logging.
  using WavefieldViewForwardType = void;  ///< Concrete forward view type.
  using WavefieldViewBackwardType = void;  ///< Concrete adjoint view type.
  using GradientType = void;  ///< Concrete Gradient type.
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_PHYSICS_TRAITS_H_

#ifndef FUNTIDES_GRADIENT_API_INCLUDE_PHYSICS_TRAITS_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_PHYSICS_TRAITS_H_

#include "sem_enums.h"

namespace gradient
{

/**
 * @brief Compile-time properties for each gradient physics type.
 *
 * All specializations MUST define:
 * - kName:                     Human-readable name for logging (const char*)
 * - WavefieldViewForwardType:  Read-only forward wavefield view for this physics
 * - WavefieldViewBackwardType: Read-only adjoint wavefield view for this physics
 * - GradientType:              Concrete gradient data type for this physics
 *
 * @tparam PHYSICS The physics type (kAcoustic, kElastic)
 *
 * Include the physics-specific trait header for full specialization:
 *   - physics_traits_acoustic.h
 *   - physics_traits_elastic.h
 */
template <utils::enums::physicType PHYSICS>
struct PhysicsTraits
{
  static constexpr const char* kName = "";
  using WavefieldViewForwardType = void;
  using WavefieldViewBackwardType = void;
  using GradientType = void;
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_PHYSICS_TRAITS_H_

#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_

#include "gradient_acoustic.h"
#include "physics_traits.h"
#include "wavefield_view_backward_acoustic.h"
#include "wavefield_view_forward_acoustic.h"

namespace gradient {

/**
 * @brief Maps the acoustic physics tag to its wavefield view and gradient types.
 */
template <>
struct PhysicsTraits<utils::enums::physicType::kAcoustic> {
  static constexpr const char* kName = "Acoustic";  ///< Human-readable physics name.
  using WavefieldViewForwardType = WavefieldViewForwardAcoustic;  ///< Forward wavefield view.
  using WavefieldViewBackwardType = WavefieldViewBackwardAcoustic;  ///< Backward (adjoint) wavefield view.
  using GradientType = GradientAcoustic;  ///< Gradient container.
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_

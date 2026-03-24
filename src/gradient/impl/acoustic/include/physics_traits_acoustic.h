#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_

#include "physics_traits.h"
#include "wavefield_view_forward_acoustic.h"
#include "wavefield_view_backward_acoustic.h"
#include "gradient_acoustic.h"

namespace gradient
{

template <>
struct PhysicsTraits<utils::enums::physicType::kAcoustic>
{
  static constexpr const char* kName = "Acoustic";
  using WavefieldViewForwardType  = WavefieldViewForwardAcoustic;
  using WavefieldViewBackwardType = WavefieldViewBackwardAcoustic;
  using GradientType              = GradientAcoustic;
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_

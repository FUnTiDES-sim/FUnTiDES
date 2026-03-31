#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_

#include "gradient_elastic.h"
#include "physics_traits.h"
#include "wavefield_view_backward_elastic.h"
#include "wavefield_view_forward_elastic.h"

namespace gradient
{

template <>
struct PhysicsTraits<utils::enums::physicType::kElastic>
{
  static constexpr const char* kName = "Elastic";
  using WavefieldViewForwardType = WavefieldViewForwardElastic;
  using WavefieldViewBackwardType = WavefieldViewBackwardElastic;
  using GradientType = GradientElastic;
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_

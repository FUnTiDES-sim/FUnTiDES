#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_

#include "gradient_elastic.h"
#include "physics_traits.h"
#include "wavefield_view_backward_elastic.h"
#include "wavefield_view_forward_elastic.h"

namespace gradient {

/**
 * @brief Type bundle of the elastic physics: forward and backward wavefield
 * views and gradient container.
 */
template <>
struct PhysicsTraits<utils::enums::physicType::kElastic> {
  static constexpr const char* kName = "Elastic";                  ///< Human-readable physics name.
  using WavefieldViewForwardType = WavefieldViewForwardElastic;    ///< Forward wavefield view.
  using WavefieldViewBackwardType = WavefieldViewBackwardElastic;  ///< Backward wavefield view.
  using GradientType = GradientElastic;                            ///< Gradient container.
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_PHYSICS_TRAITS_ELASTIC_H_

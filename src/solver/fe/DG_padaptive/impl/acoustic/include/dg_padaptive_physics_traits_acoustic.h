#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_PHYSICS_TRAITS_ACOUSTIC_H_

#include "dg_padaptive_rhs_acoustic.h"
#include "dg_padaptive_wavefield_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Specialization for p-adaptive acoustic physics.
 *
 * Acoustic wave propagation uses a single scalar pressure field for each order.
 */
struct DGPAdaptivePhysicsTraits {
  /// Human-readable name for logging
  static constexpr const char* kName = "DGPAdaptiveAcoustic";

  /// Concrete types for device access
  using WavefieldType = DGPAdaptiveWavefieldAcoustic;
  using RhsType = DGPAdaptiveRhsAcoustic;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_PHYSICS_TRAITS_ACOUSTIC_H_

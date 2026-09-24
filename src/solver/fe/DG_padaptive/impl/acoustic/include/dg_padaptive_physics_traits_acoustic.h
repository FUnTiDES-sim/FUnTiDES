#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_PHYSICS_TRAITS_ACOUSTIC_H_

#include "dg_padaptive_rhs_acoustic.h"
#include "dg_padaptive_wavefield_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Bundles the wavefield and source types of the p-adaptive acoustic solver.
 *
 * Only carries type aliases and a name; it holds no data.
 */
struct DGPAdaptivePhysicsTraits {
  /// Name used in log messages.
  static constexpr const char* kName = "DGPAdaptiveAcoustic";

  using WavefieldType = DGPAdaptiveWavefieldAcoustic;  ///< Wavefield storage type.
  using RhsType = DGPAdaptiveRhsAcoustic;              ///< Source term storage type.
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_PHYSICS_TRAITS_ACOUSTIC_H_

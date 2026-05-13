#ifndef FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_PHYSICS_TRAITS_ACOUSTIC_H_

#include "dg_wavefield_acoustic.h"
#include "rhs_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Specialization for acoustic physics.
 *
 * Acoustic wave propagation uses a single scalar pressure field.
 */
struct DGPhysicsTraits {
  /// Human-readable name for logging
  static constexpr const char* kName = "DGAcoustic";

  /// Concrete types for device access
  using WavefieldType = DGWavefieldAcoustic;
  using RhsType = RhsAcoustic;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_PHYSICS_TRAITS_ACOUSTIC_H_

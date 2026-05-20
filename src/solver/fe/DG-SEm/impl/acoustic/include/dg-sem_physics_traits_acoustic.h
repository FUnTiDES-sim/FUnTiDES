#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_PHYSICS_TRAITS_ACOUSTIC_H_

#include "dg-sem_wavefield_acoustic.h"
#include "dg-sem_rhs_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Specialization for acoustic physics.
 *
 * Acoustic wave propagation uses a single scalar pressure field.
 */
struct DGSEMPhysicsTraits {
  /// Human-readable name for logging
  static constexpr const char* kName = "DGSEMAcoustic";

  /// Concrete types for device access
  using WavefieldType = DGSEMWavefieldAcoustic;
  using RhsType = DGSEMRhsAcoustic;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_PHYSICS_TRAITS_ACOUSTIC_H_

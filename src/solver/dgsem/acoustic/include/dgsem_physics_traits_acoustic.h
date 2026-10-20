#ifndef FUNTIDES_SOLVER_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_PHYSICS_TRAITS_ACOUSTIC_H_

#include "dgsem_rhs_acoustic.h"
#include "dgsem_wavefield_acoustic.h"

namespace solver {

/**
 * @brief Type bundle for the acoustic DG-SEM solver (single scalar pressure field).
 */
struct DGSEMPhysicsTraits {
  static constexpr const char* kName = "DGSEMAcoustic";  ///< Name used in log messages.

  using WavefieldType = DGSEMWavefieldAcoustic;  ///< Wavefield container type.
  using RhsType = DGSEMRhsAcoustic;              ///< Right-hand-side (source) type.
};

}  // namespace solver
#endif  // FUNTIDES_SOLVER_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_PHYSICS_TRAITS_ACOUSTIC_H_

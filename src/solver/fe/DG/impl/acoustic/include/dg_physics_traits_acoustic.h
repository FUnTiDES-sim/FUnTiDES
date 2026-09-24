#ifndef FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_PHYSICS_TRAITS_ACOUSTIC_H_

#include "dg_wavefield_acoustic.h"
#include "rhs_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Compile-time bundle of the wavefield and source types of acoustic DG.
 *
 * Acoustic propagation uses a single scalar pressure field.
 */
struct DGPhysicsTraits {
  static constexpr const char* kName = "DGAcoustic";  ///< Name used in log messages.

  using WavefieldType = DGWavefieldAcoustic;  ///< Wavefield storage type.
  using RhsType = RhsAcoustic;                ///< Source (right-hand side) type.
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_IMPL_ACOUSTIC_INCLUDE_DG_PHYSICS_TRAITS_ACOUSTIC_H_

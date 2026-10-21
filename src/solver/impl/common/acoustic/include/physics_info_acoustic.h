#ifndef FUNTIDES_SOLVER_IMPL_COMMON_ACOUSTIC_INCLUDE_PHYSICS_INFO_ACOUSTIC_H_
#define FUNTIDES_SOLVER_IMPL_COMMON_ACOUSTIC_INCLUDE_PHYSICS_INFO_ACOUSTIC_H_

#include "physics_info.h"

namespace solver {

/**
 * @brief Description of the acoustic physics.
 *
 * Acoustic propagation carries a single scalar pressure field and a single
 * source component.
 */
template <>
struct PhysicsInfo<utils::enums::physicType::kAcoustic> {
  /// Human-readable name, for logging.
  static constexpr const char* kName = "Acoustic";

  static constexpr int kNumFields = 1;  ///< Pressure components.
  static constexpr int kNumRhs = 1;    ///< Source components per RHS.
};

}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_COMMON_ACOUSTIC_INCLUDE_PHYSICS_INFO_ACOUSTIC_H_

#ifndef FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_PHYSICS_INFO_H_
#define FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_PHYSICS_INFO_H_

#include "sem_enums.h"

namespace solver {

/**
 * @brief Method-agnostic description of a physics.
 *
 * Carries the quantities that every discretization method needs from its
 * physics and that are intrinsic to it: the number of wavefield components
 * and the number of source components per right-hand side. Solver methods
 * that only need these counts (e.g. the DG variants) consume PhysicsInfo
 * instead of a method's PhysicsTraits, which bind concrete storage types.
 *
 * The primary template is empty; each supported physics provides a
 * specialization in physics_info_<physics>.h.
 *
 * @tparam PHYSICS Physics type.
 */
template <utils::enums::physicType PHYSICS>
struct PhysicsInfo {
  static constexpr const char* kName = "";  ///< Physics name, for logging.
  static constexpr int kNumFields = 0;      ///< Wavefield components.
  static constexpr int kNumRhs = 0;         ///< Source components per RHS.
};

}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_PHYSICS_INFO_H_

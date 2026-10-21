#ifndef FUNTIDES_SOLVER_IMPL_COMMON_ELASTIC_INCLUDE_PHYSICS_INFO_ELASTIC_H_
#define FUNTIDES_SOLVER_IMPL_COMMON_ELASTIC_INCLUDE_PHYSICS_INFO_ELASTIC_H_

#include "physics_info.h"

namespace solver {

/**
 * @brief Description of the elastic physics.
 *
 * Elastic propagation is formulated in displacement: three wavefield
 * components (u_x, u_y, u_z) and three source components per RHS.
 */
template <>
struct PhysicsInfo<utils::enums::physicType::kElastic> {
  /// Human-readable name, for logging.
  static constexpr const char* kName = "Elastic";

  static constexpr int kNumFields = 3;  ///< Displacement components.
  static constexpr int kNumRhs = 3;    ///< Source components per RHS.
};

}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_COMMON_ELASTIC_INCLUDE_PHYSICS_INFO_ELASTIC_H_

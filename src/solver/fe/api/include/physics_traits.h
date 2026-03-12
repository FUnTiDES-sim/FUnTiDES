#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_
#include "sem_enums.h"

namespace solver
{
namespace fe
{

/**
 * @brief Compile-time properties for each physics type.
 *
 * All specializations MUST define:
 * - kName: Human-readable name for logging (const char*)
 * - WavefieldType: Concrete wavefield type for device access
 * - RhsType: Concrete RHS type for device access
 *
 * @tparam PHYSICS The physics type (kAcoustic, kElastic, etc.)
 *
 * Note: physicType enum is forward declared. Include sem_enums.h
 * or physics-specific trait headers for full definitions.
 */
template <enums::physicType PHYSICS>
struct PhysicsTraits
{
  static constexpr const char* kName = "";
  using WavefieldType = void;
  using RhsType = void;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_
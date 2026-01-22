#ifndef SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_
#define SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_

namespace solver
{
namespace fe
{

namespace enums
{
// Forward declaration - full definition in sem_enums.h
enum class physicType : int;
}  // namespace enums

/**
 * @brief Compile-time properties for each physics type.
 *
 * All specializations MUST define:
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
  using WavefieldType = void;
  using RhsType = void;

  static_assert(sizeof(PHYSICS) == 0,
                "PhysicsTraits must be specialized for this physics type");
};

}  // namespace fe
}  // namespace solver
#endif  // SOLVER_FE_API_INCLUDE_PHYSICS_TRAITS_H_
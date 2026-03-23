#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_
#include "physics_traits.h"
#include "rhs_acoustic.h"
#include "wavefield_acoustic.h"
#include "gradients_acoustic.h"

namespace solver
{
namespace fe
{

/**
 * @brief Specialization for acoustic physics.
 *
 * Acoustic wave propagation uses a single scalar pressure field.
 */
template <>
struct PhysicsTraits<enums::physicType::kAcoustic>
{
  /// Human-readable name for logging
  static constexpr const char* kName = "Acoustic";

  /// Concrete types for device access
  using WavefieldType = WavefieldAcoustic;
  using RhsType = RhsAcoustic;
};

/**
 * @brief Specialization for acoustic physics.
 *
 * Acoustic wave propagation with gradients uses a single scalar pressure field.
 */
template <>
struct PhysicsTraits<enums::physicType::kAcoustic>
{
  /// Human-readable name for logging
  static constexpr const char* kName = "Acoustic";

  /// Concrete types for device access
  using WavefieldType = WavefieldAcoustic;
  using RhsType = RhsAcoustic;
  using GradientsType = GradientsAcoustic;
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_PHYSICS_TRAITS_ACOUSTIC_H_

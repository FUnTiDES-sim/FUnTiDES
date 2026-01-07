#pragma once

#include "rhs_acoustic.h"
#include "rhs_elastic.h"
#include "sem_enums.h"
#include "sem_solver.h"
#include "solver_base.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace solver
{
namespace fe
{

using physicType = solver::fe::enums::physicType;

//============================================================================
// Physics Traits - Compile-time properties for each physics type
//============================================================================

/**
 * @brief Primary template for physics traits.
 * @tparam PHYSICS The physics type (kAcoustic or kElastic)
 */
template <physicType PHYSICS>
struct PhysicsTraits;

/**
 * @brief Specialization for acoustic physics.
 *
 * Acoustic wave propagation uses a single scalar pressure field.
 */
template <>
struct PhysicsTraits<enums::physicType::kAcoustic>
{
  /// Number of solution fields (1 for pressure)
  static constexpr int kNumFields = 1;

  /// Number of RHS (source) components
  static constexpr int kNumRhsComponents = 1;

  /// Human-readable name for logging
  static constexpr const char* kName = "Acoustic";

  /// Primary field name
  static constexpr const char* kFieldNames[1] = {"pressure"};

  /// Concrete types for device access
  using WavefieldType = WavefieldAcoustic;
  using RhsType = RhsAcoustic;
};

/**
 * @brief Specialization for elastic physics.
 *
 * Elastic wave propagation uses three displacement components (ux, uy, uz).
 */
template <>
struct PhysicsTraits<enums::physicType::kElastic>
{
  /// Number of solution fields (3 for displacement vector)
  static constexpr int kNumFields = 3;

  /// Number of RHS (source) components
  static constexpr int kNumRhsComponents = 3;

  /// Human-readable name for logging
  static constexpr const char* kName = "Elastic";

  /// Field names for each component
  static constexpr const char* kFieldNames[3] = {"ux", "uy", "uz"};

  /// Concrete types for device access
  using WavefieldType = WavefieldElastic;
  using RhsType = RhsElastic;
};

//============================================================================
// Unified Data Structure
//============================================================================

/**
 * @brief Unified data structure for SEM solver.
 *
 * Holds solution fields, RHS terms, and time indices for time-stepping.
 * The number of fields and RHS components is determined at compile time
 * based on the physics type.
 *
 * @tparam PHYSICS The physics type (kAcoustic or kElastic)
 */
template <physicType PHYSICS>
struct SEMsolverData : public SolverBase::DataStruct
{
  using Traits = PhysicsTraits<PHYSICS>;
  static constexpr int kNumFields = Traits::kNumFields;
  static constexpr int kNumRhs = Traits::kNumRhsComponents;

  // Use concrete types from PhysicsTraits to avoid virtual dispatch on device
  using WavefieldType = typename Traits::WavefieldType;
  using RhsType = typename Traits::RhsType;

  /**
   * @brief Constructor for acoustic physics (single field).
   */
  template <physicType P = PHYSICS,
            typename = std::enable_if_t<P == enums::physicType::kAcoustic>>
  SEMsolverData(const WavefieldAcoustic& wavefield,
                const RhsAcoustic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs)
  {
  }

  /**
   * @brief Constructor for elastic physics (three fields).
   */
  template <physicType P = PHYSICS,
            typename = std::enable_if_t<P == enums::physicType::kElastic>>
  SEMsolverData(const WavefieldElastic& wavefield,
                const RhsElastic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs)
  {
  }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  ARRAY_REAL_VIEW getRhsTerm(int i) const { return m_rhs.getTerm(i); }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  VECTOR_INT_VIEW getRhsElement() const { return m_rhs.getElement(); }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  ARRAY_REAL_VIEW getRhsWeights() const { return m_rhs.getWeights(); }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  VECTOR_REAL_VIEW getCurrentField(int i) const
  {
    return m_wavefield.getCurrentField(i);
  }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  VECTOR_REAL_VIEW getPreviousField(int i) const
  {
    return m_wavefield.getPreviousField(i);
  }

  void print() const override
  {
    std::cout << "SEMsolverData<" << Traits::kName << ">" << std::endl;
    for (int f = 0; f < kNumFields; ++f)
    {
      std::cout << "Field[" << f << "] (" << Traits::kFieldNames[f]
                << ") size: " << getCurrentField(f).extent(0) << std::endl;
    }
    for (int r = 0; r < kNumRhs; ++r)
    {
      std::cout << "RHS[" << r << "] size: " << getRhsTerm(r).extent(0)
                << std::endl;
    }
    std::cout << "RHS Element size: " << getRhsElement().extent(0) << std::endl;
    std::cout << "RHS Weights size: " << getRhsWeights().extent(0) << std::endl;
  }

  WavefieldType m_wavefield;  ///< Wavefield stored by value for GPU (lightweight view handles)
  RhsType m_rhs;              ///< RHS stored by value for GPU (lightweight view handles)
};

//============================================================================
// Backward Compatibility Type Aliases for Data Structures
//============================================================================

using SEMsolverDataAcoustic = SEMsolverData<enums::physicType::kAcoustic>;
using SEMsolverDataElastic = SEMsolverData<enums::physicType::kElastic>;

}  // namespace fe
}  // namespace solver

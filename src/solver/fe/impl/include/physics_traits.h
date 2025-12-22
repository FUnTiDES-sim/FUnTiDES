#pragma once

#include "sem_enums.h"
#include "sem_solver.h"
#include "solver_base.h"

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

  /**
   * @brief Constructor for acoustic physics (single field).
   */
  template <physicType P = PHYSICS,
            typename = std::enable_if_t<P == enums::physicType::kAcoustic>>
  SEMsolverData(int i1, int i2, ARRAY_REAL_VIEW rhsTerm,
                ARRAY_REAL_VIEW fieldGlobal, VECTOR_INT_VIEW rhsElement,
                ARRAY_REAL_VIEW rhsWeights)
      : m_i1(i1), m_i2(i2), m_rhsElement(rhsElement), m_rhsWeights(rhsWeights)
  {
    m_rhsTerms[0] = rhsTerm;
    m_fieldsGlobal[0] = fieldGlobal;
  }

  /**
   * @brief Constructor for elastic physics (three fields).
   */
  template <physicType P = PHYSICS,
            typename = std::enable_if_t<P == enums::physicType::kElastic>>
  SEMsolverData(int i1, int i2, ARRAY_REAL_VIEW rhsTermx,
                ARRAY_REAL_VIEW rhsTermy, ARRAY_REAL_VIEW rhsTermz,
                ARRAY_REAL_VIEW uxnGlobal, ARRAY_REAL_VIEW uynGlobal,
                ARRAY_REAL_VIEW uznGlobal, VECTOR_INT_VIEW rhsElement,
                ARRAY_REAL_VIEW rhsWeights)
      : m_i1(i1), m_i2(i2), m_rhsElement(rhsElement), m_rhsWeights(rhsWeights)
  {
    m_rhsTerms[0] = rhsTermx;
    m_rhsTerms[1] = rhsTermy;
    m_rhsTerms[2] = rhsTermz;
    m_fieldsGlobal[0] = uxnGlobal;
    m_fieldsGlobal[1] = uynGlobal;
    m_fieldsGlobal[2] = uznGlobal;
  }

  void print() const override
  {
    std::cout << "SEMsolverData<" << Traits::kName << ">: i1=" << m_i1
              << ", i2=" << m_i2 << std::endl;
    for (int f = 0; f < kNumFields; ++f)
    {
      std::cout << "Field[" << f << "] (" << Traits::kFieldNames[f]
                << ") size: " << m_fieldsGlobal[f].extent(0) << std::endl;
    }
    for (int r = 0; r < kNumRhs; ++r)
    {
      std::cout << "RHS[" << r << "] size: " << m_rhsTerms[r].extent(0)
                << std::endl;
    }
    std::cout << "RHS Element size: " << m_rhsElement.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_rhsWeights.extent(0) << std::endl;
  }

  int m_i1;  ///< Previous time step index
  int m_i2;  ///< Current time step index

  std::array<ARRAY_REAL_VIEW, kNumRhs> m_rhsTerms;  ///< RHS forcing terms
  std::array<ARRAY_REAL_VIEW, kNumFields> m_fieldsGlobal;  ///< Solution fields

  VECTOR_INT_VIEW m_rhsElement;  ///< Source element indices
  ARRAY_REAL_VIEW m_rhsWeights;  ///< Forcing weights per node
};

//============================================================================
// Backward Compatibility Type Aliases for Data Structures
//============================================================================

using SEMsolverDataAcoustic = SEMsolverData<enums::physicType::kAcoustic>;
using SEMsolverDataElastic = SEMsolverData<enums::physicType::kElastic>;

}  // namespace fe
}  // namespace solver

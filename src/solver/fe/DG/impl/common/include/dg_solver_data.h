#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_DATA_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_DATA_H_
#include <iostream>

#include "physics_traits_acoustic.h"
#include "solver.h"

namespace solver {
namespace fe {

using physicType = utils::enums::physicType;

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
struct DGsolverData : public Solver::DataStruct {
  using Traits = PhysicsTraits<PHYSICS>;
  static constexpr int kNumFields = Traits::WavefieldType::kNumFields;
  static constexpr int kNumRhs = Traits::RhsType::kNumRhsComponents;

  // Use concrete types from PhysicsTraits to avoid virtual dispatch on device
  using WavefieldType = typename Traits::WavefieldType;
  using RhsType = typename Traits::RhsType;

  /**
   * @brief Constructor for acoustic physics (single field).
   */
  template <physicType P = PHYSICS, typename = std::enable_if_t<P == physicType::kAcoustic>>
  DGsolverData(const WavefieldAcoustic& wavefield, const RhsAcoustic& rhs) : m_wavefield(wavefield), m_rhs(rhs) {}

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getRhsTerm(int i) const { return m_rhs.getTerm(i); }

  PROXY_HOST_DEVICE
  VECTOR_INT_VIEW getRhsElement() const { return m_rhs.getElement(); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getRhsWeights() const { return m_rhs.getWeights(); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getCurrentField(int i) const { return m_wavefield.getCurrentField(i); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getPreviousField(int i) const { return m_wavefield.getPreviousField(i); }

  void swapWavefields() { m_wavefield.swap(); }

  void print() const override {
    std::cout << "DGsolverData<" << Traits::kName << ">" << std::endl;
    for (int f = 0; f < kNumFields; ++f) {
      std::cout << "Field[" << f << "] (" << Traits::WavefieldType::kFieldNames[f]
                << ") size: " << getCurrentField(f).extent(0) << std::endl;
    }
    for (int r = 0; r < kNumRhs; ++r) {
      std::cout << "RHS[" << r << "] size: " << getRhsTerm(r).extent(0) << std::endl;
    }
    std::cout << "RHS Element size: " << getRhsElement().extent(0) << std::endl;
    std::cout << "RHS Weights size: " << getRhsWeights().extent(0) << std::endl;
  }

  bool isDistributed{false};
  WavefieldType m_wavefield;  ///< Wavefield stored by value for GPU
                              ///< (lightweight view handles)
  RhsType m_rhs;              ///< RHS stored by value for GPU (lightweight view handles)
};

//============================================================================
// Backward Compatibility Type Aliases for Data Structures
//============================================================================

using DGsolverDataAcoustic = SEMsolverData<physicType::kAcoustic>;

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_DATA_H_

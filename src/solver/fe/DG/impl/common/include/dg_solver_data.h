#ifndef FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_
#define FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_

#include <iostream>

#include "data_type.h"
#include "dg_physics_traits_acoustic.h"
#include "solver.h"

namespace solver {
namespace fe {

/**
 * @brief Data structure for the DG acoustic solver.
 *
 * Holds solution fields (per-element 2-D) and RHS terms (same layout as SEM).
 * Fields are indexed (n_elem, n_dof_per_elem); RHS follows SEM conventions.
 */
struct DGsolverDataAcoustic : public Solver::DataStruct {
  // Use concrete types from DGPhysicsTraits to avoid virtual dispatch on device
  using WavefieldType = typename DGPhysicsTraits::WavefieldType;
  using RhsType = typename DGPhysicsTraits::RhsType;

  WavefieldType m_wavefield;  ///< Wavefield stored by value for GPU (lightweight view handles)
  RhsType m_rhs;              ///< RHS stored by value for GPU (lightweight view handles)

  bool isDistributed{false};

  DGsolverDataAcoustic(const DGWavefieldAcoustic& wavefield, const RhsAcoustic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getCurrentField(int i) const { return m_wavefield.getCurrentField(i); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getPreviousField(int i) const { return m_wavefield.getPreviousField(i); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getRhsTerm(int i) const { return m_rhs.getTerm(i); }

  PROXY_HOST_DEVICE
  VECTOR_INT_VIEW getRhsElement() const { return m_rhs.getElement(); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getRhsWeights() const { return m_rhs.getWeights(); }

  void swapWavefields() { m_wavefield.swap(); }

  void print() const override {
    std::cout << "DGsolverDataAcoustic: " << m_wavefield.getPreviousField(0).extent(0) << " elems x "
              << m_wavefield.getPreviousField(0).extent(1) << " dofs" << std::endl;
  }
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_SOLVER_DATA_H_

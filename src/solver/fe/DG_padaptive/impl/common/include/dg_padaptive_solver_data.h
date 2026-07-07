#ifndef FUNTIDES_SRC_SOLVER_FE_DG_PADAPTIVE__IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_DATA_H_
#define FUNTIDES_SRC_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_DATA_H_

#include <iostream>

#include "data_type.h"
#include "dg_padaptive_rhs_acoustic.h"
#include "dg_padaptive_wavefield_acoustic.h"
#include "solver.h"

namespace solver {
namespace fe {

/**
 * @brief Data structure passed to DGPAdaptiveSolver at each time step.
 *
 * Combines two DG acoustic wavefield, one for each approximation (Lagrange basis) order,
 * and the acoustic source term.
 */
struct DGPAdaptiveSolverData : public Solver::DataStruct {
  /**
   * @param wavefield Combined the two DG wavefield with different approximation order.
   * @param rhs       source term.
   */
  DGPAdaptiveSolverData(const DGPAdaptiveWavefieldAcoustic& wavefield, const DGPAdaptiveRhsAcoustic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  void print() const override {
    m_wavefield.print();
    m_rhs.print();
  }

  /// Swap previous/current wavefields (call once per time step after
  /// computeOneStep).
  void swapWavefields() { m_wavefield.swap(); }

  DGPAdaptiveWavefieldAcoustic m_wavefield;  ///< Combined wavefield
  DGPAdaptiveRhsAcoustic m_rhs;              ///< source

  bool isDistributed{false};
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SRC_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_DATA_H_

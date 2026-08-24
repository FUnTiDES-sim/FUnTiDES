#ifndef FUNTIDES_SRC_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_DATA_H_
#define FUNTIDES_SRC_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_DATA_H_

#include "dg-sem_rhs_acoustic.h"
#include "dg-sem_wavefield_acoustic.h"
#include "solver.h"

namespace solver {
namespace fe {

/**
 * @brief Data structure passed to DGSEMsolver at each time step.
 *
 * Combines the acoustic wavefield of the DG domain, the acoustic wavefield of the SEM domain,
 * and the acoustic source term.
 */
struct DGSEMsolverData : public Solver::DataStruct {
  /**
   * @param wavefield Combined DG-SEM wavefield.
   * @param rhs       source term (either in DG or SEM).
   */
  DGSEMsolverData(const DGSEMWavefieldAcoustic& wavefield, const DGSEMRhsAcoustic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  void print() const override {
    m_wavefield.print();
    m_rhs.print();
  }

  /// Swap previous/current wavefields (call once per time step after
  /// computeOneStep).
  void swapWavefields() { m_wavefield.swap(); }

  DGSEMWavefieldAcoustic m_wavefield;  ///< Combined wavefield DG+SEM
  DGSEMRhsAcoustic m_rhs;              ///< source

  bool isDistributed{false};
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SRC_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_DATA_H_

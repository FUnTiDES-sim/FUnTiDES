#ifndef FUNTIDES_SOLVER_IMPL_DGSEM_COMMON_INCLUDE_DGSEM_SOLVER_DATA_H_
#define FUNTIDES_SOLVER_IMPL_DGSEM_COMMON_INCLUDE_DGSEM_SOLVER_DATA_H_

#include "dgsem_rhs_acoustic.h"
#include "dgsem_wavefield_acoustic.h"
#include "solver.h"

namespace solver {

/**
 * @brief Per-time-step data of the coupled DG-SEM acoustic solver.
 *
 * Bundles the acoustic wavefields of the DG and SEM domains with the acoustic source term.
 */
struct DGSEMsolverData : public Solver::DataStruct {
  /**
   * @brief Builds the data from a combined wavefield and a source term.
   * @param[in] wavefield Combined DG-SEM wavefield (copied).
   * @param[in] rhs       Source term for the DG and SEM domains (copied).
   */
  DGSEMsolverData(const DGSEMWavefieldAcoustic& wavefield, const DGSEMRhsAcoustic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  /// @brief Prints the wavefield and the source term.
  void print() const override {
    m_wavefield.print();
    m_rhs.print();
  }

  /// @brief Swaps the previous and current wavefields; call once per time step, after the step.
  void swapWavefields() { m_wavefield.swap(); }

  DGSEMWavefieldAcoustic m_wavefield;  ///< Combined DG and SEM wavefield.
  DGSEMRhsAcoustic m_rhs;              ///< Source term.

  /// @todo VERIFY: what does isDistributed select (MPI run with several ranks?) and who reads it?
  bool isDistributed{false};
};

}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_DGSEM_COMMON_INCLUDE_DGSEM_SOLVER_DATA_H_

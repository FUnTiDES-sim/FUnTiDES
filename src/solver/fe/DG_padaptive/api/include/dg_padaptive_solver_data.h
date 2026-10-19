#ifndef FUNTIDES_SRC_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_DATA_H_
#define FUNTIDES_SRC_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_DATA_H_

#include <iostream>

#include "data_type.h"
#include "dg_padaptive_rhs_acoustic.h"
#include "dg_padaptive_wavefield_acoustic.h"
#include "solver.h"

namespace solver {
namespace fe {

/**
 * @brief Per-time-step data of the p-adaptive DG acoustic solver.
 *
 * Bundles the two-order acoustic wavefield (one DG field per approximation
 * order) and the acoustic source term. Members are copied shallowly (Kokkos
 * views).
 */
struct DGPAdaptiveSolverData : public Solver::DataStruct {
  /**
   * @brief Builds the data object from a wavefield and a source term.
   * @param[in] wavefield Wavefield holding the fields of both approximation orders.
   * @param[in] rhs       Acoustic source term.
   */
  DGPAdaptiveSolverData(const DGPAdaptiveWavefieldAcoustic& wavefield, const DGPAdaptiveRhsAcoustic& rhs)
      : m_wavefield(wavefield), m_rhs(rhs) {}

  /// @brief Prints the wavefield and the source term.
  void print() const override {
    m_wavefield.print();
    m_rhs.print();
  }

  /// @brief Swaps the previous and current wavefields; call once per time step after computeOneStep.
  void swapWavefields() { m_wavefield.swap(); }

  DGPAdaptiveWavefieldAcoustic m_wavefield;  ///< Wavefield of both approximation orders.
  DGPAdaptiveRhsAcoustic m_rhs;              ///< Acoustic source term.

  bool isDistributed{false};  ///< @todo VERIFY: what does this flag mean and who reads it?
};

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SRC_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_DATA_H_

#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_WAVEFIELD_ACOUSTIC_H_

#include <data_type.h>

#include "dg_wavefield_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Combined wavefield for the p-adative DG solver.
 *
 * Holds one acoustic pressure field in each DG sub-domain at two consecutive time levels (previous and current).
 * This struct is passed to DGPAdaptiveSolverData at each time step.
 */
struct DGPAdaptiveWavefieldAcoustic {
  /// Total number of solution fields: 2 dg acoustic (p).
  static constexpr int kNumFields = 2;

  /// Field names in order: pMinDGp, pMaxDGp
  static constexpr const char* kFieldNames[2] = {"pMinDGpressure", "pMaxDGpressure"};

  DGPAdaptiveWavefieldAcoustic(arrayReal pnPMinDGPrev, arrayReal pnPMinDGCurr, arrayReal pnPMaxDGPrev, arrayReal pnPMaxDGCurr)
      : m_pMinAcoustic(pnPMinDGPrev, pnPMinDGCurr), m_pMaxAcoustic(pnPMaxDGPrev, pnPMaxDGCurr) {}

  int getNumFields() const { return kNumFields; }

  const char* const* getFieldNames() const { return kFieldNames; }

  /**
   * @brief Get the current field of order pMin.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMinCurrentField(int i) const { return m_pMinAcoustic.getCurrentField(0); }

  /**
   * @brief Get the current field of order pMax.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMaxCurrentField(int i) const { return m_pMaxAcoustic.getCurrentField(0); }

  /**
   * @brief Get the previous field of order pMin.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMinPreviousField(int i) const { return m_pMinAcoustic.getPreviousField(0); }

  /**
   * @brief Get the previous field of order pMax.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMaxPreviousField(int i) const { return m_pMaxAcoustic.getPreviousField(0); }

  void swap() {
    m_pMinAcoustic.swap();
    m_pMaxAcoustic.swap();
  }

  void print() const {
    m_pMinAcoustic.print();
    m_pMaxAcoustic.print();
  }

  DGWavefieldAcoustic m_pMinAcoustic;  ///< Acoustic pressure wavefield for order pMin
  DGWavefieldAcoustic m_pMaxAcoustic;  ///< Acoustic pressure wavefield for order pMax
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_WAVEFIELD_ACOUSTIC_H_

#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_WAVEFIELD_ACOUSTIC_H_

#include <data_type.h>

#include "dg_wavefield_acoustic.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Combined wavefield for the dg-sem coupled solver.
 *
 * Holds both the acoustic pressure field of dg and sem at two consecutive time levels (previous and current).
 * This struct is passed to DGSEMsolverDataAcoustic at each time step.
 */
struct DGSEMWavefieldAcoustic {
  /// Total number of solution fields: 1 dg acoustic (p) + 1 sem acoustic (p).
  static constexpr int kNumFields = 2;

  /// Field names in order: DGp, SEMp
  static constexpr const char* kFieldNames[2] = {"DGpressure", "SEMpressure"};

  DGSEMWavefieldAcoustic(arrayReal pnDGPrev, arrayReal pnDGCurr, vectorReal pnSEMPrev, vectorReal pnSEMCurr)
      : m_DGacoustic(pnDGPrev, pnDGCurr), m_SEMacoustic(pnSEMPrev, pnSEMCurr) {}

  int getNumFields() const { return kNumFields; }

  const char* const* getFieldNames() const { return kFieldNames; }

  /**
   * @brief Get the current field of DG.
   */
  PROXY_HOST_DEVICE
  arrayReal getDGCurrentField(int i) const { return m_DGacoustic.getCurrentField(0); }

  /**
   * @brief Get the current field of SEM.
   */
  PROXY_HOST_DEVICE
  vectorReal getSEMCurrentField(int i) const { return m_SEMacoustic.getCurrentField(0); }

  /**
   * @brief Get the previous field of DG.
   */
  PROXY_HOST_DEVICE
  arrayReal getDGPreviousField(int i) const { return m_DGacoustic.getPreviousField(0); }

  /**
   * @brief Get the previous field of SEM.
   */
  PROXY_HOST_DEVICE
  vectorReal getSEMPreviousField(int i) const { return m_SEMacoustic.getPreviousField(0); }

  void swap() {
    m_DGacoustic.swap();
    m_SEMacoustic.swap();
  }

  void print() const {
    m_DGacoustic.print();
    m_SEMacoustic.print();
  }

  DGWavefieldAcoustic m_DGacoustic;  ///< Acoustic pressure wavefield for DG
  WavefieldAcoustic m_SEMacoustic;   ///< Acoustic pressure wavefield for SEM
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_WAVEFIELD_ACOUSTIC_H_

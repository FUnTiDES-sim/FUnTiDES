#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_WAVEFIELD_ACOUSTIC_H_

#include <data_type.h>

#include "dg_wavefield_acoustic.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Pair of acoustic pressure wavefields, one for the DG domain and one for the SEM domain.
 *
 * Each sub-wavefield holds the previous and current time levels of its pressure field.
 * The two sub-wavefields are public members, accessed directly or through the getters below.
 */
struct DGSEMWavefieldAcoustic {
  /// Total number of solution fields: 1 DG pressure + 1 SEM pressure.
  static constexpr int kNumFields = 2;

  /// Field names, in order: DG pressure, SEM pressure.
  static constexpr const char* kFieldNames[2] = {"DGpressure", "SEMpressure"};

  /**
   * @brief Builds the two sub-wavefields from their previous and current pressure arrays.
   * @param[in] pnDGPrev DG pressure at the previous time level.
   * @param[in] pnDGCurr DG pressure at the current time level.
   * @param[in] pnSEMPrev SEM pressure at the previous time level.
   * @param[in] pnSEMCurr SEM pressure at the current time level.
   */
  DGSEMWavefieldAcoustic(arrayReal pnDGPrev, arrayReal pnDGCurr, vectorReal pnSEMPrev, vectorReal pnSEMCurr)
      : m_DGacoustic(pnDGPrev, pnDGCurr), m_SEMacoustic(pnSEMPrev, pnSEMCurr) {}

  /// @return Number of solution fields (kNumFields).
  int getNumFields() const { return kNumFields; }

  /// @return Array of kNumFields field names, in the order of kFieldNames.
  const char* const* getFieldNames() const { return kFieldNames; }

  /**
   * @brief Returns the DG pressure at the current time level.
   * @param[in] i Field index, currently ignored: the DG pressure is always returned.
   */
  PROXY_HOST_DEVICE
  arrayReal getDGCurrentField(int i) const { return m_DGacoustic.getCurrentField(0); }

  /**
   * @brief Returns the SEM pressure at the current time level.
   * @param[in] i Field index, currently ignored: the SEM pressure is always returned.
   */
  PROXY_HOST_DEVICE
  vectorReal getSEMCurrentField(int i) const { return m_SEMacoustic.getCurrentField(0); }

  /**
   * @brief Returns the DG pressure at the previous time level.
   * @param[in] i Field index, currently ignored: the DG pressure is always returned.
   */
  PROXY_HOST_DEVICE
  arrayReal getDGPreviousField(int i) const { return m_DGacoustic.getPreviousField(0); }

  /**
   * @brief Returns the SEM pressure at the previous time level.
   * @param[in] i Field index, currently ignored: the SEM pressure is always returned.
   */
  PROXY_HOST_DEVICE
  vectorReal getSEMPreviousField(int i) const { return m_SEMacoustic.getPreviousField(0); }

  /// @brief Swaps the previous and current time levels of both sub-wavefields.
  void swap() {
    m_DGacoustic.swap();
    m_SEMacoustic.swap();
  }

  /// @brief Prints both sub-wavefields.
  void print() const {
    m_DGacoustic.print();
    m_SEMacoustic.print();
  }

  DGWavefieldAcoustic m_DGacoustic;  ///< Acoustic pressure wavefield of the DG domain.
  WavefieldAcoustic m_SEMacoustic;   ///< Acoustic pressure wavefield of the SEM domain.
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_WAVEFIELD_ACOUSTIC_H_

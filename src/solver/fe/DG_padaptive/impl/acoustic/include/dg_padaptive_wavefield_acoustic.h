#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_WAVEFIELD_ACOUSTIC_H_

#include <data_type.h>

#include "dg_wavefield_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Acoustic pressure wavefields of the two p-adaptive DG sub-domains (orders pMin and pMax).
 *
 * Each sub-domain holds a previous and a current time level, wrapped in a DGWavefieldAcoustic.
 * The two members are public and are accessed directly by the sub-solvers.
 */
struct DGPAdaptiveWavefieldAcoustic {
  /// Number of solution fields (one pressure field per sub-domain).
  static constexpr int kNumFields = 2;

  /// Field names, in the order pMin, pMax.
  static constexpr const char* kFieldNames[2] = {"pMinDGpressure", "pMaxDGpressure"};

  /**
   * @brief Wraps the four pressure arrays.
   * @param[in] pnPMinDGPrev Pressure of the pMin sub-domain at the previous time level.
   * @param[in] pnPMinDGCurr Pressure of the pMin sub-domain at the current time level.
   * @param[in] pnPMaxDGPrev Pressure of the pMax sub-domain at the previous time level.
   * @param[in] pnPMaxDGCurr Pressure of the pMax sub-domain at the current time level.
   */
  DGPAdaptiveWavefieldAcoustic(arrayReal pnPMinDGPrev, arrayReal pnPMinDGCurr, arrayReal pnPMaxDGPrev,
                               arrayReal pnPMaxDGCurr)
      : m_pMinAcoustic(pnPMinDGPrev, pnPMinDGCurr), m_pMaxAcoustic(pnPMaxDGPrev, pnPMaxDGCurr) {}

  /// @return Number of solution fields.
  int getNumFields() const { return kNumFields; }

  /// @return Array of kNumFields field names.
  const char* const* getFieldNames() const { return kFieldNames; }

  /**
   * @brief Current pressure of the pMin sub-domain.
   * @param[in] i Unused.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMinCurrentField(int i) const { return m_pMinAcoustic.getCurrentField(0); }

  /**
   * @brief Current pressure of the pMax sub-domain.
   * @param[in] i Unused.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMaxCurrentField(int i) const { return m_pMaxAcoustic.getCurrentField(0); }

  /**
   * @brief Previous pressure of the pMin sub-domain.
   * @param[in] i Unused.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMinPreviousField(int i) const { return m_pMinAcoustic.getPreviousField(0); }

  /**
   * @brief Previous pressure of the pMax sub-domain.
   * @param[in] i Unused.
   */
  PROXY_HOST_DEVICE
  arrayReal getPMaxPreviousField(int i) const { return m_pMaxAcoustic.getPreviousField(0); }

  /// @brief Swaps the previous and current levels of both sub-domains.
  void swap() {
    m_pMinAcoustic.swap();
    m_pMaxAcoustic.swap();
  }

  /// @brief Prints both sub-domain wavefields.
  void print() const {
    m_pMinAcoustic.print();
    m_pMaxAcoustic.print();
  }

  DGWavefieldAcoustic m_pMinAcoustic;  ///< Pressure wavefield of the pMin sub-domain.
  DGWavefieldAcoustic m_pMaxAcoustic;  ///< Pressure wavefield of the pMax sub-domain.
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_WAVEFIELD_ACOUSTIC_H_

#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_RHS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_RHS_ACOUSTIC_H_

#include <data_type.h>

#include "rhs.h"
#include "rhs_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief RHS data structure for the p-adaptive DG solver.
 *
 * Holds one acoustic source in each DG sub-domain. 
 * Either source may be zero-initialised when the
 * corresponding domain is inactive. The sub-solvers receive their respective
 * @ref RhsAcoustic member directly, so @p getTerm() is
 * provided for interface compliance only.
 */
struct DGPAdaptiveRhsAcoustic : public Rhs {
  /// Number of RHS components: 2 dg acoustic (p).
  static constexpr int kNumRhsComponents = 2;

  /**
   * @param pMin_acoustic_term  2D array of acoustic source signals (n_src x n_t) for the pMin approximation order domain.
   * @param pMax_acoustic_term  2D array of acoustic source signals (n_src x n_t) for the pMax approximation order domain.
   * @param element             Indices of elements containing source points.
   * @param pMin_weights        Per-node pMin weights for source distribution.
   * @param pMax_weights        Per-node pMax weights for source distribution.
   */
  DGPAdaptiveRhsAcoustic(arrayReal pMin_acoustic_term, arrayReal pMax_acoustic_term, vectorInt element, arrayReal pMin_weights, arrayReal pMax_weights)
      : m_rhs_pMinAcoustic(pMin_acoustic_term, element, pMin_weights), m_rhs_pMaxAcoustic(pMax_acoustic_term, element, pMax_weights) {}

  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  /// @brief Returns term i: 0 = pMin, 1 = pMax.
  PROXY_HOST_DEVICE
  arrayReal getTerm(int i) const override {
    if (i == 0) return m_rhs_pMinAcoustic.getTerm(0);
    return m_rhs_pMaxAcoustic.getTerm(0);
  }

  PROXY_HOST_DEVICE
  vectorInt getElement() const { return m_rhs_pMinAcoustic.getElement(); }

  PROXY_HOST_DEVICE
  arrayReal getWeights(int i) const { 
    if (i == 0) return m_rhs_pMinAcoustic.getWeights();
    return m_rhs_pMaxAcoustic.getWeights();
  }
  PROXY_HOST_DEVICE
  arrayReal getWeights() const override {
    Kokkos::abort("getWeights need an order indicator (0 : order_min_, 1 : order_max_)");
    return {};
  }

  void print() const override {
    m_rhs_pMinAcoustic.print();
    m_rhs_pMaxAcoustic.print();
  }

  RhsAcoustic m_rhs_pMinAcoustic;   ///< pMin domain source
  RhsAcoustic m_rhs_pMaxAcoustic;   ///< pMax domain source
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_RHS_ACOUSTIC_H_

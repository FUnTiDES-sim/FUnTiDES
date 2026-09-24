#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_RHS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_RHS_ACOUSTIC_H_

#include <data_type.h>

#include "rhs.h"
#include "rhs_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Acoustic source terms of the two approximation-order domains (pMin and pMax).
 *
 * Holds one acoustic source per domain, built from the same element indices.
 * Each sub-solver reads its own source through the public members; the
 * indexed accessors are the generic interface.
 */
struct DGPAdaptiveRhsAcoustic : public Rhs {
  /// Number of RHS components (one per domain: pMin, pMax).
  static constexpr int kNumRhsComponents = 2;

  /**
   * @param pMin_acoustic_term  Source signals for the pMin domain.
   * @param pMax_acoustic_term  Source signals for the pMax domain.
   * @param element             Indices of the elements containing source points.
   * @param pMin_weights        Per-node source distribution weights, pMin domain.
   * @param pMax_weights        Per-node source distribution weights, pMax domain.
   * @todo VERIFY: exact shapes of the term and weight arrays (n_src x n_t for terms?) and units.
   */
  DGPAdaptiveRhsAcoustic(arrayReal pMin_acoustic_term, arrayReal pMax_acoustic_term, vectorInt element,
                         arrayReal pMin_weights, arrayReal pMax_weights)
      : m_rhs_pMinAcoustic(pMin_acoustic_term, element, pMin_weights),
        m_rhs_pMaxAcoustic(pMax_acoustic_term, element, pMax_weights) {}

  /// @brief Returns the number of RHS components.
  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  /**
   * @brief Returns the source signals of one domain.
   * @param[in] i  0 for pMin, any other value for pMax.
   */
  PROXY_HOST_DEVICE
  arrayReal getTerm(int i) const override {
    if (i == 0) return m_rhs_pMinAcoustic.getTerm(0);
    return m_rhs_pMaxAcoustic.getTerm(0);
  }

  /// @brief Returns the indices of the source elements (shared by both domains).
  PROXY_HOST_DEVICE
  vectorInt getElement() const { return m_rhs_pMinAcoustic.getElement(); }

  /**
   * @brief Returns the source weights of one domain.
   * @param[in] i  0 for pMin, any other value for pMax.
   */
  PROXY_HOST_DEVICE
  arrayReal getWeights(int i) const {
    if (i == 0) return m_rhs_pMinAcoustic.getWeights();
    return m_rhs_pMaxAcoustic.getWeights();
  }

  /**
   * @brief Not usable: aborts, because the weights depend on the domain.
   * @see getWeights(int)
   */
  PROXY_HOST_DEVICE
  arrayReal getWeights() const override {
    Kokkos::abort("getWeights need an order indicator (0 : order_min_, 1 : order_max_)");
    return {};
  }

  /// @brief Prints both sources.
  void print() const override {
    m_rhs_pMinAcoustic.print();
    m_rhs_pMaxAcoustic.print();
  }

  RhsAcoustic m_rhs_pMinAcoustic;  ///< Source of the pMin domain.
  RhsAcoustic m_rhs_pMaxAcoustic;  ///< Source of the pMax domain.
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_ACOUSTIC_INCLUDE_DG_PADAPTIVE_RHS_ACOUSTIC_H_

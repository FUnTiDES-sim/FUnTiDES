#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_RHS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_RHS_ACOUSTIC_H_

#include <data_type.h>

#include "rhs.h"
#include "rhs_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief RHS data structure for the dg-sem coupled solver.
 *
 * Holds one acoustic source in DG domain and one acoustic source
 * in SEM domain. Either source may be zero-initialised when the
 * corresponding domain is inactive. The sub-solvers receive their respective
 * @ref RhsAcoustic member directly, so @p getTerm() is
 * provided for interface compliance only.
 */
struct DGSEMRhsAcoustic : public Rhs {
  /// Number of RHS components: 1 dg acoustic (p) + 1 sem acoustic (p).s
  static constexpr int kNumRhsComponents = 2;

  /**
   * @param dg_acoustic_term  2D array of acoustic source signals (n_src x n_t) for the DG domain.
   * @param sem_acoustic_term  2D array of acoustic source signals (n_src x n_t).
   * @param element        Indices of elements containing source points.
   * @param weights        Per-node weights for source distribution.
   */
  DGSEMRhsAcoustic(arrayReal dg_acoustic_term, arrayReal sem_acoustic_term, vectorInt element, arrayReal weights)
      : m_rhs_DGacoustic(dg_acoustic_term, element, weights), m_rhs_SEMacoustic(sem_acoustic_term, element, weights) {}

  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  /// @brief Returns term i: 0 = DG, 1 = SEM.
  PROXY_HOST_DEVICE
  arrayReal getTerm(int i) const override {
    if (i == 0) return m_rhs_DGacoustic.getTerm(0);
    return m_rhs_SEMacoustic.getTerm(0);
  }

  PROXY_HOST_DEVICE
  vectorInt getElement() const { return m_rhs_DGacoustic.getElement(); }

  PROXY_HOST_DEVICE
  arrayReal getWeights() const { return m_rhs_DGacoustic.getWeights(); }

  void print() const override {
    m_rhs_DGacoustic.print();
    m_rhs_SEMacoustic.print();
  }

  RhsAcoustic m_rhs_DGacoustic;   ///< source from DG
  RhsAcoustic m_rhs_SEMacoustic;  ///< source from SEM
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_RHS_ACOUSTIC_H_

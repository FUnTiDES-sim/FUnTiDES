#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_RHS_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_RHS_ACOUSTIC_H_

#include <data_type.h>

#include "rhs.h"
#include "rhs_acoustic.h"

namespace solver {
namespace fe {

/**
 * @brief Right-hand side of a coupled DG / SEM acoustic problem.
 *
 * Holds two acoustic sources, one for the DG domain and one for the SEM domain.
 * Both share the same source elements and weights. Either source may be
 * zero-initialised when the corresponding domain is inactive. Each domain
 * reads its own source through the public members; getTerm() exposes both
 * through the generic Rhs interface.
 */
struct DGSEMRhsAcoustic : public Rhs {
  /// Number of RHS components: one DG acoustic term and one SEM acoustic term.
  static constexpr int kNumRhsComponents = 2;

  /**
   * @param dg_acoustic_term   Source signals of the DG domain, size n_src x n_t.
   * @param sem_acoustic_term  Source signals of the SEM domain, size n_src x n_t.
   * @param element            Indices of the elements containing the source points.
   * @param weights            Per-node weights distributing each source over its element.
   * @todo VERIFY: layout of weights (rows = sources, columns = element nodes?).
   */
  DGSEMRhsAcoustic(arrayReal dg_acoustic_term, arrayReal sem_acoustic_term, vectorInt element, arrayReal weights)
      : m_rhs_DGacoustic(dg_acoustic_term, element, weights), m_rhs_SEMacoustic(sem_acoustic_term, element, weights) {}

  /// @return Number of RHS components, always kNumRhsComponents.
  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  /**
   * @brief Returns the source signals of one domain.
   * @param i  0 selects the DG source, any other value selects the SEM source.
   * @return Source signals, size n_src x n_t.
   */
  PROXY_HOST_DEVICE
  arrayReal getTerm(int i) const override {
    if (i == 0) return m_rhs_DGacoustic.getTerm(0);
    return m_rhs_SEMacoustic.getTerm(0);
  }

  /// @return Indices of the elements containing the source points.
  PROXY_HOST_DEVICE
  vectorInt getElement() const { return m_rhs_DGacoustic.getElement(); }

  /// @return Per-node source weights, see the constructor.
  PROXY_HOST_DEVICE
  arrayReal getWeights() const { return m_rhs_DGacoustic.getWeights(); }

  /// @brief Prints both sources.
  void print() const override {
    m_rhs_DGacoustic.print();
    m_rhs_SEMacoustic.print();
  }

  RhsAcoustic m_rhs_DGacoustic;   ///< Source of the DG domain.
  RhsAcoustic m_rhs_SEMacoustic;  ///< Source of the SEM domain.
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_ACOUSTIC_INCLUDE_DG_SEM_RHS_ACOUSTIC_H_

#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_RHS_ACOUSTOELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_RHS_ACOUSTOELASTIC_H_

#include <data_type.h>

#include "rhs.h"
#include "rhs_acoustic.h"

namespace solver
{
namespace fe
{

/**
 * @brief RHS data structure for the acousto-elastic coupled solver.
 *
 * In the V1 bicouche configuration the source is located in the acoustic
 * (fluid) domain only. The elastic domain has no external forcing.
 */
struct RhsAcoustoElastic : public Rhs
{
  /// Number of RHS components (1: acoustic pressure source)
  static constexpr int kNumRhsComponents = 1;

  /**
   * @param term   2D array of source time signals (n_sources x n_samples).
   * @param element  Indices of elements containing source points.
   * @param weights  Per-node weights for source distribution within elements.
   */
  RhsAcoustoElastic(ARRAY_REAL_VIEW term, VECTOR_INT_VIEW element,
                    ARRAY_REAL_VIEW weights)
      : m_rhs_acoustic(term, element, weights)
  {
  }

  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getTerm(int /*i*/) const override
  {
    return m_rhs_acoustic.getTerm(0);
  }

  PROXY_HOST_DEVICE
  VECTOR_INT_VIEW getElement() const { return m_rhs_acoustic.getElement(); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getWeights() const { return m_rhs_acoustic.getWeights(); }

  void print() const override { m_rhs_acoustic.print(); }

  RhsAcoustic m_rhs_acoustic;  ///< Acoustic source (applied to fluid domain)
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_RHS_ACOUSTOELASTIC_H_

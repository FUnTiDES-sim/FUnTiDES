#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_RHS_ACOUSTOELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_RHS_ACOUSTOELASTIC_H_

#include <data_type.h>

#include "rhs.h"
#include "rhs_acoustic.h"
#include "rhs_elastic.h"

namespace solver
{
namespace fe
{

/**
 * @brief RHS data structure for the acousto-elastic coupled solver.
 *
 * Holds one acoustic source (pressure, fluid domain) and one elastic source
 * (force, solid domain). Either source may be zero-initialised when the
 * corresponding domain is inactive. The sub-solvers receive their respective
 * @ref RhsAcoustic / @ref RhsElastic member directly, so @p getTerm() is
 * provided for interface compliance only.
 */
struct RhsAcoustoElastic : public Rhs
{
  /// Number of RHS components: 1 acoustic (p) + 3 elastic (fx, fy, fz).
  static constexpr int kNumRhsComponents = 4;

  /**
   * @param acoustic_term  2D array of acoustic source signals (n_src x n_t).
   * @param element        Indices of elements containing source points.
   * @param weights        Per-node weights for source distribution.
   * @param elastic_termx  X-component elastic source signals (n_src x n_t).
   * @param elastic_termy  Y-component elastic source signals.
   * @param elastic_termz  Z-component elastic source signals.
   */
  RhsAcoustoElastic(ARRAY_REAL_VIEW acoustic_term, VECTOR_INT_VIEW element,
                    ARRAY_REAL_VIEW weights, ARRAY_REAL_VIEW elastic_termx,
                    ARRAY_REAL_VIEW elastic_termy,
                    ARRAY_REAL_VIEW elastic_termz)
      : m_rhs_acoustic(acoustic_term, element, weights),
        m_rhs_elastic(elastic_termx, elastic_termy, elastic_termz, element,
                      weights)
  {
  }

  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  /// @brief Returns term i: 0 = acoustic, 1/2/3 = elastic x/y/z.
  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getTerm(int i) const override
  {
    if (i == 0) return m_rhs_acoustic.getTerm(0);
    return m_rhs_elastic.getTerm(i - 1);
  }

  PROXY_HOST_DEVICE
  VECTOR_INT_VIEW getElement() const { return m_rhs_acoustic.getElement(); }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getWeights() const { return m_rhs_acoustic.getWeights(); }

  void print() const override
  {
    m_rhs_acoustic.print();
    m_rhs_elastic.print();
  }

  RhsAcoustic m_rhs_acoustic;  ///< Acoustic (fluid) source
  RhsElastic m_rhs_elastic;    ///< Elastic (solid) source
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTOELASTIC_INCLUDE_RHS_ACOUSTOELASTIC_H_

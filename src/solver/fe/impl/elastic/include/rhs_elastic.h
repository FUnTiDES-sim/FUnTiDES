#ifndef FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_RHS_ELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_RHS_ELASTIC_H_
#include <data_type.h>

#include "rhs.h"

namespace solver
{
namespace fe
{
/**
 * @brief Elastic RHS data structure.
 */
struct RhsElastic : public Rhs
{
  /// Number of RHS (source) components
  static constexpr int kNumRhsComponents = 3;

  /// Default constructor: produces an empty (no-forcing) RHS.
  /// All views are null/zero-extent; safe as long as RHS is not accessed.
  RhsElastic() = default;

  RhsElastic(ARRAY_REAL_VIEW termx, ARRAY_REAL_VIEW termy,
             ARRAY_REAL_VIEW termz, VECTOR_INT_VIEW element,
             ARRAY_REAL_VIEW weights)
      : m_termx(termx),
        m_termy(termy),
        m_termz(termz),
        m_element(element),
        m_weights(weights)
  {
  }

  int getNumRhsComponents() const override final { return kNumRhsComponents; }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getTerm(int i) const override
  {
    switch (i)
    {
      case 0:
        return m_termx;
      case 1:
        return m_termy;
      case 2:
        return m_termz;
      default:
        return m_termx;  // make it cuda happy
    }
  }

  PROXY_HOST_DEVICE
  VECTOR_INT_VIEW getElement() const { return m_element; }

  PROXY_HOST_DEVICE
  ARRAY_REAL_VIEW getWeights() const { return m_weights; }

  void print() const override
  {
    std::cout << "RHSx Term size:   " << m_termx.extent(0) << std::endl;
    std::cout << "RHSy Term size:   " << m_termy.extent(0) << std::endl;
    std::cout << "RHSz Term size:   " << m_termz.extent(0) << std::endl;
    std::cout << "RHS Element size: " << m_element.extent(0) << std::endl;
    std::cout << "RHS Weights size: " << m_weights.extent(0) << std::endl;
  }

  ARRAY_REAL_VIEW m_termx;    ///< X-component forcing term
  ARRAY_REAL_VIEW m_termy;    ///< Y-component forcing term
  ARRAY_REAL_VIEW m_termz;    ///< Z-component forcing term
  VECTOR_INT_VIEW m_element;  ///< Source element indices
  ARRAY_REAL_VIEW m_weights;  ///< Forcing weights per node
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_RHS_ELASTIC_H_
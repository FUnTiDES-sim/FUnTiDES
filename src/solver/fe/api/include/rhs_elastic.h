#ifndef RHS_ELASTIC_H_
#define RHS_ELASTIC_H_

#include <data_type.h>
#include "rhs.h"

/**
 * @brief Elastic RHS data structure.
 */
struct RhsElastic : public Rhs
{
  RhsElastic(ARRAY_REAL_VIEW termx,
             ARRAY_REAL_VIEW termy,
             ARRAY_REAL_VIEW termz,
             VECTOR_INT_VIEW element,
             ARRAY_REAL_VIEW weights)
      : m_termx(termx),
        m_termy(termy),
        m_termz(termz),
        m_element(element),
        m_weights(weights)
  {
  }

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

#endif  // RHS_ELASTIC_H_
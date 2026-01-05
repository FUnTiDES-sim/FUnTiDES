#ifndef WAVEFIELD_ELASTIC_H_
#define WAVEFIELD_ELASTIC_H_

#include <data_type.h>
#include "wavefield.h"

/**
 * @brief Elastic wavefield data structure.
 */
struct WavefieldElastic : public Wavefield
{
  WavefieldElastic(VECTOR_REAL_VIEW uxnGlobalPrev,
                   VECTOR_REAL_VIEW uxnGlobalCurr,
                   VECTOR_REAL_VIEW uynGlobalPrev,
                   VECTOR_REAL_VIEW uynGlobalCurr,
                   VECTOR_REAL_VIEW uznGlobalPrev,
                   VECTOR_REAL_VIEW uznGlobalCurr)
      : m_uxnGlobalPrev(uxnGlobalPrev),
        m_uxnGlobalCurr(uxnGlobalCurr),
        m_uynGlobalPrev(uynGlobalPrev),
        m_uynGlobalCurr(uynGlobalCurr),
        m_uznGlobalPrev(uznGlobalPrev),
        m_uznGlobalCurr(uznGlobalCurr)
  {
  }

  void advance() override
  {
    std::swap(m_uxnGlobalPrev, m_uxnGlobalCurr);
    std::swap(m_uynGlobalPrev, m_uynGlobalCurr);
    std::swap(m_uznGlobalPrev, m_uznGlobalCurr);
  }

  void print() const override
  {
    std::cout << "Ux Global Prev size: " << m_uxnGlobalPrev.extent(0) << std::endl;
    std::cout << "Ux Global Curr size: " << m_uxnGlobalCurr.extent(0) << std::endl;
    std::cout << "Uy Global Prev size: " << m_uynGlobalPrev.extent(0) << std::endl;
    std::cout << "Uy Global Curr size: " << m_uynGlobalCurr.extent(0) << std::endl;
    std::cout << "Uz Global Prev size: " << m_uznGlobalPrev.extent(0) << std::endl;
    std::cout << "Uz Global Curr size: " << m_uznGlobalCurr.extent(0) << std::endl;
  }

  VECTOR_REAL_VIEW m_uxnGlobalPrev;  ///< Displacement field in x at previous time step
  VECTOR_REAL_VIEW m_uxnGlobalCurr;  ///< Displacement field in x at current time step
  VECTOR_REAL_VIEW m_uynGlobalPrev;  ///< Displacement field in y at previous time step
  VECTOR_REAL_VIEW m_uynGlobalCurr;  ///< Displacement field in y at current time step
  VECTOR_REAL_VIEW m_uznGlobalPrev;  ///< Displacement field in z at previous time step
  VECTOR_REAL_VIEW m_uznGlobalCurr;  ///< Displacement field in z at current time step
};

#endif  // WAVEFIELD_ELASTIC_H_
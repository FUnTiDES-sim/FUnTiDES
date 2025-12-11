#ifndef WAVEFIELD_ACOUSTIC_H_
#define WAVEFIELD_ACOUSTIC_H_

#include <data_type.h>
#include "wavefield.h"

/**
 * @brief Acoustic wavefield data structure.
 */
struct WavefieldAcoustic : public Wavefield
{
  Wavefield(VECTOR_REAL_VIEW pnGlobalPrev,
            VECTOR_REAL_VIEW pnGlobalCurr)
      : m_pnGlobalPrev(pnGlobalPrev),
        m_pnGlobalCurr(pnGlobalCurr)
  {
  }

  void advance() override
  {
    std::swap(m_pnGlobalPrev, m_pnGlobalCurr);
  }

  void print() const override
  {
    std::cout << "Pn Global Prev size: " << m_pnGlobalPrev.extent(0) << std::endl;
    std::cout << "Pn Global Curr size: " << m_pnGlobalCurr.extent(0) << std::endl;
  }

  ARRAY_REAL_VIEW m_pnGlobalPrev;  ///< Pressure field at previous time step
  ARRAY_REAL_VIEW m_pnGlobalCurr;  ///< Pressure field at current time step
};

#endif  // WAVEFIELD_ACOUSTIC_H_
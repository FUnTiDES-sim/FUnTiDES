#ifndef WAVEFIELD_ACOUSTIC_H_
#define WAVEFIELD_ACOUSTIC_H_

#include <data_type.h>

#include "wavefield.h"

namespace solver
{
namespace fe
{
/**
 * @brief Acoustic wavefield data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct WavefieldAcoustic : public Wavefield
{
  WavefieldAcoustic(VECTOR_REAL_VIEW pnGlobalPrev,
                    VECTOR_REAL_VIEW pnGlobalCurr)
      : m_pnGlobalPrev(pnGlobalPrev), m_pnGlobalCurr(pnGlobalCurr)
  {
  }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  VECTOR_REAL_VIEW getCurrentField(int i) const override
  {
    return m_pnGlobalCurr;
  }

#ifdef USE_KOKKOS
  KOKKOS_FORCEINLINE_FUNCTION
#endif
  VECTOR_REAL_VIEW getPreviousField(int i) const override
  {
    return m_pnGlobalPrev;
  }

  void swap() override { std::swap(m_pnGlobalPrev, m_pnGlobalCurr); }

  void print() const override
  {
    std::cout << "Pn Global Prev size: " << m_pnGlobalPrev.extent(0)
              << std::endl;
    std::cout << "Pn Global Curr size: " << m_pnGlobalCurr.extent(0)
              << std::endl;
  }

  VECTOR_REAL_VIEW m_pnGlobalPrev;  ///< Pressure field at previous time step
  VECTOR_REAL_VIEW m_pnGlobalCurr;  ///< Pressure field at current time step
};
}  // namespace fe
}  // namespace solver

#endif  // WAVEFIELD_ACOUSTIC_H_
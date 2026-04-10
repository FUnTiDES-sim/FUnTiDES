#ifndef FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
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
  /// Number of solution fields (1 for pressure)
  static constexpr int kNumFields = 1;

  /// Primary field name
  static constexpr const char* kFieldNames[1] = {"pressure"};

  WavefieldAcoustic(VECTOR_REAL_VIEW pnGlobalPrev,
                    VECTOR_REAL_VIEW pnGlobalCurr)
      : m_pnGlobalField0(pnGlobalPrev),
        m_pnGlobalField1(pnGlobalCurr),
        m_pnGlobalPrev(&m_pnGlobalField0),
        m_pnGlobalCurr(&m_pnGlobalField1)
  {
  }

  int getNumFields() const override final { return kNumFields; }

  const char* const* getFieldNames() const override final
  {
    return kFieldNames;
  }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getCurrentField(int i) const override
  {
    return *m_pnGlobalCurr;
  }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getPreviousField(int i) const override
  {
    return *m_pnGlobalPrev;
  }

  void swap() override { std::swap(*m_pnGlobalPrev, *m_pnGlobalCurr); }

  void swapWithRotation(VECTOR_REAL_VIEW& prevPrevBuffer, int i) override
  {
    VECTOR_REAL_VIEW tmp = prevPrevBuffer;
    prevPrevBuffer       = *m_pnGlobalPrev;
    *m_pnGlobalPrev      = *m_pnGlobalCurr;
    *m_pnGlobalCurr      = tmp;
  }

  void print() const override
  {
    std::cout << "Pn Global Prev size: " << m_pnGlobalPrev->extent(0)
              << std::endl;
    std::cout << "Pn Global Curr size: " << m_pnGlobalCurr->extent(0)
              << std::endl;
  }

  VECTOR_REAL_VIEW m_pnGlobalField0;        ///< Storage for pressure field 0
  VECTOR_REAL_VIEW m_pnGlobalField1;        ///< Storage for pressure field 1
  VECTOR_REAL_VIEW* m_pnGlobalPrev;         ///< Pointer to previous field
  VECTOR_REAL_VIEW* m_pnGlobalCurr;         ///< Pointer to current field
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ACOUSTIC_INCLUDE_WAVEFIELD_ACOUSTIC_H_
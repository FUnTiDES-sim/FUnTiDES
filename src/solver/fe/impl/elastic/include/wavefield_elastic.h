#ifndef FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_WAVEFIELD_ELASTIC_H_
#define FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_WAVEFIELD_ELASTIC_H_
#include <data_type.h>

#include "wavefield.h"

namespace solver
{
namespace fe
{
/**
 * @brief Elastic wavefield data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct WavefieldElastic : public Wavefield
{
  /// Field names for each component
  static constexpr const char* kFieldNames[3] = {"ux", "uy", "uz"};
  /// Number of solution fields (3 for displacement vector)
  static constexpr int kNumFields = 3;

  WavefieldElastic(VECTOR_REAL_VIEW uxnGlobalPrev,
                   VECTOR_REAL_VIEW uxnGlobalCurr,
                   VECTOR_REAL_VIEW uynGlobalPrev,
                   VECTOR_REAL_VIEW uynGlobalCurr,
                   VECTOR_REAL_VIEW uznGlobalPrev,
                   VECTOR_REAL_VIEW uznGlobalCurr)
      : m_uxnGlobalField0(uxnGlobalPrev),
        m_uxnGlobalField1(uxnGlobalCurr),
        m_uynGlobalField0(uynGlobalPrev),
        m_uynGlobalField1(uynGlobalCurr),
        m_uznGlobalField0(uznGlobalPrev),
        m_uznGlobalField1(uznGlobalCurr),
        m_uxnGlobalPrev(&m_uxnGlobalField0),
        m_uxnGlobalCurr(&m_uxnGlobalField1),
        m_uynGlobalPrev(&m_uynGlobalField0),
        m_uynGlobalCurr(&m_uynGlobalField1),
        m_uznGlobalPrev(&m_uznGlobalField0),
        m_uznGlobalCurr(&m_uznGlobalField1)
  {
  }

  int getNumFields() const override final { return kNumFields; }

  const char* const* getFieldNames() const override final
  {
    return kFieldNames;
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getCurrentField(int i) const override
  {
    switch (i)
    {
      case 0:
        return *m_uxnGlobalCurr;
      case 1:
        return *m_uynGlobalCurr;
      case 2:
        return *m_uznGlobalCurr;
      default:
        return *m_uxnGlobalCurr;  // make it cuda happy
    }
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getPreviousField(int i) const override
  {
    switch (i)
    {
      case 0:
        return *m_uxnGlobalPrev;
      case 1:
        return *m_uynGlobalPrev;
      case 2:
        return *m_uznGlobalPrev;
      default:
        return *m_uxnGlobalCurr;  // make it cuda happy
    }
  }

  void swap() override
  {
    std::swap(*m_uxnGlobalPrev, *m_uxnGlobalCurr);
    std::swap(*m_uynGlobalPrev, *m_uynGlobalCurr);
    std::swap(*m_uznGlobalPrev, *m_uznGlobalCurr);
  }

  // NOTE: elastic has 3 components — the caller must manage one extra buffer per
  // component and call swapWithRotation once per component with the appropriate
  // field index (0=ux, 1=uy, 2=uz).
  void swapWithRotation(VECTOR_REAL_VIEW& prevPrevBuffer, int i) override
  {
    VECTOR_REAL_VIEW tmp = prevPrevBuffer;
    switch (i)
    {
      case 0:  // ux component
        prevPrevBuffer  = *m_uxnGlobalPrev;
        *m_uxnGlobalPrev = *m_uxnGlobalCurr;
        *m_uxnGlobalCurr = tmp;
        break;
      case 1:  // uy component
        prevPrevBuffer  = *m_uynGlobalPrev;
        *m_uynGlobalPrev = *m_uynGlobalCurr;
        *m_uynGlobalCurr = tmp;
        break;
      case 2:  // uz component
        prevPrevBuffer  = *m_uznGlobalPrev;
        *m_uznGlobalPrev = *m_uznGlobalCurr;
        *m_uznGlobalCurr = tmp;
        break;
      default:
        // Invalid field index - no-op to make CUDA happy
        break;
    }
  }

  void print() const override
  {
    std::cout << "Ux Global Prev size: " << m_uxnGlobalPrev->extent(0)
              << std::endl;
    std::cout << "Ux Global Curr size: " << m_uxnGlobalCurr->extent(0)
              << std::endl;
    std::cout << "Uy Global Prev size: " << m_uynGlobalPrev->extent(0)
              << std::endl;
    std::cout << "Uy Global Curr size: " << m_uynGlobalCurr->extent(0)
              << std::endl;
    std::cout << "Uz Global Prev size: " << m_uznGlobalPrev->extent(0)
              << std::endl;
    std::cout << "Uz Global Curr size: " << m_uznGlobalCurr->extent(0)
              << std::endl;
  }

  // Storage for each component
  VECTOR_REAL_VIEW m_uxnGlobalField0;  ///< Storage for ux field 0
  VECTOR_REAL_VIEW m_uxnGlobalField1;  ///< Storage for ux field 1
  VECTOR_REAL_VIEW m_uynGlobalField0;  ///< Storage for uy field 0
  VECTOR_REAL_VIEW m_uynGlobalField1;  ///< Storage for uy field 1
  VECTOR_REAL_VIEW m_uznGlobalField0;  ///< Storage for uz field 0
  VECTOR_REAL_VIEW m_uznGlobalField1;  ///< Storage for uz field 1
  
  // Pointers to current/previous fields for each component
  VECTOR_REAL_VIEW* m_uxnGlobalPrev;  ///< Pointer to ux at previous time step
  VECTOR_REAL_VIEW* m_uxnGlobalCurr;  ///< Pointer to ux at current time step
  VECTOR_REAL_VIEW* m_uynGlobalPrev;  ///< Pointer to uy at previous time step
  VECTOR_REAL_VIEW* m_uynGlobalCurr;  ///< Pointer to uy at current time step
  VECTOR_REAL_VIEW* m_uznGlobalPrev;  ///< Pointer to uz at previous time step
  VECTOR_REAL_VIEW* m_uznGlobalCurr;  ///< Pointer to uz at current time step
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_ELASTIC_INCLUDE_WAVEFIELD_ELASTIC_H_
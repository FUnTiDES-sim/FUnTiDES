#ifndef SOLVER_FE_IMPL_INCLUDE_WAVEFIELD_ELASTIC_H_
#define SOLVER_FE_IMPL_INCLUDE_WAVEFIELD_ELASTIC_H_
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
      : m_uxnGlobalPrev(uxnGlobalPrev),
        m_uxnGlobalCurr(uxnGlobalCurr),
        m_uynGlobalPrev(uynGlobalPrev),
        m_uynGlobalCurr(uynGlobalCurr),
        m_uznGlobalPrev(uznGlobalPrev),
        m_uznGlobalCurr(uznGlobalCurr)
  {
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getCurrentField(int i) const override
  {
    switch (i)
    {
      case 0:
        return m_uxnGlobalCurr;
      case 1:
        return m_uynGlobalCurr;
      case 2:
        return m_uznGlobalCurr;
      default:
        return m_uxnGlobalCurr;  // make it cuda happy
    }
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getPreviousField(int i) const override
  {
    switch (i)
    {
      case 0:
        return m_uxnGlobalPrev;
      case 1:
        return m_uynGlobalPrev;
      case 2:
        return m_uznGlobalPrev;
      default:
        return m_uxnGlobalCurr;  // make it cuda happy
    }
  }

  void swap() override
  {
    std::swap(m_uxnGlobalPrev, m_uxnGlobalCurr);
    std::swap(m_uynGlobalPrev, m_uynGlobalCurr);
    std::swap(m_uznGlobalPrev, m_uznGlobalCurr);
  }

  void print() const override
  {
    std::cout << "Ux Global Prev size: " << m_uxnGlobalPrev.extent(0)
              << std::endl;
    std::cout << "Ux Global Curr size: " << m_uxnGlobalCurr.extent(0)
              << std::endl;
    std::cout << "Uy Global Prev size: " << m_uynGlobalPrev.extent(0)
              << std::endl;
    std::cout << "Uy Global Curr size: " << m_uynGlobalCurr.extent(0)
              << std::endl;
    std::cout << "Uz Global Prev size: " << m_uznGlobalPrev.extent(0)
              << std::endl;
    std::cout << "Uz Global Curr size: " << m_uznGlobalCurr.extent(0)
              << std::endl;
  }

  VECTOR_REAL_VIEW
  m_uxnGlobalPrev;  ///< Displacement field in x at previous time step
  VECTOR_REAL_VIEW
  m_uxnGlobalCurr;  ///< Displacement field in x at current time step
  VECTOR_REAL_VIEW
  m_uynGlobalPrev;  ///< Displacement field in y at previous time step
  VECTOR_REAL_VIEW
  m_uynGlobalCurr;  ///< Displacement field in y at current time step
  VECTOR_REAL_VIEW
  m_uznGlobalPrev;  ///< Displacement field in z at previous time step
  VECTOR_REAL_VIEW
  m_uznGlobalCurr;  ///< Displacement field in z at current time step
};
}  // namespace fe
}  // namespace solver
#endif  // SOLVER_FE_IMPL_INCLUDE_WAVEFIELD_ELASTIC_H_
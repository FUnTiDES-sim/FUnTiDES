#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

#include <iostream>

#include "wavefield_view.h"

namespace gradient
{

/**
 * @brief Read-only view of an acoustic adjoint wavefield for gradient computation.
 *
 * Exposes the current adjoint pressure qn and its pre-computed second-order
 * time derivative qdt2 = (q_prev - 2*q_curr + qnext) / dt², both required
 * by the acoustic gradient kernel.
 *
 * Constructed by the caller from solver adjoint wavefield data; no solver
 * dependency here.
 *
 * Fields:
 *   getField(0) = qn    (current adjoint pressure)
 *   getField(1) = qdt2  (second-order time derivative of adjoint pressure)
 */
struct WavefieldViewBackwardAcoustic : public WavefieldView
{
  static constexpr int kNumFields = 2;
  static constexpr const char* kFieldNames[2] = {"qn", "qdt2"};

  WavefieldViewBackwardAcoustic(VECTOR_REAL_VIEW qn, VECTOR_REAL_VIEW qdt2)
      : m_qn(qn), m_qdt2(qdt2)
  {
  }

  int getNumFields() const override { return kNumFields; }

  const char* const* getFieldNames() const override { return kFieldNames; }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getField(int i) const override
  {
    switch (i)
    {
      case 0:
        return m_qn;
      case 1:
        return m_qdt2;
      default:
        return m_qn;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "WavefieldViewBackwardAcoustic:"
              << " qn size="   << m_qn.extent(0)
              << " qdt2 size=" << m_qdt2.extent(0) << "\n";
  }

  VECTOR_REAL_VIEW m_qn;    ///< Current adjoint pressure snapshot
  VECTOR_REAL_VIEW m_qdt2;  ///< Second-order time derivative of adjoint pressure
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

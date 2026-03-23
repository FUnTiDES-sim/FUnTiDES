#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

#include <iostream>

#include "common_macros.h"
#include "data_type.h"
#include "wavefield_view.h"

namespace gradient
{

/**
 * @brief Read-only view of an acoustic adjoint wavefield for gradient computation.
 *
 * Exposes the current adjoint pressure q_n and its pre-computed second-order
 * time derivative q_dt2 = (q_prev - 2*q_curr + q_next) / dt², both required
 * by the acoustic gradient kernel.
 *
 * Constructed by the caller from solver adjoint wavefield data; no solver
 * dependency here.
 *
 * Fields:
 *   getField(0) = q_n    (current adjoint pressure)
 *   getField(1) = q_dt2  (second-order time derivative of adjoint pressure)
 */
struct WavefieldViewBackwardAcoustic : public WavefieldView
{
  static constexpr int kNumFields = 2;
  static constexpr const char* kFieldNames[2] = {"q_n", "q_dt2"};

  WavefieldViewBackwardAcoustic(VECTOR_REAL_VIEW q_n, VECTOR_REAL_VIEW q_dt2)
      : m_q_n(q_n), m_q_dt2(q_dt2)
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
        return m_q_n;
      case 1:
        return m_q_dt2;
      default:
        return m_q_n;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "WavefieldViewBackwardAcoustic:"
              << " q_n size="   << m_q_n.extent(0)
              << " q_dt2 size=" << m_q_dt2.extent(0) << "\n";
  }

  VECTOR_REAL_VIEW m_q_n;    ///< Current adjoint pressure snapshot
  VECTOR_REAL_VIEW m_q_dt2;  ///< Second-order time derivative of adjoint pressure
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

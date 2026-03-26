#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

#include <iostream>
#include <string>

#include "wavefield_view.h"

namespace gradient
{

/**
 * @brief Read-only view of an acoustic adjoint wavefield for gradient
 * computation.
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
class WavefieldViewBackwardAcoustic : public WavefieldView
{
public:

  WavefieldViewBackwardAcoustic(VECTOR_REAL_VIEW qn, VECTOR_REAL_VIEW qdt2)
      : qn_(qn), qdt2_(qdt2)
  {
  }

  int getNumFields() const override { return kNumFields; }

  std::string getFieldName(int i) const override
  {
    switch (i)
    {
      case 0:
        return "qn";
      case 1:
        return "qdt2";
      default:
        return "qn";
    }
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getField(int i) const override
  {
    switch (i)
    {
      case 0:
        return qn_;
      case 1:
        return qdt2_;
      default:
        return qn_;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "WavefieldViewBackwardAcoustic:"
              << " qn size=" << qn_.extent(0)
              << " qdt2 size=" << qdt2_.extent(0) << "\n";
  }

  private:

    static constexpr int kNumFields = 2;

    VECTOR_REAL_VIEW qn_;  ///< Current adjoint pressure snapshot
    VECTOR_REAL_VIEW qdt2_; ///< Second-order time derivative of adjoint pressure
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

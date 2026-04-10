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
 * Exposes three consecutive adjoint pressure snapshots required by the
 * acoustic gradient kernel to compute the second-order time derivative on
 * the fly: qdt2 = (qnPrevPrev - 2*qnPrev + qn) / dt².
 *
 * Constructed by the caller from solver adjoint wavefield data plus an
 * externally managed third buffer; no solver dependency here.
 *
 * In a time loop, the solver's WavefieldAcoustic provides qn (current) and
 * qnPrev (previous). The caller must allocate one extra buffer qnPrevPrev
 * and rotate the three handles each step without any data copy — see the
 * three-buffer ring pattern described in the solver time-loop documentation.
 *
 * Fields:
 *   getField(0) = qn          (adjoint pressure at time n)
 *   getField(1) = qnPrev      (adjoint pressure at time n-1)
 *   getField(2) = qnPrevPrev  (adjoint pressure at time n-2)
 */
class WavefieldViewBackwardAcoustic : public WavefieldView
{
 public:
  static constexpr int kNumFields = 3;

  WavefieldViewBackwardAcoustic(VECTOR_REAL_VIEW qn, VECTOR_REAL_VIEW qnPrev,
                                VECTOR_REAL_VIEW qnPrevPrev)
      : qn_(qn), qnPrev_(qnPrev), qnPrevPrev_(qnPrevPrev)
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
        return "qnPrev";
      case 2:
        return "qnPrevPrev";
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
        return qnPrev_;
      case 2:
        return qnPrevPrev_;
      default:
        return qn_;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "WavefieldViewBackwardAcoustic:" << " qn size="
              << qn_.extent(0) << " qnPrev size=" << qnPrev_.extent(0)
              << " qnPrevPrev size=" << qnPrevPrev_.extent(0) << "\n";
  }

 private:
  VECTOR_REAL_VIEW qn_;          ///< Adjoint pressure at time n
  VECTOR_REAL_VIEW qnPrev_;      ///< Adjoint pressure at time n-1
  VECTOR_REAL_VIEW qnPrevPrev_;  ///< Adjoint pressure at time n-2
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

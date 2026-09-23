#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

#include <iostream>
#include <string>

#include "wavefield_view.h"

namespace gradient {

/**
 * @brief Read-only view of three consecutive acoustic adjoint pressure snapshots.
 *
 * The acoustic gradient kernel uses the three snapshots to build the adjoint
 * second time derivative on the fly:
 * qdt2 = (qnPrevPrev - 2*qnPrev + qn) / dt^2.
 *
 * The view only stores handles to buffers owned by the caller; it never copies
 * data. In a time loop the caller rotates the three handles at each step and
 * rebuilds the view.
 *
 * Field indices:
 * - getField(0) = qn         (adjoint pressure at time n)
 * - getField(1) = qnPrev     (adjoint pressure at time n-1)
 * - getField(2) = qnPrevPrev (adjoint pressure at time n-2)
 */
class WavefieldViewBackwardAcoustic : public WavefieldView {
 public:
  static constexpr int kNumFields = 3;  ///< Number of exposed fields (qn, qnPrev, qnPrevPrev).

  /**
   * @brief Builds a view from three externally owned snapshots.
   * @param[in] qn Adjoint pressure at time n.
   * @param[in] qnPrev Adjoint pressure at time n-1.
   * @param[in] qnPrevPrev Adjoint pressure at time n-2.
   */
  WavefieldViewBackwardAcoustic(vectorReal qn, vectorReal qnPrev, vectorReal qnPrevPrev)
      : qn_(qn), qnPrev_(qnPrev), qnPrevPrev_(qnPrevPrev) {}

  /** @brief Returns the number of exposed fields (kNumFields). */
  int getNumFields() const override { return kNumFields; }

  /**
   * @brief Returns the name of field i ("qn", "qnPrev" or "qnPrevPrev").
   * @param[in] i Field index in [0, kNumFields).
   */
  std::string getFieldName(int i) const override {
    switch (i) {
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

  /**
   * @brief Returns the handle of field i.
   * @param[in] i Field index in [0, kNumFields).
   * @return The snapshot handle; no data is copied.
   */
  // TODO: replace the switch by a template and constexpr if once C++20 is available.
  PROXY_HOST_DEVICE
  vectorReal getField(int i) const override {
    switch (i) {
      case 0:
        return qn_;
      case 1:
        return qnPrev_;
      case 2:
        return qnPrevPrev_;
      default:
        return qn_;  // the CUDA compiler requires a return on every path
    }
  }

  /** @brief Prints the size of each snapshot to stdout. */
  void print() const override {
    std::cout << "WavefieldViewBackwardAcoustic:" << " qn size=" << qn_.extent(0)
              << " qnPrev size=" << qnPrev_.extent(0) << " qnPrevPrev size=" << qnPrevPrev_.extent(0) << "\n";
  }

 private:
  vectorReal qn_;          ///< Adjoint pressure at time n
  vectorReal qnPrev_;      ///< Adjoint pressure at time n-1
  vectorReal qnPrevPrev_;  ///< Adjoint pressure at time n-2
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ACOUSTIC_H_

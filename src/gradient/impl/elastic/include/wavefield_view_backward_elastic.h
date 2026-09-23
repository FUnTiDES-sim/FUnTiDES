#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_

#include <iostream>
#include <string>
#include <vector>

#include "wavefield_view.h"

namespace gradient {

/**
 * @brief Read-only view of an elastic adjoint wavefield, holding the adjoint
 * displacement and its precomputed second time derivative.
 *
 * The view stores shallow copies of the vectors given at construction and
 * does not own their data. Field indices:
 *  - 0, 1, 2: adjoint displacement ux_n, uy_n, uz_n
 *  - 3, 4, 5: second time derivative ux_dt2, uy_dt2, uz_dt2
 */
class WavefieldViewBackwardElastic : public WavefieldView {
 public:
  static constexpr int kNumFields = 6;  ///< Number of fields exposed by the view.

  /**
   * @brief Builds the view from the six adjoint vectors.
   * @param[in] ux_n Adjoint displacement, x component.
   * @param[in] uy_n Adjoint displacement, y component.
   * @param[in] uz_n Adjoint displacement, z component.
   * @param[in] ux_dt2 Second time derivative of the adjoint displacement, x component.
   * @param[in] uy_dt2 Second time derivative of the adjoint displacement, y component.
   * @param[in] uz_dt2 Second time derivative of the adjoint displacement, z component.
   */
  WavefieldViewBackwardElastic(vectorReal ux_n, vectorReal uy_n, vectorReal uz_n, vectorReal ux_dt2, vectorReal uy_dt2,
                               vectorReal uz_dt2)
      : ux_n_(ux_n), uy_n_(uy_n), uz_n_(uz_n), ux_dt2_(ux_dt2), uy_dt2_(uy_dt2), uz_dt2_(uz_dt2) {}

  /** @brief Returns the number of fields, kNumFields. */
  int getNumFields() const override { return kNumFields; }

  /**
   * @brief Returns the name of field i.
   * @param[in] i Field index in [0, kNumFields).
   */
  std::string getFieldName(int i) const override {
    switch (i) {
      case 0:
        return "ux_n";
      case 1:
        return "uy_n";
      case 2:
        return "uz_n";
      case 3:
        return "ux_dt2";
      case 4:
        return "uy_dt2";
      case 5:
        return "uz_dt2";
      default:
        return "ux_n";
    }
  }

  // TODO use template + constexpr if when C++20 is available
  /**
   * @brief Returns field i.
   * @param[in] i Field index in [0, kNumFields).
   * @return The vector of field i.
   */
  PROXY_HOST_DEVICE
  vectorReal getField(int i) const override {
    switch (i) {
      case 0:
        return ux_n_;
      case 1:
        return uy_n_;
      case 2:
        return uz_n_;
      case 3:
        return ux_dt2_;
      case 4:
        return uy_dt2_;
      case 5:
        return uz_dt2_;
      default:
        return ux_n_;  // make it cuda happy
    }
  }

  /** @brief Prints the size of each field to stdout. */
  void print() const override {
    std::cout << "WavefieldViewBackwardElastic:" << " ux_n size=" << ux_n_.extent(0) << " uy_n size=" << uy_n_.extent(0)
              << " uz_n size=" << uz_n_.extent(0) << " ux_dt2 size=" << ux_dt2_.extent(0)
              << " uy_dt2 size=" << uy_dt2_.extent(0) << " uz_dt2 size=" << uz_dt2_.extent(0) << "\n";
  }

 private:
  vectorReal ux_n_;    ///< Adjoint displacement, x component.
  vectorReal uy_n_;    ///< Adjoint displacement, y component.
  vectorReal uz_n_;    ///< Adjoint displacement, z component.
  vectorReal ux_dt2_;  ///< Second time derivative, x component.
  vectorReal uy_dt2_;  ///< Second time derivative, y component.
  vectorReal uz_dt2_;  ///< Second time derivative, z component.
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_

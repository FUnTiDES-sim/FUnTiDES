#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_

#include <iostream>
#include <string>
#include <vector>

#include "wavefield_view.h"

namespace gradient {

/**
 * @brief Read-only view of an elastic adjoint wavefield for gradient
 * computation.
 *
 * Exposes the current adjoint displacement [ux_n, uy_n, uz_n] and their
 * pre-computed second-order time derivatives [ux_dt2, uy_dt2, uz_dt2] needed
 * by the elastic gradient kernel.
 *
 * Constructed by the caller from solver adjoint wavefield data; no solver
 * dependency here.
 *
 * Fields:
 *   getField(0) = ux_n    (current x-adjoint displacement)
 *   getField(1) = uy_n    (current y-adjoint displacement)
 *   getField(2) = uz_n    (current z-adjoint displacement)
 *   getField(3) = ux_dt2  (second-order time derivative, x-component)
 *   getField(4) = uy_dt2  (second-order time derivative, y-component)
 *   getField(5) = uz_dt2  (second-order time derivative, z-component)
 */
class WavefieldViewBackwardElastic : public WavefieldView {
 public:
  static constexpr int kNumFields = 6;

  WavefieldViewBackwardElastic(vectorReal ux_n, vectorReal uy_n, vectorReal uz_n, vectorReal ux_dt2, vectorReal uy_dt2,
                               vectorReal uz_dt2)
      : ux_n_(ux_n), uy_n_(uy_n), uz_n_(uz_n), ux_dt2_(ux_dt2), uy_dt2_(uy_dt2), uz_dt2_(uz_dt2) {}

  int getNumFields() const override { return kNumFields; }

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

  void print() const override {
    std::cout << "WavefieldViewBackwardElastic:" << " ux_n size=" << ux_n_.extent(0) << " uy_n size=" << uy_n_.extent(0)
              << " uz_n size=" << uz_n_.extent(0) << " ux_dt2 size=" << ux_dt2_.extent(0)
              << " uy_dt2 size=" << uy_dt2_.extent(0) << " uz_dt2 size=" << uz_dt2_.extent(0) << "\n";
  }

 private:
  vectorReal ux_n_;    ///< Current x-adjoint displacement snapshot
  vectorReal uy_n_;    ///< Current y-adjoint displacement snapshot
  vectorReal uz_n_;    ///< Current z-adjoint displacement snapshot
  vectorReal ux_dt2_;  ///< Second-order time derivative, x-component
  vectorReal uy_dt2_;  ///< Second-order time derivative, y-component
  vectorReal uz_dt2_;  ///< Second-order time derivative, z-component
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_

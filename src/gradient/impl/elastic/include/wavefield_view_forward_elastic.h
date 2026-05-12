#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_

#include <iostream>
#include <string>

#include "wavefield_view.h"

namespace gradient {

/**
 * @brief Read-only view of an elastic forward wavefield for gradient
 * computation.
 *
 * Exposes only the current displacement snapshots [ux_n, uy_n, uz_n] needed
 * by the elastic gradient kernel. Constructed by the caller from a solver
 * wavefield; no solver dependency here.
 *
 * Fields:
 *   getField(0) = ux_n  (current x-displacement)
 *   getField(1) = uy_n  (current y-displacement)
 *   getField(2) = uz_n  (current z-displacement)
 */
class WavefieldViewForwardElastic : public WavefieldView {
 public:
  static constexpr int kNumFields = 3;

  WavefieldViewForwardElastic(vectorReal ux_n, vectorReal uy_n, vectorReal uz_n)
      : ux_n_(ux_n), uy_n_(uy_n), uz_n_(uz_n) {}

  int getNumFields() const override { return kNumFields; }

  std::string getFieldName(int i) const override {
    switch (i) {
      case 0:
        return "ux_n";
      case 1:
        return "uy_n";
      case 2:
        return "uz_n";
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
      default:
        return ux_n_;  // make it cuda happy
    }
  }

  void print() const override {
    std::cout << "WavefieldViewForwardElastic:" << " ux_n size=" << ux_n_.extent(0) << " uy_n size=" << uy_n_.extent(0)
              << " uz_n size=" << uz_n_.extent(0) << "\n";
  }

 private:
  vectorReal ux_n_;  ///< Current x-displacement snapshot
  vectorReal uy_n_;  ///< Current y-displacement snapshot
  vectorReal uz_n_;  ///< Current z-displacement snapshot
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_

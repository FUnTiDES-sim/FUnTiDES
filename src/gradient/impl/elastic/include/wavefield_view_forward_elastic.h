#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_

#include <iostream>
#include <string>

#include "wavefield_view.h"

namespace gradient {

/**
 * @brief Read-only view of the elastic forward wavefield used by the gradient computation.
 *
 * Holds the current displacement snapshots ux_n, uy_n, uz_n. The view stores
 * the given arrays as they are and does not depend on any solver.
 *
 * Field indices: 0 = ux_n, 1 = uy_n, 2 = uz_n.
 */
class WavefieldViewForwardElastic : public WavefieldView {
 public:
  /// Number of fields exposed by the view.
  static constexpr int kNumFields = 3;

  /**
   * @brief Builds a view over the three displacement components.
   * @param[in] ux_n Current x-displacement snapshot.
   * @param[in] uy_n Current y-displacement snapshot.
   * @param[in] uz_n Current z-displacement snapshot.
   */
  WavefieldViewForwardElastic(vectorReal ux_n, vectorReal uy_n, vectorReal uz_n)
      : ux_n_(ux_n), uy_n_(uy_n), uz_n_(uz_n) {}

  /// @return Number of fields (3).
  int getNumFields() const override { return kNumFields; }

  /**
   * @brief Returns the name of a field.
   * @param[in] i Field index in [0, 2].
   * @return "ux_n", "uy_n" or "uz_n".
   */
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

  /**
   * @brief Returns a displacement snapshot.
   * @param[in] i Field index in [0, 2].
   * @return The corresponding array (shallow copy of the stored view).
   */
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

  /// @brief Prints the size of each displacement array to stdout.
  void print() const override {
    std::cout << "WavefieldViewForwardElastic:" << " ux_n size=" << ux_n_.extent(0) << " uy_n size=" << uy_n_.extent(0)
              << " uz_n size=" << uz_n_.extent(0) << "\n";
  }

 private:
  vectorReal ux_n_;  ///< Current x-displacement snapshot.
  vectorReal uy_n_;  ///< Current y-displacement snapshot.
  vectorReal uz_n_;  ///< Current z-displacement snapshot.
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_

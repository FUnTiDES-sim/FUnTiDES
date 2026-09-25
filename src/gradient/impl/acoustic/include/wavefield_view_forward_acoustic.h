#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

#include <iostream>
#include <string>

#include "wavefield_view.h"

namespace gradient {

/**
 * @brief Read-only view of an acoustic forward wavefield for gradient
 * computation.
 *
 * Exposes a single field, the current pressure snapshot pn (index 0), which is
 * the only forward quantity the gradient kernel needs. It holds a copy of the
 * vector handle passed at construction, not the data.
 */
class WavefieldViewForwardAcoustic : public WavefieldView {
 public:
  static constexpr int kNumFields = 1;  ///< Number of exposed fields.

  /**
   * @brief Wraps a pressure snapshot.
   * @param[in] pn Current pressure snapshot.
   */
  WavefieldViewForwardAcoustic(vectorReal pn) : pn_(pn) {}

  /** @brief Returns the number of exposed fields (1). */
  int getNumFields() const override { return kNumFields; }

  /** @brief Returns the name of field i ("pn"); i is ignored. */
  std::string getFieldName(int i) const override { return "pn"; }

  // TODO use template + constexpr if when C++20 is available
  /**
   * @brief Returns the pressure snapshot pn; i is ignored.
   * @param[in] i Field index (only 0 is meaningful).
   */
  PROXY_HOST_DEVICE
  vectorReal getField(int i) const override { return pn_; }

  /** @brief Prints the size of the pressure vector to stdout. */
  void print() const override { std::cout << "WavefieldViewForwardAcoustic: pn size=" << pn_.extent(0) << "\n"; }

 private:
  vectorReal pn_;  ///< Current pressure snapshot.
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

#include <iostream>
#include <string>

#include "wavefield_view.h"

namespace gradient
{

/**
 * @brief Read-only view of an acoustic forward wavefield for gradient
 * computation.
 *
 * Exposes only the current pressure snapshot pn needed by the gradient kernel.
 * Constructed by the caller from a solver wavefield; no solver dependency here.
 *
 * Fields:
 *   getField(0) = pn  (current pressure)
 */
class WavefieldViewForwardAcoustic : public WavefieldView
{
public:

  static constexpr int kNumFields = 1;

  WavefieldViewForwardAcoustic(VECTOR_REAL_VIEW pn) : pn_(pn) {}

  int getNumFields() const override { return kNumFields; }

  std::string getFieldName(int i) const override { return "pn"; }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getField(int i) const override { return pn_; }

  void print() const override
  {
    std::cout << "WavefieldViewForwardAcoustic: pn size=" << pn_.extent(0)
              << "\n";
  }

  private:

    VECTOR_REAL_VIEW pn_;  ///< Current pressure snapshot
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

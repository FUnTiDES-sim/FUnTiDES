#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

#include <iostream>

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
struct WavefieldViewForwardAcoustic : public WavefieldView
{
  static constexpr int kNumFields = 1;
  static constexpr const char* kFieldNames[1] = {"pn"};

  explicit WavefieldViewForwardAcoustic(VECTOR_REAL_VIEW pn) : m_pn(pn) {}

  int getNumFields() const override { return kNumFields; }

  const char* const* getFieldNames() const override { return kFieldNames; }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getField(int i) const override { return m_pn; }

  void print() const override
  {
    std::cout << "WavefieldViewForwardAcoustic: pn size=" << m_pn.extent(0)
              << "\n";
  }

  VECTOR_REAL_VIEW m_pn;  ///< Current pressure snapshot
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

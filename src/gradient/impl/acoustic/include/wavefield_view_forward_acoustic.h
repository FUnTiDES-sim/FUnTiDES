#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

#include <iostream>

#include "common_macros.h"
#include "data_type.h"
#include "wavefield_view.h"

namespace gradient
{

/**
 * @brief Read-only view of an acoustic forward wavefield for gradient computation.
 *
 * Exposes only the current pressure snapshot p_n needed by the gradient kernel.
 * Constructed by the caller from a solver wavefield; no solver dependency here.
 *
 * Fields:
 *   getField(0) = p_n  (current pressure)
 */
struct WavefieldViewForwardAcoustic : public WavefieldView
{
  static constexpr int kNumFields = 1;
  static constexpr const char* kFieldNames[1] = {"p_n"};

  explicit WavefieldViewForwardAcoustic(VECTOR_REAL_VIEW p_n)
      : m_p_n(p_n)
  {
  }

  int getNumFields() const override { return kNumFields; }

  const char* const* getFieldNames() const override { return kFieldNames; }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getField(int /*i*/) const override { return m_p_n; }

  void print() const override
  {
    std::cout << "WavefieldViewForwardAcoustic: p_n size=" << m_p_n.extent(0)
              << "\n";
  }

  VECTOR_REAL_VIEW m_p_n;  ///< Current pressure snapshot
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ACOUSTIC_H_

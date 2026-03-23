#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_

#include <iostream>

#include "common_macros.h"
#include "data_type.h"
#include "wavefield_view.h"

namespace gradient
{

/**
 * @brief Read-only view of an elastic forward wavefield for gradient computation.
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
struct WavefieldViewForwardElastic : public WavefieldView
{
  static constexpr int kNumFields = 3;
  static constexpr const char* kFieldNames[3] = {"ux_n", "uy_n", "uz_n"};

  WavefieldViewForwardElastic(VECTOR_REAL_VIEW ux_n,
                              VECTOR_REAL_VIEW uy_n,
                              VECTOR_REAL_VIEW uz_n)
      : m_ux_n(ux_n), m_uy_n(uy_n), m_uz_n(uz_n)
  {
  }

  int getNumFields() const override { return kNumFields; }

  const char* const* getFieldNames() const override { return kFieldNames; }

  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getField(int i) const override
  {
    switch (i)
    {
      case 0:
        return m_ux_n;
      case 1:
        return m_uy_n;
      case 2:
        return m_uz_n;
      default:
        return m_ux_n;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "WavefieldViewForwardElastic:"
              << " ux_n size=" << m_ux_n.extent(0)
              << " uy_n size=" << m_uy_n.extent(0)
              << " uz_n size=" << m_uz_n.extent(0) << "\n";
  }

  VECTOR_REAL_VIEW m_ux_n;  ///< Current x-displacement snapshot
  VECTOR_REAL_VIEW m_uy_n;  ///< Current y-displacement snapshot
  VECTOR_REAL_VIEW m_uz_n;  ///< Current z-displacement snapshot
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_FORWARD_ELASTIC_H_

#ifndef FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_

#include <iostream>

#include "wavefield_view.h"

namespace gradient
{

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
struct WavefieldViewBackwardElastic : public WavefieldView
{
  static constexpr int kNumFields = 6;
  static constexpr const char* kFieldNames[6] = {"ux_n",   "uy_n",   "uz_n",
                                                 "ux_dt2", "uy_dt2", "uz_dt2"};

  WavefieldViewBackwardElastic(VECTOR_REAL_VIEW ux_n, VECTOR_REAL_VIEW uy_n,
                               VECTOR_REAL_VIEW uz_n, VECTOR_REAL_VIEW ux_dt2,
                               VECTOR_REAL_VIEW uy_dt2, VECTOR_REAL_VIEW uz_dt2)
      : m_ux_n(ux_n),
        m_uy_n(uy_n),
        m_uz_n(uz_n),
        m_ux_dt2(ux_dt2),
        m_uy_dt2(uy_dt2),
        m_uz_dt2(uz_dt2)
  {
  }

  int getNumFields() const override { return kNumFields; }

  const char* const* getFieldNames() const override { return kFieldNames; }

  // TODO use template + constexpr if when C++20 is available
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
      case 3:
        return m_ux_dt2;
      case 4:
        return m_uy_dt2;
      case 5:
        return m_uz_dt2;
      default:
        return m_ux_n;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "WavefieldViewBackwardElastic:"
              << " ux_n size=" << m_ux_n.extent(0)
              << " uy_n size=" << m_uy_n.extent(0)
              << " uz_n size=" << m_uz_n.extent(0)
              << " ux_dt2 size=" << m_ux_dt2.extent(0)
              << " uy_dt2 size=" << m_uy_dt2.extent(0)
              << " uz_dt2 size=" << m_uz_dt2.extent(0) << "\n";
  }

  VECTOR_REAL_VIEW m_ux_n;    ///< Current x-adjoint displacement snapshot
  VECTOR_REAL_VIEW m_uy_n;    ///< Current y-adjoint displacement snapshot
  VECTOR_REAL_VIEW m_uz_n;    ///< Current z-adjoint displacement snapshot
  VECTOR_REAL_VIEW m_ux_dt2;  ///< Second-order time derivative, x-component
  VECTOR_REAL_VIEW m_uy_dt2;  ///< Second-order time derivative, y-component
  VECTOR_REAL_VIEW m_uz_dt2;  ///< Second-order time derivative, z-component
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_WAVEFIELD_VIEW_BACKWARD_ELASTIC_H_

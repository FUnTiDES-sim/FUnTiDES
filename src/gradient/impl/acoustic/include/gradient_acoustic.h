#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_

#include <iostream>

#include "gradient.h"

namespace gradient
{
/**
 * @brief Acoustic gradient data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
struct GradientAcoustic : public Gradient
{
  /// Number of gradients for model inversions (2 for Kappa and Buoyancy)
  static constexpr int kNumGrads = 2;

  /// Primary field name
  static constexpr const char* kGradsNames[2] = {"gradKappa", "gradBuoyancy"};

  GradientAcoustic(VECTOR_REAL_VIEW gradKappa, VECTOR_REAL_VIEW gradBuoyancy)
      : m_gradKappa(gradKappa), m_gradBuoyancy(gradBuoyancy)
  {
  }

  int getNumGradients() const override final { return kNumGrads; }

  const char* const* getGradientNames() const override final
  {
    return kGradsNames;
  }

  // TODO use template + constexpr if when C++20 is available
  PROXY_HOST_DEVICE
  VECTOR_REAL_VIEW getGradient(int i) const override
  {
    switch (i)
    {
      case 0:
        return m_gradKappa;
      case 1:
        return m_gradBuoyancy;
      default:
        return m_gradKappa;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "Grad Kappa size: " << m_gradKappa.extent(0) << std::endl;
    std::cout << "Grad Buoyancy size: " << m_gradBuoyancy.extent(0)
              << std::endl;
  }

  VECTOR_REAL_VIEW m_gradKappa;     ///< Gradient of Kappa
  VECTOR_REAL_VIEW m_gradBuoyancy;  ///< Gradient of Buoyancy
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_

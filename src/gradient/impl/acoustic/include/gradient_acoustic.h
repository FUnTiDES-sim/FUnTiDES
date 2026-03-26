#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_

#include <iostream>
#include <string>

#include "gradient.h"

namespace gradient
{
/**
 * @brief Acoustic gradient data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
class GradientAcoustic : public Gradient
{
public:

  GradientAcoustic(VECTOR_REAL_VIEW gradKappa, VECTOR_REAL_VIEW gradBuoyancy)
      : m_gradKappa(gradKappa), m_gradBuoyancy(gradBuoyancy)
  {
  }

  int getNumGradients() const override final { return kNumGrads; }

  // TODO use template + constexpr if when C++20 is available
  std::string getGradientName(int i) const override final
  {
    switch (i)
    {
      case 0:
        return "gradKappa";
      case 1:
        return "gradBuoyancy";
      default:
        return "gradKappa";
    }
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

  private:
    /// Number of gradients for model inversions (2 for Kappa and Buoyancy)
    static constexpr int kNumGrads = 2;
    static constexpr const std::vector<std::string> kGradientNames = {"gradKappa", "gradBuoyancy"};

    VECTOR_REAL_VIEW m_gradKappa;     ///< Gradient of Kappa
    VECTOR_REAL_VIEW m_gradBuoyancy;  ///< Gradient of Buoyancy
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_

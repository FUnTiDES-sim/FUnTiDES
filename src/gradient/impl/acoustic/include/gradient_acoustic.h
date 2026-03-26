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

  static constexpr int kNumGrads = 2;

  GradientAcoustic(VECTOR_REAL_VIEW gradKappa, VECTOR_REAL_VIEW gradBuoyancy)
      : gradKappa_(gradKappa), gradBuoyancy_(gradBuoyancy)
  {
  }

  int getNumGradients() const override final { return kNumGrads; }

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
        return gradKappa_;
      case 1:
        return gradBuoyancy_;
      default:
        return gradKappa_;  // make it cuda happy
    }
  }

  void print() const override
  {
    std::cout << "Grad Kappa size: " << gradKappa_.extent(0) << std::endl;
    std::cout << "Grad Buoyancy size: " << gradBuoyancy_.extent(0)
              << std::endl;
  }

  private:

    VECTOR_REAL_VIEW gradKappa_;     ///< Gradient of Kappa
    VECTOR_REAL_VIEW gradBuoyancy_;  ///< Gradient of Buoyancy
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_

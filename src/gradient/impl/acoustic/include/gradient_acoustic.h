#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_

#include <iostream>
#include <string>

#include "gradient.h"

namespace gradient {
/**
 * @brief Gradients of the acoustic misfit with respect to kappa and buoyancy.
 *
 * Holds two flat vectors, gradient 0 (kappa) and gradient 1 (buoyancy).
 * The vectors are shared views: copies of this object alias the same storage.
 */
class GradientAcoustic : public Gradient {
 public:
  static constexpr int kNumGrads = 2;  ///< Number of gradient fields.

  /**
   * @brief Wraps the two gradient vectors.
   * @param[in] gradKappa Gradient with respect to kappa.
   * @param[in] gradBuoyancy Gradient with respect to buoyancy.
   */
  GradientAcoustic(vectorReal gradKappa, vectorReal gradBuoyancy)
      : gradKappa_(gradKappa), gradBuoyancy_(gradBuoyancy) {}

  /// @return Number of gradient fields (kNumGrads).
  int getNumGradients() const override final { return kNumGrads; }

  /**
   * @brief Name of a gradient field.
   * @param[in] i Gradient index, 0 for kappa and 1 for buoyancy.
   * @return "gradKappa" or "gradBuoyancy"; "gradKappa" for any other index.
   */
  std::string getGradientName(int i) const override final {
    switch (i) {
      case 0:
        return "gradKappa";
      case 1:
        return "gradBuoyancy";
      default:
        return "gradKappa";
    }
  }

  // TODO use template + constexpr if when C++20 is available
  /**
   * @brief Returns a gradient vector.
   * @param[in] i Gradient index, 0 for kappa and 1 for buoyancy.
   * @return The kappa gradient for any index other than 1.
   */
  PROXY_HOST_DEVICE
  vectorReal getGradient(int i) const override {
    switch (i) {
      case 0:
        return gradKappa_;
      case 1:
        return gradBuoyancy_;
      default:
        return gradKappa_;  // make it cuda happy
    }
  }

  /// @brief Prints the size of each gradient vector to stdout.
  void print() const override {
    std::cout << "Grad Kappa size: " << gradKappa_.extent(0) << std::endl;
    std::cout << "Grad Buoyancy size: " << gradBuoyancy_.extent(0) << std::endl;
  }

 private:
  vectorReal gradKappa_;     ///< Gradient with respect to kappa.
  vectorReal gradBuoyancy_;  ///< Gradient with respect to buoyancy.
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_GRADIENT_ACOUSTIC_H_

#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_GRADIENT_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_GRADIENT_ELASTIC_H_

#include <iostream>
#include <string>

#include "gradient.h"

namespace gradient {
/**
 * @brief Holds the three elastic gradient fields (rho, lambda, mu).
 *
 * The fields are flat vectors, indexed like the model arrays they are
 * computed for. The class stores the vectors it is given (shallow copies of
 * the views), it does not allocate.
 */
class GradientElastic : public Gradient {
 public:
  static constexpr int kNumGrads = 3;  ///< Number of gradient fields.

  /**
   * @brief Wraps the three gradient fields.
   * @param[in] gradRho Gradient with respect to rho.
   * @param[in] gradLambda Gradient with respect to lambda.
   * @param[in] gradMu Gradient with respect to mu.
   */
  GradientElastic(vectorReal gradRho, vectorReal gradLambda, vectorReal gradMu)
      : gradRho_(gradRho), gradLambda_(gradLambda), gradMu_(gradMu) {}

  /** @brief Returns the number of gradient fields (3). */
  int getNumGradients() const override final { return kNumGrads; }

  /**
   * @brief Returns the name of gradient field i.
   * @param[in] i Field index: 0 = rho, 1 = lambda, 2 = mu.
   * @return "gradRho", "gradLambda" or "gradMu".
   */
  std::string getGradientName(int i) const override final {
    switch (i) {
      case 0:
        return "gradRho";
      case 1:
        return "gradLambda";
      case 2:
        return "gradMu";
      default:
        return "gradRho";
    }
  }

  // TODO use template + constexpr if when C++20 is available
  /**
   * @brief Returns gradient field i.
   * @param[in] i Field index: 0 = rho, 1 = lambda, 2 = mu.
   * @return The stored vector (shallow copy).
   */
  PROXY_HOST_DEVICE
  vectorReal getGradient(int i) const override {
    switch (i) {
      case 0:
        return gradRho_;
      case 1:
        return gradLambda_;
      case 2:
        return gradMu_;
      default:
        return gradRho_;  // make it cuda happy
    }
  }

  /** @brief Prints the size of each gradient field to stdout. */
  void print() const override {
    std::cout << "Grad Rho size: " << gradRho_.extent(0) << std::endl;
    std::cout << "Grad Lambda size: " << gradLambda_.extent(0) << std::endl;
    std::cout << "Grad Mu size: " << gradMu_.extent(0) << std::endl;
  }

 private:
  vectorReal gradRho_;     ///< Gradient with respect to rho.
  vectorReal gradLambda_;  ///< Gradient with respect to lambda.
  vectorReal gradMu_;      ///< Gradient with respect to mu.
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_GRADIENT_ELASTIC_H_

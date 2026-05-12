#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_GRADIENT_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_GRADIENT_ELASTIC_H_

#include <iostream>
#include <string>

#include "gradient.h"

namespace gradient {
/**
 * @brief Elastic gradient data structure.
 * Arrays are kept flat for easy cpp-python-fortran interop.
 */
class GradientElastic : public Gradient {
 public:
  static constexpr int kNumGrads = 3;

  GradientElastic(vectorReal gradRho, vectorReal gradLambda, vectorReal gradMu)
      : gradRho_(gradRho), gradLambda_(gradLambda), gradMu_(gradMu) {}

  int getNumGradients() const override final { return kNumGrads; }

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

  void print() const override {
    std::cout << "Grad Rho size: " << gradRho_.extent(0) << std::endl;
    std::cout << "Grad Lambda size: " << gradLambda_.extent(0) << std::endl;
    std::cout << "Grad Mu size: " << gradMu_.extent(0) << std::endl;
  }

 private:
  vectorReal gradRho_;     ///< Gradient Rho field
  vectorReal gradLambda_;  ///< Gradient Lambda field
  vectorReal gradMu_;      ///< Gradient Mu field
};
}  // namespace gradient
#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_GRADIENT_ELASTIC_H_

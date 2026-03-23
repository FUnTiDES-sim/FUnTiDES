#ifndef FUNTIDES_INVERSION_API_INCLUDE_GRADIENT_DATA_H_
#define FUNTIDES_INVERSION_API_INCLUDE_GRADIENT_DATA_H_

#include <memory>

#include "gradients.h"
#include "gradients_acoustic.h"
#include "gradients_elastic.h"

namespace inversion
{

/**
 * @brief Container for model parameter gradients.
 *
 * GradientData holds the computed gradients independently from the solver.
 * This allows gradients to be computed, stored, and accessed separately from
 * wave propagation data structures.
 *
 * Unlike SEMsolverData which contains wavefields and RHS, GradientData
 * only stores the gradient arrays themselves.
 */
struct GradientData
{
  std::shared_ptr<Gradients> gradients;

  GradientData() = default;

  /**
   * @brief Create acoustic gradient data with pre-allocated views.
   */
  GradientData(const VECTOR_REAL_VIEW& grad_kappa,
               const VECTOR_REAL_VIEW& grad_buoyancy)
      : gradients(
            std::make_shared<GradientsAcoustic>(grad_kappa, grad_buoyancy))
  {
  }

  /**
   * @brief Create elastic gradient data with pre-allocated views.
   */
  GradientData(const VECTOR_REAL_VIEW& grad_rho,
               const VECTOR_REAL_VIEW& grad_lambda,
               const VECTOR_REAL_VIEW& grad_mu)
      : gradients(std::make_shared<GradientsElastic>(grad_rho, grad_lambda,
                                                      grad_mu))
  {
  }

  /**
   * @brief Print gradient data information.
   */
  void print() const
  {
    if (gradients)
    {
      gradients->print();
    }
    else
    {
      std::cout << "GradientData: No gradients allocated\n";
    }
  }
};

}  // namespace inversion

#endif  // FUNTIDES_INVERSION_API_INCLUDE_GRADIENT_DATA_H_

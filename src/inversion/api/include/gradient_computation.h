#ifndef FUNTIDES_INVERSION_API_INCLUDE_GRADIENT_COMPUTATION_H_
#define FUNTIDES_INVERSION_API_INCLUDE_GRADIENT_COMPUTATION_H_

#include "common_macros.h"
#include "gradient_data.h"

namespace inversion
{

/**
 * @brief Abstract interface for gradient computation.
 *
 * GradientComputation provides a physics-agnostic interface for computing
 * model parameter gradients from forward and adjoint wavefields.
 *
 * This class is independent from the Solver - gradients are computed as a
 * separate post-processing step, not as part of wave propagation.
 *
 * Usage:
 *   auto grad_computer = createGradientComputation(...);
 *   GradientData grad_data;
 *   grad_computer->compute(forward_field, adjoint_fields, grad_data);
 *   auto grad_kappa = grad_data.gradients->getCurrentGrads(0);
 */
class GradientComputation
{
 public:
  virtual ~GradientComputation() = default;

  /**
   * @brief Compute gradients from forward and adjoint wavefields.
   *
   * Computes model parameter gradients (Kappa, Buoyancy for acoustic;
   * Rho, Lambda, Mu for elastic) based on forward and backward propagations.
   *
   * @param forward_field        Forward wavefield (pressure or displacement)
   * @param adjoint_fields       Adjoint wavefields (e.g., [adjoint_dt, adjoint])
   * @param grad_data            Output structure with computed gradients
   */
  virtual void compute(const VECTOR_REAL_VIEW& forward_field,
                       const std::vector<VECTOR_REAL_VIEW>& adjoint_fields,
                       GradientData& grad_data) = 0;

  /**
   * @brief Get reference to stored gradients object.
   * @return Pointer to Gradients base class instance
   */
  virtual Gradients* getGradients() = 0;

  /**
   * @brief Get polynomial order of this computation.
   * @return polynomial order (ORDER template parameter)
   */
  virtual int getOrder() const = 0;

  /**
   * @brief Check if model is discretized at nodes or element centers.
   * @return true if node-based, false if element-based
   */
  virtual bool isModelOnNodes() const = 0;

  /**
   * @brief Print information about this gradient computation.
   */
  virtual void print() const = 0;
};

}  // namespace inversion

#endif  // FUNTIDES_INVERSION_API_INCLUDE_GRADIENT_COMPUTATION_H_

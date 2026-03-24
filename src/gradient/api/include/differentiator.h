#ifndef FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_

#include "model.h"

namespace gradient
{

/**
 * @brief Abstract interface for gradient computation.
 *
 * Differentiator provides a physics-agnostic interface for computing
 * model parameter gradients from forward and adjoint wavefields.
 *
 * This class is independent from the Solver - gradients are computed as a
 * separate post-processing step, not as part of wave propagation.
 *
 * Usage:
 *   auto differentiator = createDifferentiator(...);
 *   GradientData data;
 *   differentiator->compute(mesh, data);
 *   auto grad_kappa = data.getGradient(0);
 */
class Differentiator
{
 public:
  virtual ~Differentiator() = default;

  struct DataStruct
  {
    virtual ~DataStruct() = default;
    virtual void print() const = 0;
  };

  /**
   * @brief Compute gradients from forward and adjoint wavefields.
   *
   * Computes model parameter gradients (Kappa, Buoyancy for acoustic;
   * Rho, Lambda, Mu for elastic) based on forward and backward propagations.
   *
   * @param mesh  The model mesh containing geometry and parameter info
   * @param data  Structure with required input fields and computed gradients
   */
  virtual void compute(model::ModelApi<float, int>& mesh,
                       DataStruct& data) = 0;

  /**
   * @brief Get polynomial order of this computation.
   * @return polynomial order
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

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_

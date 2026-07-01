#ifndef FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_

#include "model.h"

namespace gradient {

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
 *   float dt = 0.001;
 *   differentiator->compute(mesh, data, dt);
 *   auto grad_kappa = data.getGradient(0);
 */
class Differentiator {
 public:
  virtual ~Differentiator() = default;

  struct DataStruct {
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
   * @param dt    Time step size (for acoustic and elastic time derivatives)
   */
  virtual void compute(model::ModelApi<float, int>& mesh, DataStruct& data, float dt) const = 0;

  /**
   * @brief Initialize the geometric mass matrix (nodal volumes without model factors).
   *
   * Assembles Omega_I = sum_{e in I} w_I^e |J_I^e| without applying model
   * factors (velocities/densities). It must be called once before compute():
   * node-based gradients are normalized by it, and it is exposed via
   * getGeometricMassMatrix() for FWI preconditioning K(x_I) = G_I / Omega_I.
   *
   * @param mesh The model mesh containing geometry info
   */
  virtual void initGeometricMassMatrix(model::ModelApi<float, int>& mesh) = 0;

  /**
   * @brief Get the geometric mass matrix (nodal volumes without velocities).
   *
   * The geometric mass matrix Omega_I = sum_{e in I} w_I^e |J_I^e| is computed
   * without applying model factors (velocities/densities). It is used for FWI
   * preconditioning: K(x_I) = G_I / Omega_I.
   *
   * Call initGeometricMassMatrix() once before use.
   *
   * @return Reference to the geometric mass matrix vector
   */
  virtual vectorReal& getGeometricMassMatrix() = 0;

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

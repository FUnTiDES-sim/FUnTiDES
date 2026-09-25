#ifndef FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_
#define FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_

#include "model.h"

namespace gradient {

/**
 * @brief Computes model parameter gradients from a forward and an adjoint
 * wavefield, as a post-processing step separate from wave propagation.
 *
 * An instance is bound at creation to one physics, polynomial order, mesh type
 * and model location (nodes or elements). The caller calls compute() once per
 * time step; each call adds that step's contribution to the gradient arrays
 * held by the DataStruct, so the caller zeroes them before the first step.
 */
class Differentiator {
 public:
  virtual ~Differentiator() = default;

  /**
   * @brief Physics-specific input of compute(): forward and adjoint wavefield
   * views plus the gradient arrays to accumulate into.
   *
   * compute() only accepts the data type of its own physics.
   */
  struct DataStruct {
    virtual ~DataStruct() = default;
    /** @brief Prints a description of the held arrays to standard output. */
    virtual void print() const = 0;
  };

  /**
   * @brief Adds the gradient contribution of one time step to the gradient
   * arrays of @p data.
   *
   * Gradient arrays hold one value per global node if isModelOnNodes(), one
   * value per element otherwise; wavefield arrays hold one value per global
   * node.
   *
   * @param[in] mesh Mesh of the concrete type and polynomial order this
   *            differentiator was created for.
   * @param[in,out] data Wavefield views (read) and gradient arrays
   *                (accumulated), of the data type of this physics.
   * @param[in] dt Time step between consecutive snapshots of the adjoint view.
   * @throws std::bad_cast if @p mesh or @p data is not of the expected type.
   */
  virtual void compute(model::ModelApi<float, int>& mesh, DataStruct& data, float dt) const = 0;

  /**
   * @brief Assembles the geometric mass matrix of @p mesh (nodal volumes
   * without model factors).
   *
   * Computes Omega_I = sum over the elements e containing node I of
   * w_I^e |J_I^e|, the quadrature weight times the Jacobian determinant.
   * Call it once per mesh, before compute() and getGeometricMassMatrix().
   *
   * @param[in] mesh Mesh of the concrete type and polynomial order this
   *            differentiator was created for.
   * @throws std::bad_cast if @p mesh is not of the expected type.
   */
  virtual void initGeometricMassMatrix(model::ModelApi<float, int>& mesh) = 0;

  /**
   * @brief Returns the geometric mass matrix built by
   * initGeometricMassMatrix().
   *
   * Intended for FWI preconditioning: K(x_I) = G_I / Omega_I, where G_I is a
   * node-based gradient.
   *
   * @return Handle to the internal array, size numNodes, indexed by global
   * node; empty before initGeometricMassMatrix().
   */
  virtual vectorReal& getGeometricMassMatrix() = 0;

  /** @brief Returns the polynomial order of the elements. */
  virtual int getOrder() const = 0;

  /**
   * @brief Tells where the model, hence the gradient, is discretized.
   * @return true for one value per global node, false for one per element.
   */
  virtual bool isModelOnNodes() const = 0;

  /** @brief Prints the configuration to standard output. */
  virtual void print() const = 0;
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_API_INCLUDE_DIFFERENTIATOR_H_

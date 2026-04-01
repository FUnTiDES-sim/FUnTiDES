#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_

#include <iostream>

#include "differentiator.h"
#include "model.h"

namespace gradient
{

/**
 * @brief Elastic gradient computation for independent use.
 *
 * Computes model parameter gradients (grad_rho, grad_lambda, grad_mu) from
 * elastic forward and adjoint displacement wavefields. Completely independent
 * from the Solver.
 *
 * Features:
 * - Supports both node-based and element-based model discretization
 * - Uses standard SEM assembly with mass and stiffness matrices
 *
 * Template Parameters:
 *   ORDER                 - Polynomial order (1, 2, 3, ...)
 *   INTEGRAL_TYPE         - Integration kernel (e.g., makutu)
 *   MESH_TYPE             - Mesh topology (e.g., Cartesian)
 *   IS_MODEL_ON_NODES     - Model discretization (true=nodes, false=elements)
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class DifferentiatorElastic : public Differentiator
{
 public:
  ~DifferentiatorElastic() override = default;

  /**
   * @brief Compute elastic gradients (Rho, Lambda, Mu).
   *
   * Computes:
   *   grad_rho    = ∑_elements ∑_quadrature (u_fwd · u_adj) * mass_term
   *   grad_lambda = ∑_elements ∑_stiffness stress_strain_contraction * 0.5
   *   grad_mu     = ∑_elements ∑_stiffness stress_strain_contraction
   */
  void compute(model::ModelApi<float, int>& mesh, DataStruct& data,
               float dt) const override
  {
    throw std::runtime_error(
        "Elastic gradient computation not implemented yet");
  }

  int getOrder() const override { return kOrder; }

  bool isModelOnNodes() const override { return kIsModelOnNodes; }

  void print() const override
  {
    std::cout << "DifferentiatorElastic<ORDER=" << kOrder
              << ", INTEGRAL_TYPE=" << typeid(INTEGRAL_TYPE).name()
              << ", MESH_TYPE=" << typeid(MESH_TYPE).name()
              << ", IS_MODEL_ON_NODES=" << (kIsModelOnNodes ? "true" : "false")
              << ">\n";
  }

 private:
  static constexpr int kOrder = ORDER;
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;
  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_

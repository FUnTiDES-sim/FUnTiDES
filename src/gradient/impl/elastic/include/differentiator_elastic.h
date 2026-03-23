#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_

#include <vector>

#include "common_macros.h"
#include "data_type.h"
#include "mesh.h"
#include "model_discretization_interface.h"
#include "differentiator.h"
#include "gradient_data.h"

namespace gradient
{

/**
 * @brief Elastic gradient computation for independent use.
 *
 * Computes model parameter gradients (grad_rho, grad_lambda, grad_mu) from elastic
 * forward and adjoint displacement wavefields. Completely independent from the Solver.
 *
 * Features:
 * - Supports both node-based and element-based model discretization
 * - Uses standard SEM assembly with mass and stiffness matrices
 * - Handles 3D vector displacement fields
 * - Thread-safe accumulation with ATOMICADD
 *
 * Template Parameters:
 *   ORDER                 - Polynomial order (1, 2, 3, ...)
 *   INTEGRAL_TYPE         - Integration kernel (e.g., makutu)
 *   MESH_TYPE             - Mesh topology (e.g., CartesianStructBuilder)
 *   IS_MODEL_ON_NODES     - Model discretization (true=nodes, false=elements)
 *
 * Usage:
 *   ElasticDifferentiator<2, makutu, CartesianStructBuilder, true> gc(...);
 *   GradientData grad_data(grad_rho_view, grad_lambda_view, grad_mu_view);
 *   gc.compute(forward_displacement, adjoint_displacement, grad_data);
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class ElasticDifferentiator : public Differentiator
{
 public:
  static constexpr int kOrder = ORDER;
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;
  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  /**
   * @brief Constructor for elastic gradient computation.
   *
   * @param mesh Mesh object containing topology and reference element info
   */
  explicit ElasticDifferentiator(const MESH_TYPE& mesh)
      : m_mesh(mesh)
  {
  }

  ~ElasticDifferentiator() override = default;

  /**
   * @brief Compute elastic gradients (Rho, Lambda, Mu).
   *
   * Computes:
   *   grad_rho    = ∑_elements ∑_quadrature (u_fwd · u_adj) * mass_term
   *   grad_lambda = ∑_elements ∑_stiffness stress_strain_contraction * 0.5
   *   grad_mu     = ∑_elements ∑_stiffness stress_strain_contraction
   *
   * @param forward_field       Displacement wavefield from forward propagation
   * @param adjoint_fields      [0]=adjoint displacement (u_adj)
   * @param grad_data           Output structure containing gradient arrays
   *
   * @throws std::runtime_error if adjoint_fields size != 1
   * @throws std::runtime_error if gradient types mismatch
   * @throws std::runtime_error if forward_field size not divisible by 3
   */
  void compute(const VECTOR_REAL_VIEW& forward_field,
               const std::vector<VECTOR_REAL_VIEW>& adjoint_fields,
               Differentiator::DataStruct& grad_data) override;

  Gradients* getGradients() override { return grad_data_ptr_; }

  int getOrder() const override { return kOrder; }

  bool isModelOnNodes() const override { return kIsModelOnNodes; }

  void print() const override
  {
    std::cout << "ElasticDifferentiator<ORDER=" << kOrder
              << ", IS_MODEL_ON_NODES=" << (kIsModelOnNodes ? "true" : "false")
              << ">\n";
  }

 private:
  MESH_TYPE m_mesh;
  Gradients* grad_data_ptr_ = nullptr;  // Non-owning pointer for getGradients()

  /**
   * @brief Kernel: compute mass term contribution for grad_rho (velocity dot product).
   */
  void computeGradRhoMassTerm(
      const VECTOR_REAL_VIEW& forward_field,
      const VECTOR_REAL_VIEW& adjoint_field,
      const typename INTEGRAL_TYPE::TransformType& transformData,
      int elementNumber, float* localGradRho);

  /**
   * @brief Kernel: compute stiffness term for grad_lambda and grad_mu.
   */
  void computeGradLambdaMuStiffnessTerm(
      const VECTOR_REAL_VIEW& forward_field,
      const VECTOR_REAL_VIEW& adjoint_field,
      const typename INTEGRAL_TYPE::TransformType& transformData,
      int elementNumber, float* localGradLambda, float* localGradMu);

  /**
   * @brief Accumulate element-local gradients to global arrays.
   */
  void accumulateGradients(int elementNumber, const float* localGradRho,
                           const float* localGradLambda,
                           const float* localGradMu, VECTOR_REAL_VIEW gradRho,
                           VECTOR_REAL_VIEW gradLambda,
                           VECTOR_REAL_VIEW gradMu);
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_

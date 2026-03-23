#ifndef FUNTIDES_INVERSION_IMPL_ACOUSTIC_INCLUDE_GRADIENT_COMPUTATION_ACOUSTIC_H_
#define FUNTIDES_INVERSION_IMPL_ACOUSTIC_INCLUDE_GRADIENT_COMPUTATION_ACOUSTIC_H_

#include <vector>

#include "common_macros.h"
#include "data_type.h"
#include "mesh.h"
#include "model_discretization_interface.h"
#include "gradient_computation.h"
#include "gradient_data.h"

namespace inversion
{

/**
 * @brief Acoustic gradient computation for independent use.
 *
 * Computes model parameter gradients (grad_kappa, grad_buoyancy) from acoustic
 * forward and adjoint wavefields. Completely independent from the Solver class.
 *
 * Features:
 * - Supports both node-based and element-based model discretization
 * - Uses standard SEM assembly with mass and stiffness matrices
 * - Thread-safe accumulation with ATOMICADD
 *
 * Template Parameters:
 *   ORDER                 - Polynomial order (1, 2, 3, ...)
 *   INTEGRAL_TYPE         - Integration kernel (e.g., makutu)
 *   MESH_TYPE             - Mesh topology (e.g., CartesianStructBuilder)
 *   IS_MODEL_ON_NODES     - Model discretization (true=nodes, false=elements)
 *
 * Usage:
 *   AcousticGradientComputation<2, makutu, CartesianStructBuilder, true> gc(...);
 *   GradientData grad_data(grad_kappa_view, grad_buoyancy_view);
 *   gc.compute(forward, adjoint_dt, grad_data);
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class AcousticGradientComputation : public GradientComputation
{
 public:
  static constexpr int kOrder = ORDER;
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;
  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  /**
   * @brief Constructor for acoustic gradient computation.
   *
   * @param mesh Mesh object containing topology and reference element info
   */
  explicit AcousticGradientComputation(const MESH_TYPE& mesh)
      : m_mesh(mesh)
  {
  }

  ~AcousticGradientComputation() override = default;

  /**
   * @brief Compute acoustic gradients (Kappa, Buoyancy).
   *
   * Computes:
   *   grad_kappa   = ∑_elements ∑_quadrature q_dt² * p * mass_term
   *   grad_buoyancy = ∑_elements ∑_stiffness stiffness_term * q * p
   *
   * @param forward_field       Pressure wavefield from forward propagation
   * @param adjoint_fields      [0]=adjoint_dt (∂²q/∂t²), [1]=adjoint (q)
   * @param grad_data           Output structure containing gradient arrays
   *
   * @throws std::runtime_error if adjoint_fields size != 2
   * @throws std::runtime_error if gradient types mismatch
   */
  void compute(const VECTOR_REAL_VIEW& forward_field,
               const std::vector<VECTOR_REAL_VIEW>& adjoint_fields,
               GradientData& grad_data) override;

  Gradients* getGradients() override { return grad_data_ptr_; }

  int getOrder() const override { return kOrder; }

  bool isModelOnNodes() const override { return kIsModelOnNodes; }

  void print() const override
  {
    std::cout << "AcousticGradientComputation<ORDER=" << kOrder
              << ", IS_MODEL_ON_NODES=" << (kIsModelOnNodes ? "true" : "false")
              << ">\n";
  }

 private:
  MESH_TYPE m_mesh;
  Gradients* grad_data_ptr_ = nullptr;  // Non-owning pointer for getGradients()

  /**
   * @brief Kernel: compute mass term contribution for grad_kappa.
   */
  void computeGradKappaMassTerm(
      const VECTOR_REAL_VIEW& forward_field,
      const VECTOR_REAL_VIEW& adjoint_dt_field,
      const typename INTEGRAL_TYPE::TransformType& transformData,
      int elementNumber, float* localGradKappa);

  /**
   * @brief Kernel: compute stiffness term contribution for grad_buoyancy.
   */
  void computeGradBuoyancyStiffnessTerm(
      const VECTOR_REAL_VIEW& forward_field,
      const VECTOR_REAL_VIEW& adjoint_field,
      const typename INTEGRAL_TYPE::TransformType& transformData,
      int elementNumber, float* localGradBuoyancy);

  /**
   * @brief Accumulate element-local gradients to global arrays.
   */
  void accumulateGradients(int elementNumber, const float* localGradKappa,
                           const float* localGradBuoyancy,
                           VECTOR_REAL_VIEW gradKappa,
                           VECTOR_REAL_VIEW gradBuoyancy);
};

}  // namespace inversion

#endif  // FUNTIDES_INVERSION_IMPL_ACOUSTIC_INCLUDE_GRADIENT_COMPUTATION_ACOUSTIC_H_

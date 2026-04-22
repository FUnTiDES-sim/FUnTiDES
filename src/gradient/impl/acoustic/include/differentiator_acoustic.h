#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_

#include <iostream>

#include "differentiator.h"
#include "differentiator_data_acoustic.h"
#include "model.h"

namespace gradient
{

/**
 * @brief Acoustic gradient computation for independent use.
 * Computes model parameter gradients (grad_kappa, grad_buoyancy) from acoustic
 * forward and adjoint wavefields. Completely independent from the Solver.
 * Features:
 * - Supports both node-based and element-based model discretization
 * - Uses standard SEM assembly with mass and stiffness matrices
 * Template Parameters:
 *   ORDER                 - Polynomial order (1, 2, 3, ...)
 *   INTEGRAL_TYPE         - Integration kernel (e.g., makutu)
 *   MESH_TYPE             - Mesh topology (e.g., Cartesian)
 *   IS_MODEL_ON_NODES     - Model discretization (true=nodes, false=elements)
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class DifferentiatorAcoustic : public Differentiator
{
 public:
  static constexpr int kOrder = ORDER;
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;
  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  KOKKOS_DEFAULTED_FUNCTION ~DifferentiatorAcoustic() override = default;

  /**
   * @brief Compute acoustic gradients (Kappa, Buoyancy).
   *
   * Computes:
   *   grad_kappa   = ∑_elements ∑_quadrature q_dt² * p * mass_term
   *   grad_buoyancy = ∑_elements ∑_stiffness stiffness_term * q * p
   */
  void compute(model::ModelApi<float, int>& mesh, DataStruct& data,
               float dt) const override;

  int getOrder() const override;
  bool isModelOnNodes() const override;
  void print() const override;

  /**
   * @brief Each element writes to a unique index — no atomic add required.
   *
   * Computes qdt2 = (qnPrevPrev - 2*qnPrev + qn) / dt² on the fly.
   */
  void computeOnElements(MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn,
                         VECTOR_REAL_VIEW const qn,
                         VECTOR_REAL_VIEW const qnPrev,
                         VECTOR_REAL_VIEW const qnPrevPrev,
                         VECTOR_REAL_VIEW const gradKappa,
                         VECTOR_REAL_VIEW const gradBuoyancy) const;

  /**
   * @brief Multiple elements share boundary nodes — ATOMICADD required.
   *
   * Computes qdt2 = (qnPrevPrev - 2*qnPrev + qn) / dt² on the fly.
   * Gradients are normalized by the mass matrix diagonal to account for
   * multiple elements sharing nodes at the domain interior.
   */
  void computeOnNodes(MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn,
                      VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qnPrev,
                      VECTOR_REAL_VIEW const qnPrevPrev,
                      VECTOR_REAL_VIEW const gradKappa,
                      VECTOR_REAL_VIEW const gradBuoyancy) const;
};

}  // namespace gradient

// ============================================================================
// EXTERN TEMPLATES (avoid template bloat)
// ============================================================================
#include "Integrals.h"
#include "model_struct.h"
#include "model_unstruct.h"

#define DECLARE_EXTERN_DIFF_ACOUSTIC(ORDER, MESH_TYPE)                         \
  extern template class gradient::DifferentiatorAcoustic<                      \
      ORDER, typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, \
      MESH_TYPE, true>;                                                        \
  extern template class gradient::DifferentiatorAcoustic<                      \
      ORDER, typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, \
      MESH_TYPE, false>;

#define DECLARE_EXTERN_DIFF_ACOUSTIC_ALL_ORDERS(MESH_TYPE_MACRO) \
  DECLARE_EXTERN_DIFF_ACOUSTIC(1, MESH_TYPE_MACRO(1))            \
  DECLARE_EXTERN_DIFF_ACOUSTIC(2, MESH_TYPE_MACRO(2))            \
  DECLARE_EXTERN_DIFF_ACOUSTIC(3, MESH_TYPE_MACRO(3))

#define STRUCT_MESH_TYPE(ORDER) model::ModelStruct<float, int, ORDER>
#define UNSTRUCT_MESH_TYPE(ORDER) model::ModelUnstruct<float, int>

DECLARE_EXTERN_DIFF_ACOUSTIC_ALL_ORDERS(STRUCT_MESH_TYPE)
DECLARE_EXTERN_DIFF_ACOUSTIC_ALL_ORDERS(UNSTRUCT_MESH_TYPE)

#undef UNSTRUCT_MESH_TYPE
#undef STRUCT_MESH_TYPE
#undef DECLARE_EXTERN_DIFF_ACOUSTIC_ALL_ORDERS
#undef DECLARE_EXTERN_DIFF_ACOUSTIC

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_

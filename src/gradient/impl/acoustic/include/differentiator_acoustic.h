#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_

#include <iostream>

#include "differentiator.h"
#include "differentiator_data_acoustic.h"
#include "model.h"

namespace gradient {

/**
 * @brief Computes acoustic model gradients (grad_kappa, grad_buoyancy) from forward and
 * adjoint pressure wavefields.
 *
 * Works on its own, without a Solver. The model parameters may live on nodes or on elements.
 *
 * @tparam ORDER             Polynomial order of the elements.
 * @tparam INTEGRAL_TYPE     Element integration kernel.
 * @tparam MESH_TYPE         Mesh/model type the gradients are computed on.
 * @tparam IS_MODEL_ON_NODES True if the model is stored per node, false if per element.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
class DifferentiatorAcoustic : public Differentiator {
 public:
  static constexpr int kOrder = ORDER;                    ///< Polynomial order.
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;  ///< True if the model is stored per node.
  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);  ///< Nodes per element.

  KOKKOS_DEFAULTED_FUNCTION ~DifferentiatorAcoustic() override = default;

  /**
   * @brief Computes the acoustic gradients (kappa, buoyancy).
   *
   * grad_kappa is built from the second time derivative of the adjoint field times the
   * forward field, weighted by the mass term. grad_buoyancy is built from the product of
   * the forward and adjoint stiffness terms.
   *
   * @param[in]     mesh Mesh and model.
   * @param[in,out] data Forward and adjoint wavefield views and output gradients.
   * @param[in]     dt   Time step between the wavefield snapshots.
   */
  void compute(model::ModelApi<float, int>& mesh, DataStruct& data, float dt) const override;

  /**
   * @brief Returns the geometric mass matrix, i.e. the nodal volumes without model factors.
   *
   * Omega_I = sum over elements e containing node I of w_I^e |J_I^e|. Used to normalize
   * gradients for FWI preconditioning: K^kappa(x_I) = G^kappa_I / Omega_I.
   *
   * @return Reference to the internal vector.
   */
  vectorReal& getGeometricMassMatrix() override;

  int getOrder() const override;
  bool isModelOnNodes() const override;
  void print() const override;

  /**
   * @brief Element-based gradient kernel.
   *
   * Each element writes to its own index, so no atomic add is needed. The second time
   * derivative of q is computed on the fly as (qnPrevPrev - 2*qnPrev + qn) / dt^2.
   *
   * @param[in]  mesh            Mesh and model.
   * @param[in]  dt              Time step between snapshots.
   * @param[in]  pn              Forward field.
   * @param[in]  qn              Adjoint field at the current step.
   * @param[in]  qnPrev          Adjoint field at the previous step.
   * @param[in]  qnPrevPrev      Adjoint field two steps back.
   * @param[out] gradKappa       Gradient with respect to kappa.
   * @param[out] gradBuoyancy    Gradient with respect to buoyancy.
   */
  void computeOnElements(MESH_TYPE mesh, float dt, vectorReal const pn, vectorReal const qn, vectorReal const qnPrev,
                         vectorReal const qnPrevPrev, vectorReal const gradKappa, vectorReal const gradBuoyancy) const;

  /**
   * @brief Node-based gradient kernel.
   *
   * Elements share boundary nodes, so contributions are accumulated with ATOMICADD. The
   * second time derivative of q is computed on the fly as
   * (qnPrevPrev - 2*qnPrev + qn) / dt^2. Gradients are then normalized by the diagonal of
   * the geometric mass matrix.
   *
   * @param[in]  mesh            Mesh and model.
   * @param[in]  dt              Time step between snapshots.
   * @param[in]  pn              Forward field.
   * @param[in]  qn              Adjoint field at the current step.
   * @param[in]  qnPrev          Adjoint field at the previous step.
   * @param[in]  qnPrevPrev      Adjoint field two steps back.
   * @param[out] gradKappa       Gradient with respect to kappa.
   * @param[out] gradBuoyancy    Gradient with respect to buoyancy.
   */
  void computeOnNodes(MESH_TYPE mesh, float dt, vectorReal const pn, vectorReal const qn, vectorReal const qnPrev,
                      vectorReal const qnPrevPrev, vectorReal const gradKappa, vectorReal const gradBuoyancy) const;

  /**
   * @brief Builds the geometric mass matrix (nodal volumes without model factors).
   *
   * The result is exposed through getGeometricMassMatrix().
   *
   * @param[in] mesh The computational mesh.
   * @note Public because CUDA device lambdas in Kokkos cannot be defined in private members.
   */
  void initGeometricMassMatrix(model::ModelApi<float, int>& mesh) override;

 private:
  vectorReal geometricMassMatrix_;  ///< Nodal volumes without model factors.
};

}  // namespace gradient

// Explicit instantiations for orders 1 to 3 live in a source file; declare them extern
// here so that includers do not instantiate the class again.
#include "Integrals.h"
#include "model_struct.h"
#include "model_unstruct.h"

#define DECLARE_EXTERN_DIFF_ACOUSTIC(ORDER, MESH_TYPE)                                           \
  extern template class gradient::DifferentiatorAcoustic<                                        \
      ORDER, typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, MESH_TYPE, true>; \
  extern template class gradient::DifferentiatorAcoustic<                                        \
      ORDER, typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, MESH_TYPE, false>;

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

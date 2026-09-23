#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_

#include <Kokkos_Macros.hpp>
#include <iostream>

#include "differentiator.h"
#include "differentiator_data_elastic.h"
#include "model.h"

namespace gradient {

/**
 * @brief Elastic model gradients (rho, lambda, mu) from forward and adjoint displacement wavefields.
 *
 * Runs independently of the Solver. For isotropic media the three sensitivities are
 * (u = forward displacement, u* = adjoint displacement, u_dd* = adjoint second time derivative):
 *
 *   grad_rho    = - sum_t sum_e integral u_dd* . u dOmega
 *   grad_lambda = - sum_t sum_e integral div(u*) div(u) dOmega
 *   grad_mu     = - sum_t sum_e integral 2 eps(u*) : eps(u) dOmega
 *
 * For TTI media the strain interaction uses the full 6x6 Voigt elasticity tensor.
 * The model can be discretized on nodes or on elements (see IS_MODEL_ON_NODES).
 *
 * @tparam ORDER              Polynomial order of the spectral elements.
 * @tparam INTEGRAL_TYPE      Integration kernel class providing the element integrals.
 * @tparam MESH_TYPE          Mesh/model type traversed by the kernels.
 * @tparam IS_MODEL_ON_NODES  true if the model parameters and gradients live on nodes,
 *                            false if they live on elements.
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
class DifferentiatorElastic : public Differentiator {
 public:
  static constexpr int kOrder = ORDER;                     ///< Polynomial order.
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;  ///< Model on nodes (true) or elements (false).
  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);  ///< Nodes per hexahedron.

  ~DifferentiatorElastic() override = default;

  /**
   * @brief Compute the elastic gradients (rho, lambda, mu) from the wavefield data.
   *
   * @param[in,out] mesh  Mesh and model.
   * @param[in,out] data  Forward and adjoint wavefield views and output gradients.
   * @param[in] dt        Time step.
   */
  void compute(model::ModelApi<float, int>& mesh, DataStruct& data, float dt) const override;

  /**
   * @brief Geometric mass matrix used to normalize gradients in FWI.
   *
   * Nodal volumes Omega_I = sum over elements e containing I of w_I^e |J_I^e|, without model factors.
   *
   * @return Reference to the vector, size number of global nodes.
   *         @todo VERIFY: is it filled only after initGeometricMassMatrix()?
   */
  vectorReal& getGeometricMassMatrix() override;

  int getOrder() const override;
  bool isModelOnNodes() const override;
  void print() const override;

  /**
   * @brief Displacement gradient at a quadrature point of one element.
   *
   * Computes grad[component][spatial] = d u_component / d x_spatial.
   *
   * @param[in] qa,qb,qc  Quadrature indices in each reference direction.
   * @param[in] J         Inverse Jacobian J^{-1}[ref_dir][phys_dir].
   * @param[in] localUx,localUy,localUz  x, y and z displacement components of the element,
   *            each indexed by local node.
   * @param[out] grad     Gradient tensor grad[component][spatial].
   */
  KOKKOS_INLINE_FUNCTION
  static void computeDisplacementGradient(int qa, int qb, int qc, float const (&J)[3][3], float const* localUx,
                                          float const* localUy, float const* localUz, float (&grad)[3][3]);

  /**
   * @brief Gradient kernel for a model discretized on elements.
   *
   * @param[in] mesh  Mesh and model.
   * @param[in] dt    Time step.
   * @param[in] ux_fwd,uy_fwd,uz_fwd  Forward displacement components.
   * @param[in] ux_adj,uy_adj,uz_adj  Adjoint displacement components.
   * @param[in] ux_dt2,uy_dt2,uz_dt2  @todo VERIFY: second time derivative of the adjoint
   *            displacement, precomputed by the caller?
   * @param[out] gradRho,gradLambda,gradMu  Gradients, one value per element.
   *            @todo VERIFY: accumulated into or overwritten?
   */
  void computeOnElements(MESH_TYPE mesh, float dt, vectorReal const ux_fwd, vectorReal const uy_fwd,
                         vectorReal const uz_fwd, vectorReal const ux_adj, vectorReal const uy_adj,
                         vectorReal const uz_adj, vectorReal const ux_dt2, vectorReal const uy_dt2,
                         vectorReal const uz_dt2, vectorReal const gradRho, vectorReal const gradLambda,
                         vectorReal const gradMu) const;

  /**
   * @brief Gradient kernel for a model discretized on nodes.
   *
   * Contributions of elements sharing a node are summed with atomic adds.
   * Parameters are the same as computeOnElements(); the gradients have one value per node.
   */
  void computeOnNodes(MESH_TYPE mesh, float dt, vectorReal const ux_fwd, vectorReal const uy_fwd,
                      vectorReal const uz_fwd, vectorReal const ux_adj, vectorReal const uy_adj,
                      vectorReal const uz_adj, vectorReal const ux_dt2, vectorReal const uy_dt2,
                      vectorReal const uz_dt2, vectorReal const gradRho, vectorReal const gradLambda,
                      vectorReal const gradMu) const;

  /**
   * @brief Build the geometric mass matrix (nodal volumes without model factors).
   *
   * The result is read through getGeometricMassMatrix().
   *
   * @param[in,out] mesh  Mesh and model.
   * @note Public because Kokkos CUDA device lambdas cannot be defined in private members.
   */
  void initGeometricMassMatrix(model::ModelApi<float, int>& mesh) override;

 private:
  vectorReal geometricMassMatrix_;  ///< Nodal volumes, see getGeometricMassMatrix().
};

}  // namespace gradient

// Explicit instantiations are compiled elsewhere; these declarations avoid re-instantiation.
#include "Integrals.h"
#include "model_struct.h"
#include "model_unstruct.h"

#define DECLARE_EXTERN_DIFF_ELASTIC(ORDER, MESH_TYPE)                                            \
  extern template class gradient::DifferentiatorElastic<                                         \
      ORDER, typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, MESH_TYPE, true>; \
  extern template class gradient::DifferentiatorElastic<                                         \
      ORDER, typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, MESH_TYPE, false>;

#define DECLARE_EXTERN_DIFF_ELASTIC_ALL_ORDERS(MESH_TYPE_MACRO) \
  DECLARE_EXTERN_DIFF_ELASTIC(1, MESH_TYPE_MACRO(1))            \
  DECLARE_EXTERN_DIFF_ELASTIC(2, MESH_TYPE_MACRO(2))            \
  DECLARE_EXTERN_DIFF_ELASTIC(3, MESH_TYPE_MACRO(3))

#define STRUCT_MESH_TYPE(ORDER) model::ModelStruct<float, int, ORDER>
#define UNSTRUCT_MESH_TYPE(ORDER) model::ModelUnstruct<float, int>

DECLARE_EXTERN_DIFF_ELASTIC_ALL_ORDERS(STRUCT_MESH_TYPE)
DECLARE_EXTERN_DIFF_ELASTIC_ALL_ORDERS(UNSTRUCT_MESH_TYPE)

#undef UNSTRUCT_MESH_TYPE
#undef STRUCT_MESH_TYPE
#undef DECLARE_EXTERN_DIFF_ELASTIC_ALL_ORDERS
#undef DECLARE_EXTERN_DIFF_ELASTIC

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_

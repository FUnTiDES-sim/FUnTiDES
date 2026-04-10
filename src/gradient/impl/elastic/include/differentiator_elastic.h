#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_H_

#include <iostream>

#include "differentiator.h"
#include "differentiator_data_elastic.h"
#include "model.h"
#include <Kokkos_Macros.hpp>

namespace gradient
{

/**
 * @brief Elastic gradient computation for independent use.
 * Template Parameters:
 * ORDER                 - Polynomial order (1, 2, 3, ...)
 * INTEGRAL_TYPE         - Integration kernel (e.g., makutu)
 * MESH_TYPE             - Mesh topology (e.g., Cartesian)
 * IS_MODEL_ON_NODES     - Model discretization (true=nodes, false=elements)
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class DifferentiatorElastic : public Differentiator
{
 public:
  static constexpr int kOrder = ORDER;
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;
  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  ~DifferentiatorElastic() override = default;

  void compute(model::ModelApi<float, int>& mesh, DataStruct& data,
               float dt) const override;

  int getOrder() const override;
  bool isModelOnNodes() const override;
  void print() const override;

  KOKKOS_INLINE_FUNCTION
  static void computeDisplacementGradient(
      int qa, int qb, int qc, float const (&J)[3][3], float const* localUx,
      float const* localUy, float const* localUz, float (&grad)[3][3]);

  void computeOnElements(
      MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const ux_fwd,
      VECTOR_REAL_VIEW const uy_fwd, VECTOR_REAL_VIEW const uz_fwd,
      VECTOR_REAL_VIEW const ux_adj, VECTOR_REAL_VIEW const uy_adj,
      VECTOR_REAL_VIEW const uz_adj, VECTOR_REAL_VIEW const ux_dt2,
      VECTOR_REAL_VIEW const uy_dt2, VECTOR_REAL_VIEW const uz_dt2,
      VECTOR_REAL_VIEW const gradRho, VECTOR_REAL_VIEW const gradLambda,
      VECTOR_REAL_VIEW const gradMu) const;

  void computeOnNodes(
      MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const ux_fwd,
      VECTOR_REAL_VIEW const uy_fwd, VECTOR_REAL_VIEW const uz_fwd,
      VECTOR_REAL_VIEW const ux_adj, VECTOR_REAL_VIEW const uy_adj,
      VECTOR_REAL_VIEW const uz_adj, VECTOR_REAL_VIEW const ux_dt2,
      VECTOR_REAL_VIEW const uy_dt2, VECTOR_REAL_VIEW const uz_dt2,
      VECTOR_REAL_VIEW const gradRho, VECTOR_REAL_VIEW const gradLambda,
      VECTOR_REAL_VIEW const gradMu) const;
};

}  // namespace gradient

// ============================================================================
// EXTERN TEMPLATES (avoid template bloat)
// ============================================================================
#include "Integrals.h"
#include "model_struct.h"
#include "model_unstruct.h"

#define DECLARE_EXTERN_DIFF_ELASTIC(ORDER, MESH_TYPE) \
  extern template class gradient::DifferentiatorElastic<ORDER, \
      typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, \
      MESH_TYPE, true>; \
  extern template class gradient::DifferentiatorElastic<ORDER, \
      typename IntegralTypeSelector<ORDER, IntegralType::MAKUTU>::type, \
      MESH_TYPE, false>;

#define DECLARE_EXTERN_DIFF_ELASTIC_ALL_ORDERS(MESH_TYPE_MACRO) \
  DECLARE_EXTERN_DIFF_ELASTIC(1, MESH_TYPE_MACRO(1)) \
  DECLARE_EXTERN_DIFF_ELASTIC(2, MESH_TYPE_MACRO(2)) \
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

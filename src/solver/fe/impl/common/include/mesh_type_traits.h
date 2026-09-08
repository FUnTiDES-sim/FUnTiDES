#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MESH_TYPE_TRAITS_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MESH_TYPE_TRAITS_H_

namespace model {
template <typename FloatType, typename ScalarType, int Order>
class ModelStruct;
}  // namespace model

namespace solver {
namespace fe {

/**
 * @brief Tells whether a mesh type maps every element with a constant Jacobian.
 *
 * When true, the Jacobian and its determinant may be evaluated once per element
 * instead of once per quadrature point. The default is false so that any mesh
 * type not explicitly known to be affine keeps the general per-quadrature-point
 * path: a wrong true here would silently degrade the operator, a wrong false
 * only costs time.
 *
 * @tparam MESH_TYPE The model type held by the solver.
 */
template <typename MESH_TYPE>
struct HasConstantJacobian {
  static constexpr bool value = false;
};

/**
 * @brief ModelStruct is a regular Cartesian grid: axis-aligned elements of
 *   identical size, hence a Jacobian that is constant over each element.
 */
template <typename FloatType, typename ScalarType, int Order>
struct HasConstantJacobian<model::ModelStruct<FloatType, ScalarType, Order>> {
  static constexpr bool value = true;
};

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MESH_TYPE_TRAITS_H_

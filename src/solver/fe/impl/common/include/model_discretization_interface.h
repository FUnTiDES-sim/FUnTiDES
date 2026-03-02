#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MODEL_DISCRETIZATION_INTERFACE_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MODEL_DISCRETIZATION_INTERFACE_H_
#include <type_traits>

#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

namespace solver
{
namespace fe
{
/// @brief Namespace for model-discretization interface utilities
namespace model_discretization_interface
{

/**
 * @brief Enumeration of the different types of transform data structures
 */
enum class transform_types
{
  linear_transform,  /// simple linear transform struct used in makutu kernel
  invalid_transform  /// invalid transform type
};

namespace detail
{
// Helper to detect if a type has a `data` member of type float[8][3]
template <typename T, typename = void>
struct is_makutu_transform : std::false_type
{
};

template <typename T>
struct is_makutu_transform<
    T, std::enable_if_t<std::is_same_v<decltype(T::data), float[8][3]>>>
    : std::true_type
{
};
}  // namespace detail

/**
 * @brief Template struct to select the transform type enumeration value
 *        corresponding to a given transform data structure type.
 * @tparam TRANSFORM_TYPE The transform data structure type.
 */
template <typename T, typename = void>
struct transform_type_selector
{
  /// define the transform type as invalid by default
  static constexpr transform_types type = transform_types::invalid_transform;
};

// Specialization for makutu TransformType (detected by float data[8][3] member)
template <typename T>
struct transform_type_selector<
    T, std::enable_if_t<detail::is_makutu_transform<T>::value>>
{
  static constexpr transform_types type = transform_types::linear_transform;
};

/**
 * @brief Gathers the transform data for a given element from the mesh into the
 *        provided transform data structure.
 */
template <typename MESH_TYPE, typename TRANSFORM_TYPE>
static constexpr PROXY_HOST_DEVICE void gatherTransformData(
    const int& elementNumber, const MESH_TYPE& mesh,
    TRANSFORM_TYPE& transformData)
{
  using TT = std::remove_cv_t<TRANSFORM_TYPE>;

  if constexpr (transform_type_selector<TT>::type ==
                transform_types::linear_transform)
  {
    typename MESH_TYPE::IndexType elementIndex =
        mesh.elementIndex(elementNumber);

    int I = 0;
    for (int k = 0; k < 2; ++k)
    {
      for (int j = 0; j < 2; ++j)
      {
        for (int i = 0; i < 2; ++i)
        {
          typename MESH_TYPE::IndexType const vertexIndex =
              mesh.globalVertexIndex(elementIndex, i, j, k);
          mesh.vertexCoords(vertexIndex, transformData.data[I]);
          ++I;
        }
      }
    }
  }
}
}  // namespace model_discretization_interface
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_MODEL_DISCRETIZATION_INTERFACE_H_

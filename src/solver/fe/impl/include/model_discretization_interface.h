#pragma once

#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include <type_traits>

#ifdef ENABLE_Shiva
#include "shiva/geometry/mapping/LinearTransform.hpp"
#include "shiva/geometry/mapping/UniformScaling.hpp"
#endif

/// @brief Namespace for model-discretization interface utilities
namespace model_discretization_interface
{

/**
 * @brief Enumeration of the different types of transform data structures
 */
enum class transform_types
{
#ifdef ENABLE_Shiva
  shiva_linear_transform,           /// shiva linear transform
  shiva_uniform_scaling_transform,  /// shiva uniform scaling transform (single h value)
  shiva_scaling_transform,          /// shiva scaling transform (hx, hy, hz)
#endif
  linear_transform,  /// simple linear transform struct used in makutu kernel
  invalid_transform  /// invalid transform type
};

namespace detail
{
  // Helper to detect if a type has a `data` member of type float[8][3]
  template <typename T, typename = void>
  struct is_makutu_transform : std::false_type {};

  template <typename T>
  struct is_makutu_transform<T,
      std::enable_if_t<std::is_same_v<decltype(T::data), float[8][3]>>>
      : std::true_type {};
}

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
struct transform_type_selector<T,
    std::enable_if_t<detail::is_makutu_transform<T>::value>>
{
  static constexpr transform_types type = transform_types::linear_transform;
};

#ifdef ENABLE_Shiva

// Specialization for shiva::geometry::LinearTransform
template <typename REAL_TYPE, typename INTERPOLATED_SHAPE>
struct transform_type_selector<
    shiva::geometry::LinearTransform<REAL_TYPE, INTERPOLATED_SHAPE>,
    void>
{
  static constexpr transform_types type =
      transform_types::shiva_linear_transform;
};

// Specialization for shiva::geometry::UniformScaling
template <typename REAL_TYPE>
struct transform_type_selector<
    shiva::geometry::UniformScaling<REAL_TYPE, void>,
    void>
{
  static constexpr transform_types type =
      transform_types::shiva_uniform_scaling_transform;
};

#endif

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

#ifdef ENABLE_Shiva
  if constexpr (transform_type_selector<TT>::type ==
                transform_types::shiva_linear_transform)
  {
    typename MESH_TYPE::IndexType const elementIndex =
        mesh.elementIndex(elementNumber);
    typename TRANSFORM_TYPE::DataType& cellCoordData = transformData.getData();
    for (int k = 0; k < 2; ++k)
    {
      for (int j = 0; j < 2; ++j)
      {
        for (int i = 0; i < 2; ++i)
        {
          typename MESH_TYPE::IndexType const vertexIndex =
              mesh.globalVertexIndex(elementIndex, i, j, k);
          float* const coords = &cellCoordData(i, j, k, 0);
          mesh.vertexCoords(vertexIndex, coords);
        }
      }
    }
  }
  else
#endif
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

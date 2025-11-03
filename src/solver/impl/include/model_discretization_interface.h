#pragma once

#include "shiva/geometry/mapping/LinearTransform.hpp"
#include "shiva/geometry/mapping/UniformScaling.hpp"
#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

/// @brief Namespace for model-discretization interface utilities
namespace model_discretization_interface
{

/**
 * @brief Enumeration of the different types of transform data structures
 */
enum class transform_types
{
  shiva_linear_transform,           /// shiva linear transform
  shiva_uniform_scaling_transform,  /// shiva uniform scaling transform (single h value)
  shiva_scaling_transform,          /// shiva scaling transform (hx, hy, hz)
  linear_transform,                 /// simple linear transform struct used in makutu kernel
  invalid_transform                 /// invalid transform type
};

/**
 * @brief Template struct to select the transform type enumeration value
 *        corresponding to a given transform data structure type.
 * @tparam TRANSFORM_TYPE The transform data structure type.
 */
template< typename >
struct transform_type_selector
{
  /// define the transform type as invalid by default
  static constexpr transform_types type = transform_types::invalid_transform;
};

/**
 * @brief Specialization of transform_type_selector for shiva::geometry::LinearTransform
 * @tparam REAL_TYPE The floating point type used in the LinearTransform
 * @tparam INTERPOLATED_SHAPE The interpolated shape type used in the LinearTransform
 */
template< typename REAL_TYPE,
          typename INTERPOLATED_SHAPE >
struct transform_type_selector< shiva::geometry::LinearTransform< REAL_TYPE, INTERPOLATED_SHAPE > >
{
  /// define the transform type as shiva_linear_transform
  static constexpr transform_types type = transform_types::shiva_linear_transform;
};


/**
 * @brief Specialization of transform_type_selector for shiva::geometry::UniformScaling
 * @tparam REAL_TYPE The floating point type used in the UniformScaling
 */
template< typename REAL_TYPE >
struct transform_type_selector< shiva::geometry::UniformScaling< REAL_TYPE, void > >
{
  /// define the transform type as shiva_uniform_scaling_transform
  static constexpr transform_types type = transform_types::shiva_uniform_scaling_transform;
};



template<>
struct transform_type_selector< typename Q1_Hexahedron_Lagrange_GaussLobatto::TransformType >
{
  static constexpr transform_types type = transform_types::linear_transform;
};
template<>
struct transform_type_selector< typename Q2_Hexahedron_Lagrange_GaussLobatto::TransformType >
{
  static constexpr transform_types type = transform_types::linear_transform;
};
template<>
struct transform_type_selector< typename Q3_Hexahedron_Lagrange_GaussLobatto::TransformType >
{
  static constexpr transform_types type = transform_types::linear_transform;
};
template<>
struct transform_type_selector< typename Q4_Hexahedron_Lagrange_GaussLobatto::TransformType >
{
  static constexpr transform_types type = transform_types::linear_transform;
};


/**
 * @brief Gathers the transform data for a given element from the mesh into the
 *        provided transform data structure.
 * @tparam MESH_TYPE The type of the mesh (e.g., structured, unstructured).
 * @tparam TRANSFORM_TYPE The type of the transform data structure.
 * @param elementNumber The global element number.
 * @param mesh The mesh object.
 * @param transformData The transform data structure to populate.
 */
template< typename MESH_TYPE, typename TRANSFORM_TYPE >
static constexpr
PROXY_HOST_DEVICE
void
gatherTransformData( const int & elementNumber,
                     const MESH_TYPE & mesh,
                     TRANSFORM_TYPE & transformData )
{
  using TT = std::remove_cv_t<TRANSFORM_TYPE>;

  if constexpr ( transform_type_selector<TT>::type == transform_types::shiva_linear_transform )
  {
    typename MESH_TYPE::IndexType const elementIndex = mesh.elementIndex( elementNumber );
    typename TRANSFORM_TYPE::DataType & cellCoordData = transformData.getData();
    for ( int k = 0; k < 2; ++k )
    {
      for ( int j = 0; j < 2; ++j )
      {
        for ( int i = 0; i < 2; ++i )
        {
          typename MESH_TYPE::IndexType const vertexIndex = mesh.globalVertexIndex( elementIndex, i, j, k);
          float * const coords = &cellCoordData(i, j, k, 0);
          mesh.vertexCoords( vertexIndex, coords );
        }
      }
    }
  }
  else if constexpr ( transform_type_selector<TT>::type == transform_types::linear_transform )
  {
    typename MESH_TYPE::IndexType elementIndex = mesh.elementIndex( elementNumber );

    int I = 0;
    for (int k = 0; k < 2; ++k)
    {
      for (int j = 0; j < 2; ++j)
      {
        for (int i=0; i<2; ++i)
        {
          typename MESH_TYPE::IndexType const vertexIndex = mesh.globalVertexIndex( elementIndex, i, j, k);
          mesh.vertexCoords( vertexIndex, transformData.data[I] );
          ++I;
        }
      }
    }
  }
}



} // namespace model_discretization_interface
#pragma once

#include "shiva/geometry/mapping/LinearTransform.hpp"
#include "shiva/geometry/mapping/UniformScaling.hpp"
#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

namespace model_discretization_interface
{


enum class transform_types
{
  shiva_linear_transform,
  shiva_uniform_scaling_transform,
  shiva_scaling_transform,
  linear_transform,
  invalid_transform
};

template< typename >
struct transform_type_selector
{
  static constexpr transform_types type = transform_types::invalid_transform;
};

template< typename REAL_TYPE,
          typename INTERPOLATED_SHAPE >
struct transform_type_selector< shiva::geometry::LinearTransform< REAL_TYPE, INTERPOLATED_SHAPE > >
{
  static constexpr transform_types type = transform_types::shiva_linear_transform;
};


template< typename REAL_TYPE >
struct transform_type_selector< shiva::geometry::UniformScaling< REAL_TYPE, void > >
{
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
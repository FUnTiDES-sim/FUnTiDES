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
    typename TRANSFORM_TYPE::DataType & cellCoordData = transformData.getData();
    for ( int k = 0; k < 2; ++k )
    {
      for ( int j = 0; j < 2; ++j )
      {
        for ( int i = 0; i < 2; ++i )
        {
          int const nodeIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
          cellCoordData( i, j, k, 0 ) = mesh.nodeCoord(nodeIdx, 0);
          cellCoordData( i, j, k, 1 ) = mesh.nodeCoord(nodeIdx, 1);
          cellCoordData( i, j, k, 2 ) = mesh.nodeCoord(nodeIdx, 2);
        }
      }
    }
  }
  else if constexpr ( transform_type_selector<TT>::type == transform_types::linear_transform )
  {
    int I = 0;
    for (int k = 0; k < 2; ++k)
    {
      for (int j = 0; j < 2; ++j)
      {
        for (int i=0; i<2; ++i)
        {
          int nodeIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
          transformData.data[I][0] = mesh.nodeCoord(nodeIdx, 0);
          transformData.data[I][1] = mesh.nodeCoord(nodeIdx, 1);
          transformData.data[I][2] = mesh.nodeCoord(nodeIdx, 2);
          ++I;
        }
      }
    }
  }
}



} // namespace model_discretization_interface
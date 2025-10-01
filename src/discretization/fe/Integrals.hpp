#pragma once

#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
// #include "finiteElement/shiva/SEMQkGLIntegralsShiva.hpp"

template <int ORDER, int METHOD_TYPE>
struct IntegralTypeSelector;

namespace IntegralType
{
enum
{
  MAKUTU,
  SHIVA
};
}

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::MAKUTU>
{
  using type =
      typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
};

// template <int ORDER>
// struct IntegralTypeSelector<ORDER, IntegralType::SHIVA>
// {
//   using TransformType = LinearTransform<
//       float, InterpolatedShape<float, Cube<float>,
//                                LagrangeBasis<float, 1, EqualSpacing>,
//                                LagrangeBasis<float, 1, EqualSpacing>,
//                                LagrangeBasis<float, 1, EqualSpacing> > >;

//   using ParentElementType =
//       ParentElement<float, Cube<float>,
//                     LagrangeBasis<float, ORDER, EqualSpacing>,
//                     LagrangeBasis<float, ORDER, EqualSpacing>,
//                     LagrangeBasis<float, ORDER, EqualSpacing> >;

//   using type = SEMQkGLIntegralsShiva<ORDER, TransformType,
//   ParentElementType>;
// };

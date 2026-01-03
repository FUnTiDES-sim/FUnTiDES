#pragma once

#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include "finiteElement/tensorial/Qk_Hexahedron_Tensorial.hpp"

#ifdef ENABLE_Shiva
#include "finiteElement/shiva/SEMQkGLIntegralsShiva.hpp"
#endif

template <int ORDER, int METHOD_TYPE>
struct IntegralTypeSelector;

namespace IntegralType
{
enum
{
  MAKUTU,
  SHIVA,
  TENSORIAL
};
}

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::MAKUTU>
{
  using type =
      typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
};

#ifdef ENABLE_Shiva
template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::SHIVA>
{
  using TransformType = LinearTransform<
      float, InterpolatedShape<float, Cube<float>,
                               LagrangeBasis<float, 1, EqualSpacing>,
                               LagrangeBasis<float, 1, EqualSpacing>,
                               LagrangeBasis<float, 1, EqualSpacing> > >;

  using ParentElementType =
      ParentElement<float, Cube<float>,
                    LagrangeBasis<float, ORDER, EqualSpacing>,
                    LagrangeBasis<float, ORDER, EqualSpacing>,
                    LagrangeBasis<float, ORDER, EqualSpacing> >;

  using type = SEMQkGLIntegralsShiva<ORDER, TransformType, ParentElementType>;
};
#endif

// Tensorial selector - uses O(n^4) algorithm
template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::TENSORIAL>
{
  using type =
      typename Qk_Hexahedron_Lagrange_GaussLobatto_Tensorial_Selector<ORDER>::type;
};

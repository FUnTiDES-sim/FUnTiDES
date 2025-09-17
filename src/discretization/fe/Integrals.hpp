#pragma once

#include "finiteElement/classic/SEMQkGLIntegralsClassic.hpp"
#include "finiteElement/geos/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include "finiteElement/optim/SEMQkGLIntegralsOptim.hpp"
#include "finiteElement/shiva/SEMQkGLIntegralsShiva.hpp"

template <int ORDER, int METHOD_TYPE>
struct IntegralTypeSelector;

namespace IntegralType
{
enum
{
  CLASSIC,
  OPTIM,
  GEOS,
  SHIVA
};
}

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::CLASSIC>
{
  using type = SEMQkGLIntegralsClassic<ORDER>;
};

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::OPTIM>
{
  using type = SEMQkGLIntegralsOptim<ORDER, float, float>;
};

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::GEOS>
{
  using type =
      typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
};

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

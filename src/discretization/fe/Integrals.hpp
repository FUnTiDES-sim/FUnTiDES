#pragma once

#include "finiteElement/classic/SEMQkGLIntegralsClassic.hpp"
#include "finiteElement/optim/SEMQkGLIntegralsOptim.hpp"
#include "finiteElement/geos/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include "finiteElement/shiva/SEMQkGLIntegralsShiva.hpp"

template < int ORDER, int METHOD_TYPE >
struct IntegralTypeSelector
{
  static_assert( true, "Unsupported method type" );  
};


template< int ORDER >
struct IntegralTypeSelector< ORDER, 0 >
{
  using type = SEMQkGLIntegralsClassic<ORDER>; // Default integral type for methodType 0
};

template< int ORDER >
struct IntegralTypeSelector< ORDER, 1 >
{
  using type = SEMQkGLIntegralsOptim< ORDER, float, float >; // Example for methodType 1
};

template< int ORDER >
struct IntegralTypeSelector< ORDER, 2 >
{
  using type = Q2_Hexahedron_Lagrange_GaussLobatto; // Example for methodType 1
};


template< int ORDER >
struct IntegralTypeSelector< ORDER, 3 >
{
    using TransformType =
    LinearTransform< float,
                     InterpolatedShape< float,
                                        Cube< float >,
                                        LagrangeBasis< float, 1, EqualSpacing >,
                                        LagrangeBasis< float, 1, EqualSpacing >,
                                        LagrangeBasis< float, 1, EqualSpacing > > >;

  using ParentElementType =
    ParentElement< float,
                   Cube< float >,
                   LagrangeBasis< float, ORDER, EqualSpacing >,
                   LagrangeBasis< float, ORDER, EqualSpacing >,
                   LagrangeBasis< float, ORDER, EqualSpacing > >;

  using type = SEMQkGLIntegralsShiva< ORDER, TransformType, ParentElementType >;
};

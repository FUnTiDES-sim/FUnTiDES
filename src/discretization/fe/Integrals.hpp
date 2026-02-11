#pragma once

#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

template <int ORDER, int METHOD_TYPE>
struct IntegralTypeSelector;

namespace IntegralType
{
enum
{
  MAKUTU
};
}

template <int ORDER>
struct IntegralTypeSelector<ORDER, IntegralType::MAKUTU>
{
  using type =
      typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
};

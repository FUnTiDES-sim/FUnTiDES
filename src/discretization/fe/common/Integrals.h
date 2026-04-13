#ifndef FUNTIDES_DISCRETIZATION_FE_INTEGRALS_H_
#define FUNTIDES_DISCRETIZATION_FE_INTEGRALS_H_

#include "../makutu/include/Qk_Hexahedron_Lagrange_GaussLobatto.h"

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

#endif  // FUNTIDES_DISCRETIZATION_FE_INTEGRALS_H_

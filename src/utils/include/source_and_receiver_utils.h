#ifndef SOURCEANDRECEIVERUTILS_HPP
#define SOURCEANDRECEIVERUTILS_HPP

#include "../../src/discretization/fe/SEMKernels/src/finiteElement/geos/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include "data_type.h"

using namespace std::chrono;

namespace SourceAndReceiverUtils
{

template <int ORDER>
void ComputeRHSWeights(real_t const (&cornerCoords)[8][3], float coordsReal[3],
                       ARRAY_REAL_VIEW& rhsWeights)
{
  constexpr int numNodes =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;

  // Compute coordinates on reference element
  double coordsRef[3]{};
  float invJ[3][3] = {{0}};
  Qk_Hexahedron_Lagrange_GaussLobatto_Selector<
      ORDER>::type::invJacobianTransformation(0, cornerCoords, invJ);
  for (int i = 0; i < 3; i++)
  {
    coordsRef[i] = -1.0;
    for (int j = 0; j < 3; j++)
    {
      coordsRef[i] += invJ[i][j] * (coordsReal[j] - cornerCoords[0][j]);
    }
  }

  // ComputeRhsWeights
  double N[numNodes] = {0};
  Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::calcN(coordsRef,
                                                                   N);
  for (int i = 0; i < numNodes; i++)
  {
    rhsWeights(0, i) = N[i];
  }
}

}  // namespace SourceAndReceiverUtils

#endif  // SOURCEANDRECEIVERUTILS_HPP

#ifndef SOURCEANDRECEIVERUTILS_HPP
#define SOURCEANDRECEIVERUTILS_HPP

#include <array>

#include "data_type.h"
#include "finiteElement/makutu/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

using namespace std::chrono;

namespace SourceAndReceiverUtils
{

template <int ORDER>
void ComputeRHSWeights(real_t const (&cornerCoords)[8][3],
                       std::array<float, 3> coordsReal,
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

template <int ORDER>
void ComputeSourceMomentWeights(real_t const (&cornerCoords)[8][3],
                                std::array<float, 3> coordsReal,
                                ARRAY_REAL_VIEW& rhsWeightsX,
                                ARRAY_REAL_VIEW& rhsWeightsY,
                                ARRAY_REAL_VIEW& rhsWeightsZ,
                                double M0)
{
  constexpr int numNodes =
      Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;

  real_t coordsRef[3]{};
  float invJ[3][3] = {{0}};
  Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::
      invJacobianTransformation(0, cornerCoords, invJ);
  for (int i = 0; i < 3; i++)
  {
    coordsRef[i] = -1.0f;
    for (int j = 0; j < 3; j++)
    {
      coordsRef[i] += invJ[i][j] * (coordsReal[j] - cornerCoords[0][j]);
    }
  }

  real_t gradN[numNodes][3] = {0};
  Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::
      calcGradNWithCorners(coordsRef, cornerCoords, gradN);

  for (int i = 0; i < numNodes; i++)
  {
    rhsWeightsX(0, i) = M0 * gradN[i][0];
    rhsWeightsY(0, i) = M0 * gradN[i][1];
    rhsWeightsZ(0, i) = M0 * gradN[i][2];
  }
}

}  // namespace SourceAndReceiverUtils

#endif  // SOURCEANDRECEIVERUTILS_HPP

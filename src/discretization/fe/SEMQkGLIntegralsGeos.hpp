#ifndef SEMQKINTEGRALSGEOS_HPP_
#define SEMQKINTEGRALSGEOS_HPP_

#include <SEMdata.hpp>
#include <SEMmacros.hpp>
#include <dataType.hpp>
#include "Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

// class SEMQkGLIntegralsGeos
// {
// private:
//   int order;
//   struct SEMinfo infos;

//   constexpr static double sqrt5 = 2.2360679774997897;
//   // order of polynomial approximation
//   constexpr static int ORDER = SEMinfo::myOrderNumber;
//   // Half the number of support points, rounded down. Precomputed for efficiency
//   constexpr static int halfNodes = ORDER / 2;

// public:
//   PROXY_HOST_DEVICE SEMQkGLIntegralsGeos() {};
//   PROXY_HOST_DEVICE ~SEMQkGLIntegralsGeos(){};

//   template <int ORDER>
//   PROXY_HOST_DEVICE void computeMassMatrixAndStiffnessVector(
//       const int &elementNumber, const int &nPointsPerElement,
//       ARRAY_REAL_VIEW const &nodesCoordsX, ARRAY_REAL_VIEW const &nodesCoordsY,
//       ARRAY_REAL_VIEW const &nodesCoordsZ, float massMatrixLocal[],
//       float pnLocal[], float Y[]) const {
//     float X[8][3];
//     int I = 0;

//     for (int k = 0; k < 2; k++) {
//       for (int j = 0; j < 2; j++) {
//         for (int i = 0; i < 2; i++) {
//           int l = i + j * 2 + k * 4;
//           X[I][0] = nodesCoordsX(elementNumber, l);
//           X[I][1] = nodesCoordsZ(elementNumber, l);
//           X[I][2] = nodesCoordsY(elementNumber, l);
//           I++;
//         }
//       }
//     }

//     for (int q = 0; q < nPointsPerElement; q++) {
//       Y[q] = 0;
//     }

//     for (int q = 0; q < nPointsPerElement; q++) {
//       massMatrixLocal[q] = computeMassTerm<ORDER>(q, X);
//       computeStiffnessTerm(
//           q, X, [&](const int i, const int j, const double val) {
//             float localIncrement = val * pnLocal[j];
//             Y[i] += localIncrement;
//           });
//     }
//   }
// };





#endif // SEMQKINTEGRALSGEOS_HPP_

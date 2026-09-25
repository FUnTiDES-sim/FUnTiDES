/**
 * @file source_and_receiver_utils.h
 * @brief Utilities to compute interpolation weights for point sources/receivers and
 * Distributed Acoustic Sensing (DAS) sample weights within a single spectral element.
 */
#ifndef FUNTIDES_ACQUISITION_INCLUDE_SOURCE_AND_RECEIVER_UTILS_H_
#define FUNTIDES_ACQUISITION_INCLUDE_SOURCE_AND_RECEIVER_UTILS_H_
#include <array>
#include <cmath>

#include "Qk_Hexahedron_Lagrange_GaussLobatto.h"
#include "data_type.h"

using namespace std::chrono;

namespace SourceAndReceiverUtils {

/**
 * @brief Maps a polynomial order to the corresponding 1D Lagrange basis type.
 * @tparam ORDER polynomial order.
 */
template <int ORDER>
struct LagrangeBasisSelector;

template <>
struct LagrangeBasisSelector<1> {
  using type = LagrangeBasis1;
};
template <>
struct LagrangeBasisSelector<2> {
  using type = LagrangeBasis2;
};
template <>
struct LagrangeBasisSelector<3> {
  using type = LagrangeBasis3GL;
};
template <>
struct LagrangeBasisSelector<4> {
  using type = LagrangeBasis4GL;
};
template <>
struct LagrangeBasisSelector<5> {
  using type = LagrangeBasis5GL;
};
template <>
struct LagrangeBasisSelector<6> {
  using type = LagrangeBasis6GL;
};

/**
 * @brief Receiver type for Distributed Acoustic Sensing (DAS) computations.
 */
enum class DASType {
  kNone,               ///< Standard point receiver (no DAS).
  kDipole,             ///< Displacement difference between two endpoints divided by gauge length.
  kStrainIntegration,  ///< Strain integrated along the fiber direction.
};

/**
 * @brief Compute the fiber direction unit vector from dip and azimuth angles.
 * @param[in] dip_deg dip angle in degrees, measured from horizontal.
 * @param[in] azimuth_deg azimuth angle in degrees, measured in the XY plane from the X axis.
 * @return unit vector {x, y, z} giving the fiber direction.
 */
inline std::array<float, 3> ComputeDASVector(float dip_deg, float azimuth_deg) {
  const float dip = dip_deg * static_cast<float>(M_PI) / 180.0f;
  const float az = azimuth_deg * static_cast<float>(M_PI) / 180.0f;
  return {std::cos(dip) * std::cos(az), std::cos(dip) * std::sin(az), std::sin(dip)};
}

/**
 * @brief Compute the interpolation (shape function) weights for a point located inside a
 * single hexahedral element, for use as right hand side weights.
 * @tparam ORDER polynomial order of the element.
 * @param[in] cornerCoords 8 corner coordinates of the element.
 * @param[in] coordsReal physical coordinates of the point.
 * @param[out] rhsWeights receives the basis function values at the point, indexed by node;
 * existing content is overwritten, not accumulated.
 * @todo VERIFY: what does row index 0 of rhsWeights represent (e.g. a fixed field component)?
 */
template <int ORDER>
void ComputeRHSWeights(real_t const (&cornerCoords)[8][3], std::array<float, 3> coordsReal,
                       arrayReal::host_mirror_type& rhsWeights) {
  constexpr int numNodes = Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::numNodes;

  // Map physical coordinates to the reference element via the inverse Jacobian at node 0.
  double coordsRef[3]{};
  float invJ[3][3] = {{0}};
  Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::invJacobianTransformation(0, cornerCoords, invJ);
  for (int i = 0; i < 3; i++) {
    coordsRef[i] = -1.0;
    for (int j = 0; j < 3; j++) {
      coordsRef[i] += invJ[i][j] * (coordsReal[j] - cornerCoords[0][j]);
    }
  }

  double N[numNodes] = {0};
  Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type::calcN(coordsRef, N);
  for (int i = 0; i < numNodes; i++) {
    rhsWeights(0, i) = N[i];
  }
}

/**
 * @brief Compute the DAS weights contributed by a single quadrature sample point within an
 * element, and accumulate them into dasWeights.
 * In dipole mode, weights are the basis function values scaled by integrationConstant.
 * In strain integration mode, weights are the basis function gradients projected onto
 * direction, scaled by integrationConstant.
 * @tparam ORDER polynomial order of the element.
 * @param[in] cornerCoords 8 corner coordinates of the element.
 * @param[in] coordsReal physical coordinates of the sample point.
 * @param[in] direction fiber direction unit vector.
 * @param[in] integrationConstant quadrature weight for this sample point.
 * @param[in] dasType dipole or strain integration mode.
 * @param[in,out] dasWeights array of size numNodes; contributions are added to the existing
 * values, so the caller must initialize it before the first call.
 */
template <int ORDER>
void ComputeDASWeightsForSample(real_t const (&cornerCoords)[8][3], std::array<float, 3> coordsReal,
                                const std::array<float, 3>& direction, float integrationConstant, DASType dasType,
                                float* dasWeights) {
  using FEType = typename Qk_Hexahedron_Lagrange_GaussLobatto_Selector<ORDER>::type;
  using Basis = typename LagrangeBasisSelector<ORDER>::type;
  constexpr int numNodes = FEType::numNodes;

  // Map physical coordinates to reference element
  real_t invJ[3][3] = {{0}};
  FEType::invJacobianTransformation(0, cornerCoords, invJ);
  double coordsRef[3]{};
  for (int i = 0; i < 3; i++) {
    coordsRef[i] = -1.0;
    for (int j = 0; j < 3; j++) {
      coordsRef[i] += invJ[i][j] * (coordsReal[j] - cornerCoords[0][j]);
    }
  }

  if (dasType == DASType::kStrainIntegration) {
    // Compute parent gradients (dN/dxi) at coordsRef for each node and transform them to
    // physical gradients using the inverse Jacobian. Basis::value/gradient return double;
    // results are cast explicitly to avoid narrowing warnings in the FE kernel.
    for (int c = 0; c < ORDER + 1; ++c) {
      for (int b = 0; b < ORDER + 1; ++b) {
        for (int a_idx = 0; a_idx < ORDER + 1; ++a_idx) {
          int nodeIndex = a_idx + b * (ORDER + 1) + c * (ORDER + 1) * (ORDER + 1);

          double Na = Basis::value(a_idx, coordsRef[0]);
          double Nb = Basis::value(b, coordsRef[1]);
          double Nc = Basis::value(c, coordsRef[2]);
          double dNa = Basis::gradient(a_idx, coordsRef[0]);
          double dNb = Basis::gradient(b, coordsRef[1]);
          double dNc = Basis::gradient(c, coordsRef[2]);

          real_t dNdXi[3] = {static_cast<real_t>(dNa * Nb * Nc), static_cast<real_t>(Na * dNb * Nc),
                             static_cast<real_t>(Na * Nb * dNc)};

          // Transform to physical gradients: gradN_i = sum_j invJ[j][i] *
          // dNdXi[j]
          real_t gx = 0, gy = 0, gz = 0;
          for (int j = 0; j < 3; ++j) {
            gx += invJ[j][0] * dNdXi[j];
            gy += invJ[j][1] * dNdXi[j];
            gz += invJ[j][2] * dNdXi[j];
          }
          dasWeights[nodeIndex] += (gx * direction[0] + gy * direction[1] + gz * direction[2]) * integrationConstant;
        }
      }
    }
  } else {
    // Dipole mode: weight is the basis function value at the sample point.
    double N[numNodes] = {0};
    FEType::calcN(coordsRef, N);
    for (int a = 0; a < numNodes; ++a) {
      dasWeights[a] += static_cast<float>(N[a]) * integrationConstant;
    }
  }
}

}  // namespace SourceAndReceiverUtils
#endif  // FUNTIDES_ACQUISITION_INCLUDE_SOURCE_AND_RECEIVER_UTILS_H_

#ifndef FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_PENALTY_H_
#define FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_PENALTY_H_

#include <cmath>

#include "common_macros.h"
#include "data_type.h"

namespace solver {
namespace fe {

/**
 * @brief Compute the face area of a bilinear quadrilateral from its 4 corners.
 *
 * Uses the cross product of the two diagonals: area = 0.5 * |d1 x d2|.
 *
 * @param X  Corner coordinates [4][3].
 * @return   Face area.
 */
// PROXY_HOST_DEVICE real_t computeFaceArea(real_t const (&X)[4][3]) {
//   real_t d1[3], d2[3];
//   for (int k = 0; k < 3; ++k) {
//     d1[k] = X[2][k] - X[0][k];
//     d2[k] = X[3][k] - X[1][k];
//   }
//   real_t const cx = d1[1] * d2[2] - d1[2] * d2[1];
//   real_t const cy = d1[2] * d2[0] - d1[0] * d2[2];
//   real_t const cz = d1[0] * d2[1] - d1[1] * d2[0];
//   return 0.5f * sqrt(cx * cx + cy * cy + cz * cz);
// }

PROXY_HOST_DEVICE real_t computeFaceArea(const real_t (&X)[4][3]) {
  real_t area = 0.0;

  for (int i = 0; i < 4; ++i) {
    int j = (i + 1) % 4;

    real_t edge1[3], edge2[3];

    for (int k = 0; k < 3; ++k) {
      edge1[k] = X[i][k] - X[0][k];
      edge2[k] = X[j][k] - X[0][k];
    }

    real_t cx = edge1[1] * edge2[2] - edge1[2] * edge2[1];
    real_t cy = edge1[2] * edge2[0] - edge1[0] * edge2[2];
    real_t cz = edge1[0] * edge2[1] - edge1[1] * edge2[0];

    area += 0.5f * sqrt(cx * cx + cy * cy + cz * cz);
  }

  return area;
}

/**
 * @brief Estimate the volume of a trilinear hexahedron from its 8 corners.
 *
 * Approximates vol ≈ 8 * |det(J_center)| where J is evaluated at (0,0,0)
 * in the reference element [-1,1]^3.
 *
 * Corner ordering: X[iv + 2*jv + 4*kv] for iv,jv,kv in {0,1}.
 *
 * @param X  Corner coordinates [8][3].
 * @return   Approximate element volume.
 */
PROXY_HOST_DEVICE real_t computeHexVolume(real_t const (&X)[8][3]) {
  real_t dxi[3], deta[3], dzeta[3];
  for (int k = 0; k < 3; ++k) {
    dxi[k] = 0.25f * (-X[0][k] + X[1][k] - X[2][k] + X[3][k] - X[4][k] + X[5][k] - X[6][k] + X[7][k]);
    deta[k] = 0.25f * (-X[0][k] - X[1][k] + X[2][k] + X[3][k] - X[4][k] - X[5][k] + X[6][k] + X[7][k]);
    dzeta[k] = 0.25f * (-X[0][k] - X[1][k] - X[2][k] - X[3][k] + X[4][k] + X[5][k] + X[6][k] + X[7][k]);
  }
  real_t const det = dxi[0] * (deta[1] * dzeta[2] - deta[2] * dzeta[1]) -
                     dxi[1] * (deta[0] * dzeta[2] - deta[2] * dzeta[0]) +
                     dxi[2] * (deta[0] * dzeta[1] - deta[1] * dzeta[0]);
  return 8.0f * fabs(det);
}

/**
 * @brief Compute the SIPG penalty parameter gamma for a face.
 *
 * Formula: gamma = penalty_factor * p * (p+1) / h_f
 * where h_f = vol_e / area_f is the characteristic face length scale,
 * and penalty_factor is a user-supplied constant (>= 1, typically 10).
 *
 * @tparam ORDER      Polynomial order p.
 * @param faceCoords  Face corner coordinates [4][3].
 * @param X8          Element corner coordinates [8][3].
 * @param penalty_factor  Dimensionless constant (must be large enough for
 *                        coercivity; problem-dependent).
 * @return            Penalty parameter gamma.
 */
template <int ORDER>
PROXY_HOST_DEVICE real_t computeSIPGPenalty(real_t const (&faceCoords)[4][3], real_t const (&X8)[8][3],
                                            real_t const penalty_factor) {
  real_t const area = computeFaceArea(faceCoords);
  real_t const vol = computeHexVolume(X8);
  real_t const h_f = vol / area;
  return penalty_factor * static_cast<real_t>((ORDER + 1) * (ORDER + 1)) / h_f;
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_PENALTY_H_

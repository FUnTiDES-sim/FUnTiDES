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
 * Splits the quad on the X[0]-X[2] diagonal and sums the two triangle areas,
 * 0.5 * |(Xa - X0) x (Xb - X0)|.
 *
 * @param X  Corner coordinates [4][3].
 * @return   Face area.
 */
PROXY_HOST_DEVICE real_t computeFaceArea(const real_t (&X)[4][3]) {
  real_t d1[3], d2[3], d3[3];
  for (int k = 0; k < 3; ++k) {
    d1[k] = X[1][k] - X[0][k];
    d2[k] = X[2][k] - X[0][k];
    d3[k] = X[3][k] - X[0][k];
  }

  real_t const ax = d1[1] * d2[2] - d1[2] * d2[1];
  real_t const ay = d1[2] * d2[0] - d1[0] * d2[2];
  real_t const az = d1[0] * d2[1] - d1[1] * d2[0];

  real_t const bx = d2[1] * d3[2] - d2[2] * d3[1];
  real_t const by = d2[2] * d3[0] - d2[0] * d3[2];
  real_t const bz = d2[0] * d3[1] - d2[1] * d3[0];

  return 0.5f * (sqrt(ax * ax + ay * ay + az * az) + sqrt(bx * bx + by * by + bz * bz));
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
 * @brief Compute the SIPG penalty parameter gamma for a face of known area.
 *
 * Formula: gamma = penalty_factor * (p+1)^2 / h_f
 * where h_f = vol_e / area_f is the characteristic face length scale,
 * and penalty_factor is a user-supplied constant (>= 1, typically 10).
 *
 * Takes the face area rather than the corners so that callers computing gamma
 * on both sides of the same face (the two elements share area_f, only X8
 * differs) pay for it once.
 *
 * @tparam ORDER      Polynomial order p.
 * @param area        Face area, from computeFaceArea().
 * @param X8          Element corner coordinates [8][3].
 * @param penalty_factor  Dimensionless constant (must be large enough for
 *                        coercivity; problem-dependent).
 * @return            Penalty parameter gamma.
 */
template <int ORDER>
PROXY_HOST_DEVICE real_t computeSIPGPenaltyFromArea(real_t const area, real_t const (&X8)[8][3],
                                                    real_t const penalty_factor) {
  real_t const h_f = computeHexVolume(X8) / area;
  return penalty_factor * static_cast<real_t>((ORDER + 1) * (ORDER + 1)) / h_f;
}

/**
 * @brief Compute the SIPG penalty parameter gamma for a face.
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
  return computeSIPGPenaltyFromArea<ORDER>(computeFaceArea(faceCoords), X8, penalty_factor);
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SRC_SOLVER_FE_DG_IMPL_COMMON_INCLUDE_DG_PENALTY_H_

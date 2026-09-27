#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_ELASTIC_FLUX_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_ELASTIC_FLUX_H_

#include "common_macros.h"
#include "data_type_kokkos.h"

namespace solver {
namespace fe {
namespace flux {

/**
 * @brief Physical displacement gradient H[i][t] = du_t/dx_i.
 * @param[in] J_inv Inverse Jacobian, J_inv[r][i] = dxi_r/dx_i.
 * @param[in] grad_u_ref Reference gradient, grad_u_ref[r][t] = du_t/dxi_r.
 * @param[out] H Physical gradient.
 */
PROXY_HOST_DEVICE void physicalGradient(float const (&J_inv)[3][3], float const (&grad_u_ref)[3][3], float (&H)[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int t = 0; t < 3; ++t)
      H[i][t] = J_inv[0][i] * grad_u_ref[0][t] + J_inv[1][i] * grad_u_ref[1][t] + J_inv[2][i] * grad_u_ref[2][t];
}

/**
 * @brief Pull the stress back to the reference element: flux[p][s] = sum_i J_inv[p][i] * sigma[i][s].
 * @param[in] J_inv Inverse Jacobian, J_inv[p][i] = dxi_p/dx_i.
 * @param[in] sigma Cauchy stress, symmetric.
 * @param[out] flux Reference-element flux.
 */
PROXY_HOST_DEVICE void pullBackStress(float const (&J_inv)[3][3], float const (&sigma)[3][3], float (&flux)[3][3]) {
  for (int p = 0; p < 3; ++p)
    for (int s = 0; s < 3; ++s)
      flux[p][s] = J_inv[p][0] * sigma[0][s] + J_inv[p][1] * sigma[1][s] + J_inv[p][2] * sigma[2][s];
}

/**
 * @brief Isotropic elastic flux: reference gradient to reference-element flux.
 *
 * The Cauchy stress sigma = mu * (H + H^T) + lambda * tr(H) * I is built as a
 * symmetric tensor, then pulled back, so the map grad_u_ref -> flux carries the
 * major symmetry of the stiffness operator by construction.
 *
 * @param[in] J_inv Inverse Jacobian, J_inv[r][i] = dxi_r/dx_i.
 * @param[in] mu Lame parameter mu.
 * @param[in] lambda Lame parameter lambda.
 * @param[in] grad_u_ref Reference displacement gradient, grad_u_ref[r][t] = du_t/dxi_r.
 * @param[out] flux Reference-element flux.
 */
PROXY_HOST_DEVICE void elasticFluxIso(float const (&J_inv)[3][3], float mu, float lambda,
                                      float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
  float H[3][3];
  physicalGradient(J_inv, grad_u_ref, H);
  float const div = H[0][0] + H[1][1] + H[2][2];
  float sigma[3][3];
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) sigma[a][b] = mu * (H[a][b] + H[b][a]);
  sigma[0][0] += lambda * div;
  sigma[1][1] += lambda * div;
  sigma[2][2] += lambda * div;
  pullBackStress(J_inv, sigma, flux);
}

/**
 * @brief VTI elastic flux: reference gradient to reference-element flux.
 *
 * The stress is built as a symmetric tensor from the six independent VTI
 * stiffness coefficients (Voigt notation), then pulled back.
 *
 * @param[in] J_inv Inverse Jacobian, J_inv[r][i] = dxi_r/dx_i.
 * @param[in] c11 Stiffness coefficient C11.
 * @param[in] c12 Stiffness coefficient C12.
 * @param[in] c13 Stiffness coefficient C13.
 * @param[in] c33 Stiffness coefficient C33.
 * @param[in] c44 Stiffness coefficient C44.
 * @param[in] c66 Stiffness coefficient C66.
 * @param[in] grad_u_ref Reference displacement gradient, grad_u_ref[r][t] = du_t/dxi_r.
 * @param[out] flux Reference-element flux.
 * @todo VERIFY: which axis is the VTI symmetry axis (the formulas use z)?
 */
PROXY_HOST_DEVICE void elasticFluxVti(float const (&J_inv)[3][3], float c11, float c12, float c13, float c33, float c44,
                                      float c66, float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
  float H[3][3];
  physicalGradient(J_inv, grad_u_ref, H);
  float sigma[3][3];
  sigma[0][0] = c11 * H[0][0] + c12 * H[1][1] + c13 * H[2][2];
  sigma[1][1] = c12 * H[0][0] + c11 * H[1][1] + c13 * H[2][2];
  sigma[2][2] = c13 * (H[0][0] + H[1][1]) + c33 * H[2][2];
  sigma[1][2] = sigma[2][1] = c44 * (H[1][2] + H[2][1]);
  sigma[0][2] = sigma[2][0] = c44 * (H[0][2] + H[2][0]);
  sigma[0][1] = sigma[1][0] = c66 * (H[0][1] + H[1][0]);
  pullBackStress(J_inv, sigma, flux);
}

/**
 * @brief TTI elastic flux: reference gradient to reference-element flux.
 *
 * Uses the full 6x6 stiffness matrix in Voigt notation with engineering shear
 * strains (order xx, yy, zz, yz, xz, xy). The stress tensor is assembled
 * symmetric, so the map grad_u_ref -> flux carries the major symmetry of the
 * stiffness operator provided CTTI is symmetric.
 *
 * @param[in] J_inv Inverse Jacobian, J_inv[r][i] = dxi_r/dx_i.
 * @param[in] CTTI Symmetric 6x6 stiffness matrix in Voigt notation.
 * @param[in] grad_u_ref Reference displacement gradient, grad_u_ref[r][t] = du_t/dxi_r.
 * @param[out] flux Reference-element flux.
 */
PROXY_HOST_DEVICE void elasticFluxTti(float const (&J_inv)[3][3], float const (&CTTI)[6][6],
                                      float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
  float H[3][3];
  physicalGradient(J_inv, grad_u_ref, H);
  float const eps[6] = {H[0][0], H[1][1], H[2][2], H[1][2] + H[2][1], H[0][2] + H[2][0], H[0][1] + H[1][0]};
  float sv[6];
  for (int a = 0; a < 6; ++a) {
    sv[a] = 0.0f;
    for (int b = 0; b < 6; ++b) sv[a] += CTTI[a][b] * eps[b];
  }
  float const sigma[3][3] = {{sv[0], sv[5], sv[4]}, {sv[5], sv[1], sv[3]}, {sv[4], sv[3], sv[2]}};
  pullBackStress(J_inv, sigma, flux);
}

/// Number of floats of the compact TTI description used by elasticFluxTtiCompact().
constexpr int kTtiCompactSize = 8;

/**
 * @brief Compact TTI description from the VTI coefficients and the symmetry axis.
 *
 * A TTI medium is a VTI medium whose symmetry axis is the unit vector n, so its stress is
 *   sigma = (lambda tr(e) + alpha e_n) I + 2 mu_t e + (alpha tr(e) + beta e_n) n n^T
 *         + 2 (mu_l - mu_t) (n w^T + w n^T),
 * with e the strain, w = e n and e_n = n.w. Five coefficients and n replace the 21 entries
 * of the rotated 6x6 tensor.
 *
 * @param[in] c11 VTI stiffness coefficient C11.
 * @param[in] c13 VTI stiffness coefficient C13.
 * @param[in] c33 VTI stiffness coefficient C33.
 * @param[in] c44 VTI stiffness coefficient C44.
 * @param[in] c66 VTI stiffness coefficient C66 (C12 = C11 - 2 C66).
 * @param[in] n Unit symmetry axis.
 * @param[out] p lambda, mu_t, alpha, mu_l - mu_t, beta, n.
 */
PROXY_HOST_DEVICE void ttiCompactFromVti(float c11, float c13, float c33, float c44, float c66, float const (&n)[3],
                                         float (&p)[kTtiCompactSize]) {
  float const lambda = c11 - 2.0f * c66;
  p[0] = lambda;
  p[1] = c66;
  p[2] = c13 - lambda;
  p[3] = c44 - c66;
  p[4] = c11 + c33 - 2.0f * c13 - 4.0f * c44;
  p[5] = n[0];
  p[6] = n[1];
  p[7] = n[2];
}

/**
 * @brief TTI elastic flux from the compact description built by ttiCompactFromVti().
 *
 * Same result as elasticFluxTti() with the equivalent rotated 6x6 tensor.
 *
 * @param[in] J_inv Inverse Jacobian, J_inv[r][i] = dxi_r/dx_i.
 * @param[in] p Compact TTI description.
 * @param[in] grad_u_ref Reference displacement gradient, grad_u_ref[r][t] = du_t/dxi_r.
 * @param[out] flux Reference-element flux.
 */
PROXY_HOST_DEVICE void elasticFluxTtiCompact(float const (&J_inv)[3][3], float const (&p)[kTtiCompactSize],
                                             float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
  float H[3][3];
  physicalGradient(J_inv, grad_u_ref, H);
  float const e[3][3] = {{H[0][0], 0.5f * (H[0][1] + H[1][0]), 0.5f * (H[0][2] + H[2][0])},
                         {0.5f * (H[0][1] + H[1][0]), H[1][1], 0.5f * (H[1][2] + H[2][1])},
                         {0.5f * (H[0][2] + H[2][0]), 0.5f * (H[1][2] + H[2][1]), H[2][2]}};
  float const n[3] = {p[5], p[6], p[7]};
  float w[3];
  for (int i = 0; i < 3; ++i) w[i] = e[i][0] * n[0] + e[i][1] * n[1] + e[i][2] * n[2];
  float const tr = e[0][0] + e[1][1] + e[2][2];
  float const en = n[0] * w[0] + n[1] * w[1] + n[2] * w[2];

  float const diag = p[0] * tr + p[2] * en;
  float const two_mu_t = 2.0f * p[1];
  float const axial = p[2] * tr + p[4] * en;
  float const two_dmu = 2.0f * p[3];
  float sigma[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j) {
      float const s = two_mu_t * e[i][j] + axial * n[i] * n[j] + two_dmu * (n[i] * w[j] + w[i] * n[j]);
      sigma[i][j] = sigma[j][i] = (i == j) ? s + diag : s;
    }
  pullBackStress(J_inv, sigma, flux);
}

}  // namespace flux
}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_ELASTIC_FLUX_H_

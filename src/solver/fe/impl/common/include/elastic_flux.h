#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_ELASTIC_FLUX_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_ELASTIC_FLUX_H_

#include "common_macros.h"
#include "data_type_kokkos.h"

namespace solver {
namespace fe {
namespace flux {

/**
 * @brief Physical displacement gradient H[i][t] = du_t/dx_i.
 * @param J_inv Inverse Jacobian, J_inv[r][i] = dxi_r/dx_i.
 * @param grad_u_ref Reference gradient, grad_u_ref[r][t] = du_t/dxi_r.
 * @param H Output physical gradient.
 */
PROXY_HOST_DEVICE void physicalGradient(float const (&J_inv)[3][3], float const (&grad_u_ref)[3][3], float (&H)[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int t = 0; t < 3; ++t)
      H[i][t] = J_inv[0][i] * grad_u_ref[0][t] + J_inv[1][i] * grad_u_ref[1][t] + J_inv[2][i] * grad_u_ref[2][t];
}

/**
 * @brief Pull the stress back to the reference element: flux[p][s] = sum_i J_inv[p][i] * sigma[i][s].
 * @param J_inv Inverse Jacobian, J_inv[p][i] = dxi_p/dx_i.
 * @param sigma Cauchy stress (symmetric).
 * @param flux Output reference-element flux.
 */
PROXY_HOST_DEVICE void pullBackStress(float const (&J_inv)[3][3], float const (&sigma)[3][3], float (&flux)[3][3]) {
  for (int p = 0; p < 3; ++p)
    for (int s = 0; s < 3; ++s)
      flux[p][s] = J_inv[p][0] * sigma[0][s] + J_inv[p][1] * sigma[1][s] + J_inv[p][2] * sigma[2][s];
}

/**
 * @brief Isotropic elastic flux: reference gradient to reference-element flux.
 *
 * Assembles the Cauchy stress once as a symmetric tensor
 * sigma = mu * (H + H^T) + lambda * tr(H) * I, then pulls it back. Because sigma
 * is built symmetric, the resulting linear map grad_u_ref -> flux carries the
 * major symmetry of the stiffness operator by construction.
 *
 * @param J_inv Inverse Jacobian.
 * @param mu Lame parameter mu.
 * @param lambda Lame parameter lambda.
 * @param grad_u_ref Reference displacement gradient, grad_u_ref[r][t] = du_t/dxi_r.
 * @param flux Output reference-element flux.
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
 * Stress from the six independent VTI stiffness coefficients (Voigt), built as a
 * symmetric tensor before the pull-back.
 *
 * @param J_inv Inverse Jacobian.
 * @param c11 c12 c13 c33 c44 c66 VTI stiffness coefficients.
 * @param grad_u_ref Reference displacement gradient.
 * @param flux Output reference-element flux.
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
 * Full 6x6 stiffness matrix in Voigt notation with engineering shear strains.
 * CTTI is expected symmetric; the assembled stress tensor is symmetric, so the
 * map grad_u_ref -> flux inherits the major symmetry of the stiffness operator.
 *
 * @param J_inv Inverse Jacobian.
 * @param CTTI Symmetric 6x6 stiffness matrix (Voigt).
 * @param grad_u_ref Reference displacement gradient.
 * @param flux Output reference-element flux.
 */
PROXY_HOST_DEVICE void elasticFluxTti(float const (&J_inv)[3][3], float const (&CTTI)[6][6],
                                      float const (&grad_u_ref)[3][3], float (&flux)[3][3]) {
  float H[3][3];
  physicalGradient(J_inv, grad_u_ref, H);
  // Strain in Voigt notation (engineering shear strains).
  float const eps[6] = {H[0][0], H[1][1], H[2][2], H[1][2] + H[2][1], H[0][2] + H[2][0], H[0][1] + H[1][0]};
  float sv[6];
  for (int a = 0; a < 6; ++a) {
    sv[a] = 0.0f;
    for (int b = 0; b < 6; ++b) sv[a] += CTTI[a][b] * eps[b];
  }
  float const sigma[3][3] = {{sv[0], sv[5], sv[4]}, {sv[5], sv[1], sv[3]}, {sv[4], sv[3], sv[2]}};
  pullBackStress(J_inv, sigma, flux);
}

}  // namespace flux
}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_ELASTIC_FLUX_H_

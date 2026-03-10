#ifndef FUNTIDES_MODEL_MESH_IMPL_COMMON_INCLUDE_ELASTICITY_UTILS_H_
#define FUNTIDES_MODEL_MESH_IMPL_COMMON_INCLUDE_ELASTICITY_UTILS_H_
#include <data_type.h>

/**
 * @brief Compute isotropic elasticity coefficients (Lamé parameters)
 */
template <typename FloatType>
PROXY_HOST_DEVICE void computeIsotropicCoefficients(FloatType vp, FloatType vs,
                                                    FloatType rho,
                                                    FloatType& lambda,
                                                    FloatType& mu)
{
  mu = rho * vs * vs;
  lambda = rho * (vp * vp - FloatType(2.0) * vs * vs);
}

/**
 * @brief Build isotropic elasticity tensor in Voigt notation
 */
template <typename FloatType>
PROXY_HOST_DEVICE void buildIsotropicTensor(FloatType lambda, FloatType mu,
                                            FloatType C[6][6])
{
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) C[i][j] = FloatType(0.0);

  FloatType const lambda_plus_2mu = lambda + FloatType(2.0) * mu;
  C[0][0] = lambda_plus_2mu;
  C[1][1] = lambda_plus_2mu;
  C[2][2] = lambda_plus_2mu;
  C[3][3] = mu;
  C[4][4] = mu;
  C[5][5] = mu;

  C[0][1] = lambda;
  C[1][0] = lambda;
  C[0][2] = lambda;
  C[2][0] = lambda;
  C[1][2] = lambda;
  C[2][1] = lambda;
}

/**
 * @brief Compute VTI elasticity coefficients (no rotation)
 */
template <typename FloatType>
PROXY_HOST_DEVICE void computeVTICoefficients(
    FloatType vp, FloatType vs, FloatType rho, FloatType delta,
    FloatType epsilon, FloatType gamma, FloatType& c11, FloatType& c12,
    FloatType& c13, FloatType& c33, FloatType& c44, FloatType& c66)
{
  FloatType const rho_vp2 = rho * vp * vp;
  FloatType const rho_vs2 = rho * vs * vs;

  c33 = rho_vp2;
  c44 = rho_vs2;
  c11 = rho_vp2 * (FloatType(1.0) + FloatType(2.0) * epsilon);
  c66 = rho_vs2 * (FloatType(1.0) + FloatType(2.0) * gamma);

  FloatType const vp2_vs2 = vp * vp - vs * vs;
  FloatType const sqrt_arg =
      vp2_vs2 * vp2_vs2 + FloatType(2.0) * rho_vp2 * delta * vp2_vs2;
  c13 = rho * sqrt(sqrt_arg) - rho_vs2;
  c12 = c11 - FloatType(2.0) * c66;
}

/**
 * @brief Build VTI elasticity tensor in Voigt notation (no rotation)
 */
template <typename FloatType>
PROXY_HOST_DEVICE void buildVTITensor(FloatType c11, FloatType c12,
                                      FloatType c13, FloatType c33,
                                      FloatType c44, FloatType c66,
                                      FloatType C[6][6])
{
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) C[i][j] = FloatType(0.0);

  C[0][0] = c11;
  C[1][1] = c11;
  C[2][2] = c33;
  C[0][1] = c12;
  C[1][0] = c12;
  C[0][2] = c13;
  C[2][0] = c13;
  C[1][2] = c13;
  C[2][1] = c13;
  C[3][3] = c44;
  C[4][4] = c44;
  C[5][5] = c66;
}

/**
 * @brief Compute the full TTI elasticity tensor (VTI + rotation)
 */
template <typename FloatType>
PROXY_HOST_DEVICE void computeCTensor(FloatType vp, FloatType vs, FloatType rho,
                                      FloatType delta, FloatType epsilon,
                                      FloatType gamma, FloatType theta,
                                      FloatType phi, FloatType CTTI[6][6])
{
  FloatType CVTI[6][6] = {0.0};
  CVTI[0][0] = rho * vp * vp * (1.0 + 2.0 * epsilon);
  CVTI[1][1] = CVTI[0][0];
  CVTI[2][2] = rho * vp * vp;
  CVTI[3][3] = rho * vs * vs;
  CVTI[4][4] = CVTI[3][3];
  CVTI[5][5] = rho * vs * vs * (1.0 + 2.0 * gamma);
  CVTI[0][1] = CVTI[0][0] - 2.0 * CVTI[5][5];
  CVTI[1][0] = CVTI[0][1];

  FloatType vp2 = vp * vp;
  FloatType vs2 = vs * vs;
  FloatType diff = vp2 - vs2;
  CVTI[0][2] = rho * sqrt(diff * diff + 2.0 * vp2 * delta * diff) - rho * vs2;
  CVTI[1][2] = CVTI[0][2];
  CVTI[2][0] = CVTI[0][2];
  CVTI[2][1] = CVTI[0][2];

  constexpr FloatType PI = FloatType(3.14159265358979323846);
  FloatType theta_rad = theta * PI / FloatType(180.0);
  FloatType phi_rad = phi * PI / FloatType(180.0);

  FloatType ctheta = cos(theta_rad);
  FloatType stheta = sin(theta_rad);
  FloatType cphi = cos(phi_rad);
  FloatType sphi = sin(phi_rad);

  FloatType R[3][3];
  R[0][0] = ctheta * cphi;
  R[0][1] = ctheta * sphi;
  R[0][2] = -stheta;
  R[1][0] = -sphi;
  R[1][1] = cphi;
  R[1][2] = 0.0;
  R[2][0] = stheta * cphi;
  R[2][1] = stheta * sphi;
  R[2][2] = ctheta;

  FloatType M[6][6] = {0.0};
  M[0][0] = R[0][0] * R[0][0];
  M[0][1] = R[0][1] * R[0][1];
  M[0][2] = R[0][2] * R[0][2];
  M[1][0] = R[1][0] * R[1][0];
  M[1][1] = R[1][1] * R[1][1];
  M[1][2] = R[1][2] * R[1][2];
  M[2][0] = R[2][0] * R[2][0];
  M[2][1] = R[2][1] * R[2][1];
  M[2][2] = R[2][2] * R[2][2];
  M[0][3] = R[0][1] * R[0][2];
  M[0][4] = R[0][0] * R[0][2];
  M[0][5] = R[0][0] * R[0][1];
  M[1][3] = R[1][1] * R[1][2];
  M[1][4] = R[1][0] * R[1][2];
  M[1][5] = R[1][0] * R[1][1];
  M[2][3] = R[2][1] * R[2][2];
  M[2][4] = R[2][0] * R[2][2];
  M[2][5] = R[2][0] * R[2][1];
  M[3][0] = 2.0 * R[1][0] * R[2][0];
  M[3][1] = 2.0 * R[1][1] * R[2][1];
  M[3][2] = 2.0 * R[1][2] * R[2][2];
  M[3][3] = R[1][1] * R[2][2] + R[1][2] * R[2][1];
  M[3][4] = R[1][0] * R[2][2] + R[1][2] * R[2][0];
  M[3][5] = R[1][0] * R[2][1] + R[1][1] * R[2][0];
  M[4][0] = 2.0 * R[0][0] * R[2][0];
  M[4][1] = 2.0 * R[0][1] * R[2][1];
  M[4][2] = 2.0 * R[0][2] * R[2][2];
  M[4][3] = R[0][1] * R[2][2] + R[0][2] * R[2][1];
  M[4][4] = R[0][0] * R[2][2] + R[0][2] * R[2][0];
  M[4][5] = R[0][0] * R[2][1] + R[0][1] * R[2][0];
  M[5][0] = 2.0 * R[0][0] * R[1][0];
  M[5][1] = 2.0 * R[0][1] * R[1][1];
  M[5][2] = 2.0 * R[0][2] * R[1][2];
  M[5][3] = R[0][1] * R[1][2] + R[0][2] * R[1][1];
  M[5][4] = R[0][0] * R[1][2] + R[0][2] * R[1][0];
  M[5][5] = R[0][0] * R[1][1] + R[0][1] * R[1][0];

  FloatType temp[6][6] = {0.0};
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      for (int k = 0; k < 6; k++) temp[i][j] += M[i][k] * CVTI[k][j];

  for (int i = 0; i < 6; i++)
    for (int j = i; j < 6; j++)
    {
      CTTI[i][j] = 0.0;
      for (int k = 0; k < 6; k++) CTTI[i][j] += temp[i][k] * M[j][k];
      if (i != j) CTTI[j][i] = CTTI[i][j];
    }
}
#endif  // FUNTIDES_MODEL_MESH_IMPL_COMMON_INCLUDE_ELASTICITY_UTILS_H_
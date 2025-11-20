#ifndef ELASTICITY_UTILS_H
#define ELASTICITY_UTILS_H

#include <data_type.h>

/**
 * @brief Compute the elasticity tensor C from Thomsen parameters
 * @param vp P-wave velocity
 * @param vs S-wave velocity
 * @param rho Density
 * @param delta Thomsen parameter delta
 * @param epsilon Thomsen parameter epsilon
 * @param gamma Thomsen parameter gamma
 * @param theta Rotation angle theta (degrees)
 * @param phi Rotation angle phi (degrees)
 * @param[out] CTTI Output 6x6 elasticity tensor (Voigt notation)
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
    for (int k = 0; k < 6; k++)
      for (int j = 0; j < 6; j++) temp[i][j] += M[i][k] * CVTI[k][j];

  for (int i = 0; i < 6; i++)
    for (int j = i; j < 6; j++)
    {
      CTTI[i][j] = 0.0;
      for (int k = 0; k < 6; k++) CTTI[i][j] += temp[i][k] * M[j][k];
      if (i != j) CTTI[j][i] = CTTI[i][j];
    }
}

#endif  // ELASTICITY_UTILS_H

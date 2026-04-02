#pragma once

#include <cmath>
#include <tuple>

#include "macros.hpp"

template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE T determinant(T const (&m)[3][3])
{
  return +m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
         m[0][1] * (m[1][0] * m[2][2] - m[2][0] * m[1][2]) +
         m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
}

template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE typename T::value_type
determinant(T const& m)
{
  return +m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
         m(0, 1) * (m(1, 0) * m(2, 2) - m(2, 0) * m(1, 2)) +
         m(0, 2) * (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));
}

/**
 * @brief Calculates the linear index for support/quadrature points from ijk
 *   coordinates.
 * @param r order of polynomial approximation
 * @param i The index in the xi0 direction (0,r)
 * @param j The index in the xi1 direction (0,r)
 * @param k The index in the xi2 direction (0,r)
 * @return The linear index of the support/quadrature point (0-(r+1)^3)
 */
template <int ORDER>
static constexpr inline SEMKERNELS_HOST_DEVICE int linearIndex(const int i,
                                                               const int j,
                                                               const int k)
{
  return i + (ORDER + 1) * j + (ORDER + 1) * (ORDER + 1) * k;
}

/**
 * @brief Calculate the Cartesian/TensorProduct index given the linear index
 *   of a support point.
 * @param linearIndex The linear index of support point
 * @param r order of polynomial approximation
 * @param i0 The Cartesian index of the support point in the xi0 direction.
 * @param i1 The Cartesian index of the support point in the xi1 direction.
 * @param i2 The Cartesian index of the support point in the xi2 direction.
 */
template <int ORDER>
static constexpr inline SEMKERNELS_HOST_DEVICE std::tuple<int, int, int>
tripleIndex(int const linearIndex)
{
  return {(linearIndex % ((ORDER + 1) * (ORDER + 1))) % (ORDER + 1),
          (linearIndex % ((ORDER + 1) * (ORDER + 1))) / (ORDER + 1),
          (linearIndex / ((ORDER + 1) * (ORDER + 1)))};
}

template <typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[3])
{
  return B[0] * B[1] - B[2] * B[2];
}

template <typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[6])
{
  return B[0] * B[1] * B[2] + B[5] * B[4] * B[3] * 2 - B[0] * B[3] * B[3] -
         B[1] * B[4] * B[4] - B[2] * B[5] * B[5];
}

template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3], T const (&J)[3][3])
{
  Jinv[0][0] = J[1][1] * J[2][2] - J[1][2] * J[2][1];
  Jinv[0][1] = J[0][2] * J[2][1] - J[0][1] * J[2][2];
  Jinv[0][2] = J[0][1] * J[1][2] - J[0][2] * J[1][1];
  T const det =
      J[0][0] * Jinv[0][0] + J[1][0] * Jinv[0][1] + J[2][0] * Jinv[0][2];

  T const invDet = T(1) / det;

  Jinv[0][0] *= invDet;
  Jinv[0][1] *= invDet;
  Jinv[0][2] *= invDet;
  Jinv[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * invDet;
  Jinv[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * invDet;
  Jinv[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * invDet;
  Jinv[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * invDet;
  Jinv[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * invDet;
  Jinv[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * invDet;

  return det;
}

template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3])
{
  T const J[3][3] = {{Jinv[0][0], Jinv[0][1], Jinv[0][2]},
                     {Jinv[1][0], Jinv[1][1], Jinv[1][2]},
                     {Jinv[2][0], Jinv[2][1], Jinv[2][2]}};
  return invert3x3(Jinv, J);
}

/**
 * @brief Invert the symmetric matrix @p srcSymMatrix and store the result in @p
 * dstSymMatrix.
 * @param dstSymMatrix The 3x3 symmetric matrix to write the inverse to.
 * @param srcSymMatrix The 3x3 symmetric matrix to take the inverse of.
 * @return The determinant.
 * @note @p srcSymMatrix can contain integers but @p dstMatrix must contain
 * floating point values.
 */
template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE void symInvert(
    T (&dstSymMatrix)[6], T const (&srcSymMatrix)[6])
{
  dstSymMatrix[0] =
      srcSymMatrix[1] * srcSymMatrix[2] - srcSymMatrix[3] * srcSymMatrix[3];
  dstSymMatrix[5] =
      srcSymMatrix[4] * srcSymMatrix[3] - srcSymMatrix[5] * srcSymMatrix[2];
  dstSymMatrix[4] =
      srcSymMatrix[5] * srcSymMatrix[3] - srcSymMatrix[4] * srcSymMatrix[1];

  T det = srcSymMatrix[0] * dstSymMatrix[0] +
          srcSymMatrix[5] * dstSymMatrix[5] + srcSymMatrix[4] * dstSymMatrix[4];

  T const invDet = 1.0 / det;

  dstSymMatrix[0] *= invDet;
  dstSymMatrix[5] *= invDet;
  dstSymMatrix[4] *= invDet;
  dstSymMatrix[1] =
      (srcSymMatrix[0] * srcSymMatrix[2] - srcSymMatrix[4] * srcSymMatrix[4]) *
      invDet;
  dstSymMatrix[3] =
      (srcSymMatrix[5] * srcSymMatrix[4] - srcSymMatrix[0] * srcSymMatrix[3]) *
      invDet;
  dstSymMatrix[2] =
      (srcSymMatrix[0] * srcSymMatrix[1] - srcSymMatrix[5] * srcSymMatrix[5]) *
      invDet;
}

/**
 * @brief Invert the symmetric matrix @p symMatrix overwritting it.
 * @param symMatrix The 3x3 symmetric matrix to take the inverse of and
 * overwrite.
 * @return The determinant.
 * @note @p symMatrix can contain integers but @p dstMatrix must contain
 * floating point values.
 */
template <typename T>
static inline SEMKERNELS_HOST_DEVICE void symInvert(T (&symMatrix)[6])
{
  T temp[6];
  symInvert(temp, symMatrix);

  symMatrix[0] = temp[0];
  symMatrix[1] = temp[1];
  symMatrix[2] = temp[2];
  symMatrix[3] = temp[3];
  symMatrix[4] = temp[4];
  symMatrix[5] = temp[5];
}

template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE void computeB(T const (&J)[3][3],
                                                             T (&B)[6])
{
  B[0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]);
  B[1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]);
  B[2] = (J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]);
  B[3] = (J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]);
  B[4] = (J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]);
  B[5] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]);

  symInvert(B);
}

template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE void computeB(
    T const& J, typename T::value_type (&B)[6])
{
  B[0] = (J(0, 0) * J(0, 0) + J(1, 0) * J(1, 0) + J(2, 0) * J(2, 0));
  B[1] = (J(0, 1) * J(0, 1) + J(1, 1) * J(1, 1) + J(2, 1) * J(2, 1));
  B[2] = (J(0, 2) * J(0, 2) + J(1, 2) * J(1, 2) + J(2, 2) * J(2, 2));
  B[3] = (J(0, 1) * J(0, 2) + J(1, 1) * J(1, 2) + J(2, 1) * J(2, 2));
  B[4] = (J(0, 0) * J(0, 2) + J(1, 0) * J(1, 2) + J(2, 0) * J(2, 2));
  B[5] = (J(0, 0) * J(0, 1) + J(1, 0) * J(1, 1) + J(2, 0) * J(2, 1));
  symInvert(B);
}

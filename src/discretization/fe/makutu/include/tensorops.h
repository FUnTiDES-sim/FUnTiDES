#ifndef TENSOROPS_H_
#define TENSOROPS_H_

// #include <commonMacro.hpp>
/**
 * @brief Inverts a 3x3 matrix and returns its determinant.
 *        The matrix is modified in place.
 * @param J The 3x3 matrix to invert.
 * @return The determinant of the original matrix.
 */
template <typename T>
PROXY_HOST_DEVICE T invert3x3(T (&J)[3][3]);

template <typename T>
PROXY_HOST_DEVICE T invert3x3(T (&J)[3][3])
{
  // Compute the determinant
  T det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
          J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
          J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

  T invDet = 1.0 / det;

  T inv[3][3];

  inv[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * invDet;
  inv[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * invDet;
  inv[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * invDet;

  inv[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * invDet;
  inv[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * invDet;
  inv[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * invDet;

  inv[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * invDet;
  inv[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * invDet;
  inv[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * invDet;

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) J[i][j] = inv[i][j];

  return det;
}

/**
 * @brief Invert the source matrix @p J and store the result in @p Jinv.
 * @tparam Jinv The type of @p Jinv.
 * @tparam J The type of @p J.
 * @param Jinv The 3v3 matrix to write the inverse to.
 * @param srcMatrix The 3x3 matrix to take the inverse of.
 * @return The determinant.
 * @note @p srcMatrix can contain integers but @p dstMatrix must contain
 * floating point values.
 */
template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3], T (&J)[3][3]);

template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3], T (&J)[3][3])
{
  Jinv[0][0] = J[1][1] * J[2][2] - J[1][2] * J[2][1];
  Jinv[0][1] = J[0][2] * J[2][1] - J[0][1] * J[2][2];
  Jinv[0][2] = J[0][1] * J[1][2] - J[0][2] * J[1][1];

  auto const det =
      J[0][0] * Jinv[0][0] + J[1][0] * Jinv[0][1] + J[2][0] * Jinv[0][2];

  auto const invDet = T(1) / det;

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

// template< int N, typename T >
// PROXY_HOST_DEVICE
// T determinant(const T (&A)[N][N]);

// template<typename T>
// PROXY_HOST_DEVICE
// T determinant(const T (&A)[3][3])
// {
//   return
//       A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
//       A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
//       A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
// }

template <int N, typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[3]);

template <typename T>
PROXY_HOST_DEVICE T symDeterminant<2>(T (&B)[3])
{
  return B[0] * B[1] - B[2] * B[2];
}

template <int N, typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[6]);

template <typename T>
PROXY_HOST_DEVICE T symDeterminant<3>(T (&B)[6])
{
  return B[0] * B[1] * B[2] + B[5] * B[4] * B[3] * 2 - B[0] * B[3] * B[3] -
         B[1] * B[4] * B[4] - B[2] * B[5] * B[5];
}

/**
 * @brief Invert the symmetric matrix @p J and store the result in @p dst.
 * @tparam DST_SYM_MATRIX The type of @p dst.
 * @tparam SRC_SYM_MATRIX The type of @p J.
 * @param dst The 3x3 symmetric matrix to write the inverse to.
 * @param J The 3x3 symmetric matrix to take the inverse of.
 * @return The determinant.
 * @note @p J can contain integers but @p dstMatrix must contain floating point
 * values.
 */
template <typename T>
PROXY_HOST_DEVICE static auto symInvert(T (&dst)[6], T const (&J)[6])
{
  using FloatingPoint = std::decay_t<decltype(dst[0])>;

  dst[0] = J[1] * J[2] - J[3] * J[3];
  dst[5] = J[4] * J[3] - J[5] * J[2];
  dst[4] = J[5] * J[3] - J[4] * J[1];

  auto const det = J[0] * dst[0] + J[5] * dst[5] + J[4] * dst[4];
  FloatingPoint const invDet = FloatingPoint(1) / det;

  dst[0] *= invDet;
  dst[5] *= invDet;
  dst[4] *= invDet;
  dst[1] = (J[0] * J[2] - J[4] * J[4]) * invDet;
  dst[3] = (J[5] * J[4] - J[0] * J[3]) * invDet;
  dst[2] = (J[0] * J[1] - J[5] * J[5]) * invDet;

  return det;
}

template <typename T>
PROXY_HOST_DEVICE static auto symInvert(T (&J)[6])
{
  std::remove_reference_t<decltype(J[0])> temp[6];
  auto const det = symInvert(temp, J);
  // std::copy< 6 >( J, temp );
  for (int i = 0; i < 6; i++) J[i] = temp[i];

  return det;
}

#endif  // TENSOROPS_H_

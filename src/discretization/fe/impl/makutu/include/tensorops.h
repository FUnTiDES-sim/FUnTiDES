#ifndef TENSOROPS_H_
#define TENSOROPS_H_

/**
 * @file tensorops.h
 * @brief Small dense and symmetric 3x3 matrix helpers.
 * @note Included by no file, and duplicates helpers of mathUtilites.h with
 * different signatures. It relies on PROXY_HOST_DEVICE and std type traits
 * without including them.
 */

/**
 * @brief Inverts a 3x3 matrix in place.
 * @tparam T Floating-point scalar type.
 * @param[in,out] J The matrix; holds its inverse on return. Must be invertible:
 * a zero determinant is not detected.
 * @return The determinant of the input matrix.
 */
template <typename T>
PROXY_HOST_DEVICE T invert3x3(T (&J)[3][3]);

template <typename T>
PROXY_HOST_DEVICE T invert3x3(T (&J)[3][3]) {
  T det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
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
 * @brief Inverts a 3x3 matrix into another one.
 * @tparam T Floating-point scalar type.
 * @param[out] Jinv The inverse of @p J.
 * @param[in] J The matrix to invert (not modified, although not const). Must
 * be invertible: a zero determinant is not detected.
 * @return The determinant of @p J.
 */
template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3], T (&J)[3][3]);

template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3], T (&J)[3][3]) {
  Jinv[0][0] = J[1][1] * J[2][2] - J[1][2] * J[2][1];
  Jinv[0][1] = J[0][2] * J[2][1] - J[0][1] * J[2][2];
  Jinv[0][2] = J[0][1] * J[1][2] - J[0][2] * J[1][1];

  auto const det = J[0][0] * Jinv[0][0] + J[1][0] * Jinv[0][1] + J[2][0] * Jinv[0][2];

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

/**
 * @brief Determinant of a symmetric 2x2 matrix in Voigt storage.
 * @tparam N Matrix size; 2 here.
 * @tparam T Floating-point scalar type.
 * @param[in] B Matrix entries in the order (B00, B11, B01).
 * @return The determinant of @p B.
 * @note The definition below is a partial specialization of a function
 * template, which C++ does not allow: it would not compile if this header
 * were included.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <int N, typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[3]);

template <typename T>
PROXY_HOST_DEVICE T symDeterminant<2>(T (&B)[3]) {
  return B[0] * B[1] - B[2] * B[2];
}

/**
 * @brief Determinant of a symmetric 3x3 matrix in Voigt storage.
 * @tparam N Matrix size; 3 here.
 * @tparam T Floating-point scalar type.
 * @param[in] B Matrix entries in the order (B00, B11, B22, B12, B02, B01).
 * @return The determinant of @p B.
 * @note Same partial-specialization limitation as the 2x2 overload.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <int N, typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[6]);

template <typename T>
PROXY_HOST_DEVICE T symDeterminant<3>(T (&B)[6]) {
  return B[0] * B[1] * B[2] + B[5] * B[4] * B[3] * 2 - B[0] * B[3] * B[3] - B[1] * B[4] * B[4] - B[2] * B[5] * B[5];
}

/**
 * @brief Inverts a symmetric 3x3 matrix in Voigt storage.
 * @tparam T Floating-point scalar type.
 * @param[out] dst The inverse of @p J, in the same storage.
 * @param[in] J The matrix to invert. Must be invertible: a zero determinant is
 * not detected.
 * @return The determinant of @p J.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <typename T>
PROXY_HOST_DEVICE static auto symInvert(T (&dst)[6], T const (&J)[6]) {
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

/**
 * @brief Inverts a symmetric 3x3 matrix in Voigt storage, in place.
 * @tparam T Floating-point scalar type.
 * @param[in,out] J The matrix; holds its inverse on return.
 * @return The determinant of the input matrix.
 */
template <typename T>
PROXY_HOST_DEVICE static auto symInvert(T (&J)[6]) {
  std::remove_reference_t<decltype(J[0])> temp[6];
  auto const det = symInvert(temp, J);
  for (int i = 0; i < 6; i++) J[i] = temp[i];

  return det;
}

#endif  // TENSOROPS_H_

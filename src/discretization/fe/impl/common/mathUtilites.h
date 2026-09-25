/**
 * @file
 * @brief Small 3x3 matrix helpers and element index maps used by the
 * hexahedral discretization kernels.
 *
 * Symmetric matrices use the Voigt storage of docs/design.md.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
#pragma once

#include <cmath>
#include <tuple>

#include "macros.h"

/**
 * @brief Determinant of a 3x3 matrix stored as a C array, m[row][column].
 * @tparam T Scalar type of the matrix entries.
 */
template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE T determinant(T const (&m)[3][3]) {
  return +m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - m[0][1] * (m[1][0] * m[2][2] - m[2][0] * m[1][2]) +
         m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
}

/**
 * @brief Determinant of a 3x3 matrix accessed as m(row, column).
 * @tparam T Matrix type with operator()(int, int) and a value_type member.
 */
template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE typename T::value_type determinant(T const& m) {
  return +m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) - m(0, 1) * (m(1, 0) * m(2, 2) - m(2, 0) * m(1, 2)) +
         m(0, 2) * (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));
}

/**
 * @brief Element-local index of the node (i, j, k) of a hexahedron.
 * @tparam ORDER Polynomial order; i, j and k lie in [0, ORDER].
 * @param[in] i Index along the first parent axis (xi0).
 * @param[in] j Index along the second parent axis (xi1).
 * @param[in] k Index along the third parent axis (xi2).
 * @return i + (ORDER+1)*j + (ORDER+1)^2*k, in [0, (ORDER+1)^3).
 * @see docs/design.md, "Hexahedron local numbering".
 */
template <int ORDER>
static constexpr inline SEMKERNELS_HOST_DEVICE int linearIndex(const int i, const int j, const int k) {
  return i + (ORDER + 1) * j + (ORDER + 1) * (ORDER + 1) * k;
}

/**
 * @brief Inverse of linearIndex(): the (i, j, k) indices of an element-local
 * node.
 * @tparam ORDER Polynomial order.
 * @param[in] linearIndex Element-local node index, in [0, (ORDER+1)^3).
 * @return The tuple (i, j, k), each in [0, ORDER].
 * @see docs/design.md, "Hexahedron local numbering".
 */
template <int ORDER>
static constexpr inline SEMKERNELS_HOST_DEVICE std::tuple<int, int, int> tripleIndex(int const linearIndex) {
  return {(linearIndex % ((ORDER + 1) * (ORDER + 1))) % (ORDER + 1),
          (linearIndex % ((ORDER + 1) * (ORDER + 1))) / (ORDER + 1), (linearIndex / ((ORDER + 1) * (ORDER + 1)))};
}

/**
 * @brief Determinant of a symmetric 2x2 matrix in Voigt storage
 * (B00, B11, B01).
 * @tparam T Scalar type of the matrix entries.
 */
template <typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[3]) {
  return B[0] * B[1] - B[2] * B[2];
}

/**
 * @brief Determinant of a symmetric 3x3 matrix in Voigt storage.
 * @tparam T Scalar type of the matrix entries.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <typename T>
PROXY_HOST_DEVICE T symDeterminant(T (&B)[6]) {
  return B[0] * B[1] * B[2] + B[5] * B[4] * B[3] * 2 - B[0] * B[3] * B[3] - B[1] * B[4] * B[4] - B[2] * B[5] * B[5];
}

/**
 * @brief Inverts a 3x3 matrix.
 * @tparam T Scalar type of the matrix entries.
 * @param[out] Jinv Inverse of @p J; must not alias @p J.
 * @param[in] J Matrix to invert, J[row][column]. It must be invertible: a zero
 * determinant is not detected.
 * @return The determinant of @p J.
 */
template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3], T const (&J)[3][3]) {
  Jinv[0][0] = J[1][1] * J[2][2] - J[1][2] * J[2][1];
  Jinv[0][1] = J[0][2] * J[2][1] - J[0][1] * J[2][2];
  Jinv[0][2] = J[0][1] * J[1][2] - J[0][2] * J[1][1];
  T const det = J[0][0] * Jinv[0][0] + J[1][0] * Jinv[0][1] + J[2][0] * Jinv[0][2];

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

/**
 * @brief Inverts a 3x3 matrix in place.
 * @tparam T Scalar type of the matrix entries.
 * @param[in,out] Jinv Matrix to invert, replaced by its inverse. It must be
 * invertible: a zero determinant is not detected.
 * @return The determinant of the input matrix.
 */
template <typename T>
PROXY_HOST_DEVICE auto invert3x3(T (&Jinv)[3][3]) {
  T const J[3][3] = {
      {Jinv[0][0], Jinv[0][1], Jinv[0][2]}, {Jinv[1][0], Jinv[1][1], Jinv[1][2]}, {Jinv[2][0], Jinv[2][1], Jinv[2][2]}};
  return invert3x3(Jinv, J);
}

/**
 * @brief Inverts a symmetric 3x3 matrix in Voigt storage.
 * @tparam T Scalar type of the matrix entries.
 * @param[out] dstSymMatrix Inverse of @p srcSymMatrix; must not alias it.
 * @param[in] srcSymMatrix Matrix to invert. It must be invertible: a zero
 * determinant is not detected.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE void symInvert(T (&dstSymMatrix)[6], T const (&srcSymMatrix)[6]) {
  dstSymMatrix[0] = srcSymMatrix[1] * srcSymMatrix[2] - srcSymMatrix[3] * srcSymMatrix[3];
  dstSymMatrix[5] = srcSymMatrix[4] * srcSymMatrix[3] - srcSymMatrix[5] * srcSymMatrix[2];
  dstSymMatrix[4] = srcSymMatrix[5] * srcSymMatrix[3] - srcSymMatrix[4] * srcSymMatrix[1];

  T det = srcSymMatrix[0] * dstSymMatrix[0] + srcSymMatrix[5] * dstSymMatrix[5] + srcSymMatrix[4] * dstSymMatrix[4];

  T const invDet = 1.0 / det;

  dstSymMatrix[0] *= invDet;
  dstSymMatrix[5] *= invDet;
  dstSymMatrix[4] *= invDet;
  dstSymMatrix[1] = (srcSymMatrix[0] * srcSymMatrix[2] - srcSymMatrix[4] * srcSymMatrix[4]) * invDet;
  dstSymMatrix[3] = (srcSymMatrix[5] * srcSymMatrix[4] - srcSymMatrix[0] * srcSymMatrix[3]) * invDet;
  dstSymMatrix[2] = (srcSymMatrix[0] * srcSymMatrix[1] - srcSymMatrix[5] * srcSymMatrix[5]) * invDet;
}

/**
 * @brief Inverts a symmetric 3x3 matrix in Voigt storage, in place.
 * @tparam T Scalar type of the matrix entries.
 * @param[in,out] symMatrix Matrix to invert, replaced by its inverse. It must
 * be invertible: a zero determinant is not detected.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <typename T>
static inline SEMKERNELS_HOST_DEVICE void symInvert(T (&symMatrix)[6]) {
  T temp[6];
  symInvert(temp, symMatrix);

  symMatrix[0] = temp[0];
  symMatrix[1] = temp[1];
  symMatrix[2] = temp[2];
  symMatrix[3] = temp[3];
  symMatrix[4] = temp[4];
  symMatrix[5] = temp[5];
}

/**
 * @brief Computes B = (J^T J)^-1 from a Jacobian matrix.
 * @tparam T Scalar type of the matrix entries.
 * @param[in] J Jacobian matrix, J[row][column]. J^T J must be invertible.
 * @param[out] B (J^T J)^-1 in Voigt storage.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE void computeB(T const (&J)[3][3], T (&B)[6]) {
  B[0] = (J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0]);
  B[1] = (J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]);
  B[2] = (J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]);
  B[3] = (J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]);
  B[4] = (J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]);
  B[5] = (J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]);

  symInvert(B);
}

/**
 * @brief Computes B = (J^T J)^-1 from a Jacobian matrix accessed as J(row,
 * column).
 * @tparam T Matrix type with operator()(int, int) and a value_type member.
 * @param[in] J Jacobian matrix. J^T J must be invertible.
 * @param[out] B (J^T J)^-1 in Voigt storage.
 * @see docs/design.md, "Symmetric 3x3 matrices (Voigt storage)".
 */
template <typename T>
static constexpr inline SEMKERNELS_HOST_DEVICE void computeB(T const& J, typename T::value_type (&B)[6]) {
  B[0] = (J(0, 0) * J(0, 0) + J(1, 0) * J(1, 0) + J(2, 0) * J(2, 0));
  B[1] = (J(0, 1) * J(0, 1) + J(1, 1) * J(1, 1) + J(2, 1) * J(2, 1));
  B[2] = (J(0, 2) * J(0, 2) + J(1, 2) * J(1, 2) + J(2, 2) * J(2, 2));
  B[3] = (J(0, 1) * J(0, 2) + J(1, 1) * J(1, 2) + J(2, 1) * J(2, 2));
  B[4] = (J(0, 0) * J(0, 2) + J(1, 0) * J(1, 2) + J(2, 0) * J(2, 2));
  B[5] = (J(0, 0) * J(0, 1) + J(1, 0) * J(1, 1) + J(2, 0) * J(2, 1));
  symInvert(B);
}

#pragma once
// Common utilities shared by all Qk Hexahedron test group headers.
// TestedBases must be defined by the per-order driver before the group header is included.

#include <gtest/gtest.h>

#include <cmath>
#include <set>

#include "Qk_Hexahedron_Lagrange_GaussLobatto.h"

// Tolerances for floating point comparisons
// Adapted for float precision (real_t is float by default)
constexpr double TOL = 1.0e-10;           // For exact tests
constexpr double TOL_NUMERICAL = 1.0e-4;  // For numerical integration (accumulation errors)

#ifdef USE_DOUBLE
constexpr double TOL_MATRIX_INVERSION = 1e-12;  // For double precision
#else
constexpr double TOL_MATRIX_INVERSION = 1e-6;  // For float precision
#endif

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * @brief Creates an arbitrary cube [x0, x0+size]^3
 */
template <typename BASIS>
void createArbitraryCube(real_t (&X)[8][3], real_t x0, real_t y0, real_t z0, real_t size) {
  X[0][0] = x0;
  X[0][1] = y0;
  X[0][2] = z0;
  X[1][0] = x0 + size;
  X[1][1] = y0;
  X[1][2] = z0;
  X[2][0] = x0;
  X[2][1] = y0 + size;
  X[2][2] = z0;
  X[3][0] = x0 + size;
  X[3][1] = y0 + size;
  X[3][2] = z0;
  X[4][0] = x0;
  X[4][1] = y0;
  X[4][2] = z0 + size;
  X[5][0] = x0 + size;
  X[5][1] = y0;
  X[5][2] = z0 + size;
  X[6][0] = x0;
  X[6][1] = y0 + size;
  X[6][2] = z0 + size;
  X[7][0] = x0 + size;
  X[7][1] = y0 + size;
  X[7][2] = z0 + size;
}

/**
 * @brief Creates a reference hexahedron [-1,1]^3
 */
template <typename BASIS>
void createReferenceHex(real_t (&X)[8][3]) {
  createArbitraryCube<BASIS>(X, -1.0, -1.0, -1.0, 2.0);
}

/**
 * @brief Creates a unit hexahedron [0,1]^3
 */
template <typename BASIS>
void createUnitHex(real_t (&X)[8][3]) {
  createArbitraryCube<BASIS>(X, 0.0, 0.0, 0.0, 1.0);
}

/**
 * @brief Helper to compute matrix-matrix product
 */
inline void matMul3x3(real_t const (&A)[3][3], real_t const (&B)[3][3], real_t (&C)[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      C[i][j] = 0.0;
      for (int k = 0; k < 3; ++k) C[i][j] += A[i][k] * B[k][j];
    }
}

/**
 * @brief Check if matrix is identity
 */
inline bool isIdentity(real_t const (&M)[3][3], double tol = TOL) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      real_t expected = (i == j) ? 1.0 : 0.0;
      if (std::abs(M[i][j] - expected) > tol) return false;
    }
  return true;
}

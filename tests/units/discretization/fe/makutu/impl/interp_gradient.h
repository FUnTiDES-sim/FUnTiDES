#pragma once
#include "common.h"

template <typename QK_BASIS>
class InterpolationTest : public ::testing::Test {};

template <typename QK_BASIS>
class GradientTest : public ::testing::Test {};

template <typename QK_BASIS>
class BasisGradientTest : public ::testing::Test {};

TYPED_TEST_SUITE(InterpolationTest, TestedBases);
TYPED_TEST_SUITE(GradientTest, TestedBases);
TYPED_TEST_SUITE(BasisGradientTest, TestedBases);

// ============================================================================
// INTERPOLATION TESTS
// ============================================================================

TYPED_TEST(InterpolationTest, TrilinearInterpAtCorners) {
  using QK = TypeParam;

  real_t X[8][3];
  createUnitHex<QK>(X);

  struct CornerTest {
    real_t alpha, beta, gamma;
    int idx;
  };

  CornerTest corners[] = {{0.0, 0.0, 0.0, 0}, {1.0, 0.0, 0.0, 1}, {0.0, 1.0, 0.0, 2}, {1.0, 1.0, 0.0, 3},
                          {0.0, 0.0, 1.0, 4}, {1.0, 0.0, 1.0, 5}, {0.0, 1.0, 1.0, 6}, {1.0, 1.0, 1.0, 7}};

  for (int c = 0; c < 8; ++c) {
    real_t coords[3];
    QK::trilinearInterp(corners[c].alpha, corners[c].beta, corners[c].gamma, X, coords);

    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR(coords[i], X[corners[c].idx][i], TOL)
          << "Trilinear interpolation at corner " << corners[c].idx << " should match corner coordinates";
    }
  }
}

TYPED_TEST(InterpolationTest, TrilinearInterpAtCenter_ArbitraryCube) {
  using QK = TypeParam;

  real_t x0 = 3.5, y0 = -2.0, z0 = 1.5, size = 2.0;
  real_t X[8][3];
  createArbitraryCube<QK>(X, x0, y0, z0, size);

  real_t coords[3];
  QK::trilinearInterp(0.5, 0.5, 0.5, X, coords);

  real_t halfSize = size * 0.5f;
  real_t expectedCenter[3] = {x0 + halfSize, y0 + halfSize, z0 + halfSize};

  EXPECT_NEAR(coords[0], expectedCenter[0], TOL);
  EXPECT_NEAR(coords[1], expectedCenter[1], TOL);
  EXPECT_NEAR(coords[2], expectedCenter[2], TOL);
}

TYPED_TEST(InterpolationTest, ComputeLocalCoordsConsistency) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t Xmesh[8][3];
  createArbitraryCube<QK>(Xmesh, -1.5, 0.8, 2.2, 1.3);

  real_t X[numNodes][3];
  QK::computeLocalCoords(Xmesh, X);

  for (int k = 0; k < 8; ++k) {
    int nodeIdx = QK::meshIndexToLinearIndex3D(k);
    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR(X[nodeIdx][i], Xmesh[k][i], TOL) << "Corner node " << k << " should match mesh corner";
    }
  }
}

TYPED_TEST(InterpolationTest, InterpolationCoefficientsSum) {
  using QK = TypeParam;

  for (int q = 0; q < QK::num1dNodes; ++q) {
    real_t c0 = QK::interpolationCoord(q, 0);
    real_t c1 = QK::interpolationCoord(q, 1);

    EXPECT_NEAR(c0 + c1, 1.0, TOL) << "Interpolation coefficients should sum to 1 at quadrature point " << q;

    EXPECT_GE(c0, 0.0) << "Interpolation coefficient should be non-negative";
    EXPECT_LE(c0, 1.0) << "Interpolation coefficient should be <= 1";
    EXPECT_GE(c1, 0.0) << "Interpolation coefficient should be non-negative";
    EXPECT_LE(c1, 1.0) << "Interpolation coefficient should be <= 1";
  }
}

TYPED_TEST(InterpolationTest, InterpolationCoefficientsAtBoundaries) {
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  real_t c0_first = QK::interpolationCoord(0, 0);
  real_t c1_first = QK::interpolationCoord(0, 1);

  EXPECT_NEAR(c0_first, 1.0, TOL) << "At first node, k=0 coefficient should be 1";
  EXPECT_NEAR(c1_first, 0.0, TOL) << "At first node, k=1 coefficient should be 0";

  real_t c0_last = QK::interpolationCoord(num1d - 1, 0);
  real_t c1_last = QK::interpolationCoord(num1d - 1, 1);

  EXPECT_NEAR(c0_last, 0.0, TOL) << "At last node, k=0 coefficient should be 0";
  EXPECT_NEAR(c1_last, 1.0, TOL) << "At last node, k=1 coefficient should be 1";
}

// ============================================================================
// GRADIENT OPERATOR TESTS
// ============================================================================

TYPED_TEST(GradientTest, GradientOfLinearFieldIsConstant_ArbitraryCube) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -2.0, 1.5, 0.3, 1.7);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8) {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  } else {
    QK::computeLocalCoords(X, Xfull);
  }

  real_t u[numNodes][3];
  for (int i = 0; i < numNodes; ++i) {
    u[i][0] = Xfull[i][0];
    u[i][1] = 2.0 * Xfull[i][1];
    u[i][2] = 3.0 * Xfull[i][2];
  }

  real_t expectedGrad[3][3] = {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}};

  int numTestPoints = (QK::numQuadraturePoints < 5) ? QK::numQuadraturePoints : 5;
  for (int q = 0; q < numTestPoints; ++q) {
    real_t J[3][3] = {{0}};
    real_t invJ[3][3] = {{0}};
    real_t grad[3][3] = {{0}};

    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);
    QK::jacobianTransformation(qa, qb, qc, X, J);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) invJ[i][j] = J[i][j];
    invert3x3(invJ);

    QK::gradient(q, invJ, u, grad);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        EXPECT_NEAR(grad[i][j], expectedGrad[i][j], TOL_NUMERICAL)
            << "Gradient of linear field should be constant at all quadrature points";
      }
    }
  }
}

TYPED_TEST(GradientTest, SymmetricGradientIsSymmetric) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.5, -1.2, 3.0, 2.1);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8) {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  } else {
    QK::computeLocalCoords(X, Xfull);
  }

  real_t u[numNodes][3];
  for (int i = 0; i < numNodes; ++i) {
    u[i][0] = Xfull[i][0] + 0.5 * Xfull[i][1];
    u[i][1] = Xfull[i][1] + 0.3 * Xfull[i][2];
    u[i][2] = Xfull[i][2] + 0.2 * Xfull[i][0];
  }

  int numTestPoints = (QK::numQuadraturePoints < 10) ? QK::numQuadraturePoints : 10;
  for (int q = 0; q < numTestPoints; ++q) {
    real_t J[3][3] = {{0}};
    real_t invJ[3][3] = {{0}};

    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);
    QK::jacobianTransformation(qa, qb, qc, X, J);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) invJ[i][j] = J[i][j];
    invert3x3(invJ);

    real_t symGrad[6] = {0};
    QK::symmetricGradient(q, invJ, u, symGrad);

    EXPECT_TRUE(std::isfinite(symGrad[0]));
    EXPECT_TRUE(std::isfinite(symGrad[1]));
    EXPECT_TRUE(std::isfinite(symGrad[2]));
    EXPECT_TRUE(std::isfinite(symGrad[3]));
    EXPECT_TRUE(std::isfinite(symGrad[4]));
    EXPECT_TRUE(std::isfinite(symGrad[5]));
  }
}

// ============================================================================
// BASIS GRADIENT TESTS
// ============================================================================

TYPED_TEST(BasisGradientTest, BasisGradientSymmetryProperty) {
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int q = 0; q < num1d; ++q) {
    for (int p = QK::halfNodes + 1; p < num1d; ++p) {
      real_t g1 = QK::basisGradientAt(q, p);
      real_t g2 = QK::basisGradientAt(num1d - 1 - q, num1d - 1 - p);

      EXPECT_NEAR(g1, -g2, TOL) << "Basis gradient should satisfy symmetry property: " << "grad(" << q << "," << p
                                << ") = -grad(" << (num1d - 1 - q) << "," << (num1d - 1 - p) << ")";
    }
  }
}

TYPED_TEST(BasisGradientTest, BasisGradientZeroAtSameNode) {
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int q = 1; q < num1d - 1; ++q) {
    real_t grad = QK::basisGradientAt(q, q);
    EXPECT_NEAR(grad, 0.0, TOL) << "Basis gradient should be zero at its own interior node";
  }
}

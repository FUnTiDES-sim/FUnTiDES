#pragma once
#include "common.h"

template <typename QK_BASIS>
class FaceOperationsTest : public ::testing::Test {};

template <typename QK_BASIS>
class InterfaceFluxTest : public ::testing::Test {};

template <typename QK_BASIS>
class VirtualMethodTest : public ::testing::Test {};

TYPED_TEST_SUITE(FaceOperationsTest, TestedBases);
TYPED_TEST_SUITE(InterfaceFluxTest, TestedBases);
TYPED_TEST_SUITE(VirtualMethodTest, TestedBases);

// ============================================================================
// 2D FACE OPERATIONS TESTS
// ============================================================================

TYPED_TEST(FaceOperationsTest, Jacobian2DRankTwo) {
  using QK = TypeParam;

  real_t X[4][3];
  X[0][0] = 0.0;
  X[0][1] = 0.0;
  X[0][2] = 0.0;
  X[1][0] = 1.0;
  X[1][1] = 0.0;
  X[1][2] = 0.0;
  X[2][0] = 0.0;
  X[2][1] = 1.0;
  X[2][2] = 0.0;
  X[3][0] = 1.0;
  X[3][1] = 1.0;
  X[3][2] = 0.0;

  int qa = QK::num1dNodes / 2;
  int qb = QK::num1dNodes / 2;

  real_t J[3][2] = {{0}};
  QK::jacobianTransformation2d(qa, qb, X, J);

  EXPECT_NEAR(J[0][0], 0.5, TOL) << "J[0][0] should be ~0.5 for unit square";
  EXPECT_NEAR(J[1][0], 0.0, TOL) << "J[1][0] should be ~0";
  EXPECT_NEAR(J[2][0], 0.0, TOL) << "J[2][0] should be ~0";

  EXPECT_NEAR(J[0][1], 0.0, TOL) << "J[0][1] should be ~0";
  EXPECT_NEAR(J[1][1], 0.5, TOL) << "J[1][1] should be ~0.5 for unit square";
  EXPECT_NEAR(J[2][1], 0.0, TOL) << "J[2][1] should be ~0";
}

TYPED_TEST(FaceOperationsTest, DampingTermPositive_ArbitrarySquare) {
  using QK = TypeParam;
  constexpr int numNodesPerFace = QK::numNodesPerFace;

  real_t x0 = 2.5, y0 = -1.0, size = 1.5;
  real_t X[4][3];
  X[0][0] = x0;
  X[0][1] = y0;
  X[0][2] = 0.0;
  X[1][0] = x0 + size;
  X[1][1] = y0;
  X[1][2] = 0.0;
  X[2][0] = x0;
  X[2][1] = y0 + size;
  X[2][2] = 0.0;
  X[3][0] = x0 + size;
  X[3][1] = y0 + size;
  X[3][2] = 0.0;

  real_t totalDamping = 0.0;
  real_t expectedArea = size * size;

  for (int q = 0; q < numNodesPerFace; ++q) {
    real_t damping = QK::computeDampingTerm(q, X);

    EXPECT_GT(damping, 0.0) << "Damping term should be positive at node " << q;

    totalDamping += damping;
  }

  EXPECT_NEAR(totalDamping, expectedArea, TOL_NUMERICAL) << "Sum of damping terms should equal face area";
}

TYPED_TEST(FaceOperationsTest, DampingTermScaling) {
  using QK = TypeParam;

  real_t X1[4][3], X2[4][3];
  X1[0][0] = 0.0;
  X1[0][1] = 0.0;
  X1[0][2] = 0.0;
  X1[1][0] = 1.0;
  X1[1][1] = 0.0;
  X1[1][2] = 0.0;
  X1[2][0] = 0.0;
  X1[2][1] = 1.0;
  X1[2][2] = 0.0;
  X1[3][0] = 1.0;
  X1[3][1] = 1.0;
  X1[3][2] = 0.0;

  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 3; ++i) X2[k][i] = 2.0 * X1[k][i];

  int q = QK::numNodesPerFace / 2;

  real_t d1 = QK::computeDampingTerm(q, X1);
  real_t d2 = QK::computeDampingTerm(q, X2);

  EXPECT_NEAR(d2 / d1, 4.0, TOL_NUMERICAL) << "Damping term should scale quadratically with element size";
}

// ============================================================================
// INTERFACE FLUX TESTS
// ============================================================================

TYPED_TEST(InterfaceFluxTest, InterfaceFluxIsZero) {
  using QK = TypeParam;
  constexpr int numNodesPerFace = QK::numNodesPerFace;

  real_t X8[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

  real_t X[4][3];
  X[0][0] = 0.0;
  X[0][1] = 0.0;
  X[0][2] = 0.0;
  X[1][0] = 1.0;
  X[1][1] = 0.0;
  X[1][2] = 0.0;
  X[2][0] = 0.0;
  X[2][1] = 1.0;
  X[2][2] = 0.0;
  X[3][0] = 1.0;
  X[3][1] = 1.0;
  X[3][2] = 0.0;

  real_t CKK[numNodesPerFace][numNodesPerFace][3] = {{{0}}};
  for (int faceId = 0; faceId < 6; ++faceId) {
    QK::computeInterfaceFluxTerm(X, X8, faceId, [&](int i, int j, int k, real_t Cijk) { CKK[i][j][k] += Cijk; });

    real_t SumGrad[3] = {0};
    real_t Sum;
    for (int i = 0; i < numNodesPerFace; ++i) {
      for (int j = 0; j < numNodesPerFace; ++j) {
        for (int k = 0; k < 3; ++k) {
          SumGrad[k] += CKK[i][j][k];
        }
      }
    }
    Sum = SumGrad[0] + SumGrad[1] + SumGrad[2];
    EXPECT_NEAR(Sum, 0.0, TOL_NUMERICAL) << "Sum of all CKK coefficients should be zero";
  }
}

// The interface flux operator must reproduce the EXACT normal derivative of a field, not only its
// tangential part. A field varying only along the face normal is the case that isolates this: such
// a field is constant over the face dofs, so any operator restricted to face dofs sees a constant
// and returns zero by partition of unity -- which is precisely what InterfaceFluxIsZero asserts.
// Reproducing it requires the dofs of the line perpendicular to the face, hence the funcNormal
// channel.
//
// Unit cube, p = the physical coordinate running along the face normal. Then grad(p) = e_kDir and
//   sum_i p_i C_i,j,kDir = int_F (dp/dx_kDir) phi_j = int_F phi_j = computeDampingTerm(j, X).
// The face normal's own sign is applied by the caller, not by the operator, so it does not enter
// the expected value.
TYPED_TEST(InterfaceFluxTest, ReproducesNormalDerivativeOfLinearField) {
  using QK = TypeParam;
  constexpr int numNodesPerFace = QK::numNodesPerFace;
  constexpr int num1dNodes = QK::num1dNodes;

  real_t X8[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  real_t X[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};

  // Nodal values of that coordinate along one direction, in [0, 1] on the unit cube.
  real_t coord1d[num1dNodes];
  for (int m = 0; m < num1dNodes; ++m) coord1d[m] = QK::interpolationCoord(m, 1);

  for (int faceId = 0; faceId < 6; ++faceId) {
    const int kDir = faceId / 2;
    const int kQFixed = (faceId % 2 == 0) ? 0 : num1dNodes - 1;

    real_t acc[numNodesPerFace] = {0};
    for (int q = 0; q < numNodesPerFace; ++q) {
      QK::computeInterfaceFluxTermAt(
          q, X, X8, faceId,
          // Tangential channel: i is a face dof, where p takes its face-constant value.
          [&](int i, int j, int k, real_t Cijk) {
            (void)i;
            if (k == kDir) acc[j] += coord1d[kQFixed] * Cijk;
          },
          // Normal channel: m is the depth along the line through face dof j, where p varies.
          [&](int m, int j, int k, real_t Cijk) {
            if (k == kDir) acc[j] += coord1d[m] * Cijk;
          });
    }

    for (int j = 0; j < numNodesPerFace; ++j) {
      EXPECT_NEAR(acc[j], QK::computeDampingTerm(j, X), TOL_NUMERICAL) << "faceId=" << faceId << ", face dof " << j;
    }
  }
}

// ============================================================================
// VIRTUAL METHOD TESTS
// ============================================================================

TYPED_TEST(VirtualMethodTest, GetNumQuadraturePointsMatchesStatic) {
  using QK = TypeParam;
  QK elem;
  EXPECT_EQ(elem.getNumQuadraturePoints(), QK::numQuadraturePoints);
  EXPECT_EQ(elem.getNumQuadraturePoints(), QK::num1dNodes * QK::num1dNodes * QK::num1dNodes);
}

TYPED_TEST(VirtualMethodTest, GetNumSupportPointsMatchesStatic) {
  using QK = TypeParam;
  QK elem;
  EXPECT_EQ(elem.getNumSupportPoints(), QK::numNodes);
  EXPECT_EQ(elem.getNumSupportPoints(), QK::num1dNodes * QK::num1dNodes * QK::num1dNodes);
}

TYPED_TEST(VirtualMethodTest, GetMaxSupportPointsMatchesStatic) {
  using QK = TypeParam;
  const QK elem;
  EXPECT_EQ(elem.getMaxSupportPoints(), QK::maxSupportPoints);
  EXPECT_EQ(elem.getMaxSupportPoints(), QK::numNodes);
}

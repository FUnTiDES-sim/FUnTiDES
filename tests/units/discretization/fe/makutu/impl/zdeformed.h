#pragma once
// Tests of the helpers for elements deformed along z only: they must agree with the
// generic geometry and stiffness routines on such an element.

#include "common.h"

template <typename T>
class ZDeformedTest : public ::testing::Test {};

TYPED_TEST_SUITE(ZDeformedTest, TestedBases);

/**
 * @brief Box whose 8 vertices are shifted along z by different amounts.
 * @param[out] X Vertex coordinates, vertex k = ka + 2 kb + 4 kc.
 * @param[out] Z z coordinate of the 8 vertices.
 * @param[out] J00 Half the element size along x.
 * @param[out] J11 Half the element size along y.
 */
template <typename BASIS>
void createZDeformedHex(real_t (&X)[8][3], real_t (&Z)[8], real_t &J00, real_t &J11) {
  createArbitraryCube<BASIS>(X, -3.0, 1.5, -0.5, 1.8);
  for (int k = 0; k < 8; ++k) {
    X[k][1] = real_t(1.5) + real_t(1.25) * (X[k][1] - real_t(1.5));
    X[k][2] += real_t(0.07) * static_cast<real_t>(k) - real_t(0.2) * static_cast<real_t>(k % 3);
    Z[k] = X[k][2];
  }
  J00 = real_t(0.5) * (X[1][0] - X[0][0]);
  J11 = real_t(0.5) * (X[2][1] - X[0][1]);
}

/// Largest deviation of invJ * J from the identity.
inline real_t identityDefect(real_t const (&invJ)[3][3], real_t const (&J)[3][3]) {
  real_t P[3][3];
  matMul3x3(invJ, J, P);
  real_t defect = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) defect = std::fmax(defect, std::fabs(P[i][j] - (i == j ? real_t(1) : real_t(0))));
  return defect;
}

inline real_t det3x3(real_t const (&J)[3][3]) {
  return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
         J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

TYPED_TEST(ZDeformedTest, BasisTableMatchesBasis) {
  using QK = TypeParam;
  constexpr int n = QK::num1dNodes;

  EXPECT_EQ(QK::kBasisTableSize, n * n + 2 * n);
  for (int q = 0; q < n; ++q) {
    for (int p = 0; p < n; ++p) EXPECT_NEAR(QK::basisTableEntry(q * n + p), QK::basisGradientAt(q, p), TOL);
    EXPECT_NEAR(QK::basisTableEntry(n * n + q), static_cast<real_t>(QK::BasisType::weight(q)), TOL);
    EXPECT_NEAR(QK::basisTableEntry(n * n + n + q), QK::interpolationCoord(q, 1), TOL);
  }
}

TYPED_TEST(ZDeformedTest, TrilinearJacobianMatchesJacobianTransformation) {
  using QK = TypeParam;

  real_t X[8][3], Z[8], J00, J11;
  createZDeformedHex<QK>(X, Z, J00, J11);

  for (int q = 0; q < QK::numQuadraturePoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J_ref[3][3] = {{0}};
    QK::jacobianTransformation(qa, qb, qc, X, J_ref);

    real_t const alpha[3] = {QK::interpolationCoord(qa, 1), QK::interpolationCoord(qb, 1),
                             QK::interpolationCoord(qc, 1)};
    real_t J[3][3];
    QK::trilinearJacobian(alpha, X, J);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) EXPECT_NEAR(J[i][j], J_ref[i][j], TOL_NUMERICAL) << "point " << q;
  }
}

TYPED_TEST(ZDeformedTest, InverseJacobianMatchesFullJacobian) {
  using QK = TypeParam;

  real_t X[8][3], Z[8], J00, J11;
  createZDeformedHex<QK>(X, Z, J00, J11);
  real_t dZ[3][4];
  QK::zDeformedHalfDifferences(Z, dZ);

  for (int q = 0; q < QK::numQuadraturePoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J[3][3] = {{0}};
    QK::jacobianTransformation(qa, qb, qc, X, J);
    real_t const det_ref = det3x3(J);

    real_t invJ[3][3];
    real_t const det = QK::invJacobianZDeformed(qa, qb, qc, J00, J11, Z, invJ);
    EXPECT_LT(identityDefect(invJ, J), TOL_NUMERICAL) << "point " << q;
    EXPECT_NEAR(det, det_ref, TOL_NUMERICAL * std::fabs(det_ref)) << "point " << q;

    real_t const alpha[3] = {QK::interpolationCoord(qa, 1), QK::interpolationCoord(qb, 1),
                             QK::interpolationCoord(qc, 1)};
    real_t invJ_bilinear[3][3];
    real_t const det_bilinear = QK::invJacobianZDeformedBilinear(alpha, J00, J11, dZ, invJ_bilinear);
    EXPECT_LT(identityDefect(invJ_bilinear, J), TOL_NUMERICAL) << "point " << q;
    EXPECT_NEAR(det_bilinear, det_ref, TOL_NUMERICAL * std::fabs(det_ref)) << "point " << q;
  }
}

TYPED_TEST(ZDeformedTest, BMatrixMatchesGenericBMatrix) {
  using QK = TypeParam;

  real_t X[8][3], Z[8], J00, J11;
  createZDeformedHex<QK>(X, Z, J00, J11);

  for (int q = 0; q < QK::numQuadraturePoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J[3][3] = {{0}};
    real_t B_ref[6] = {0};
    QK::computeBMatrix(qa, qb, qc, X, J, B_ref);

    real_t B[6];
    QK::computeBMatrixZDeformed(qa, qb, qc, J00, J11, Z, B);

    for (int i = 0; i < 6; ++i)
      EXPECT_NEAR(B[i], B_ref[i], TOL_NUMERICAL * std::fmax(real_t(1), std::fabs(B_ref[i])))
          << "point " << q << ", component " << i;
  }
}

TYPED_TEST(ZDeformedTest, ScatterSumFactMatchesGenericSumFact) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3], Z[8], J00, J11;
  createZDeformedHex<QK>(X, Z, J00, J11);

  real_t p[numNodes];
  for (int i = 0; i < numNodes; ++i) p[i] = std::sin(static_cast<real_t>(i));

  real_t f_ref[numNodes] = {0};
  QK::computeStiffnessTermSumFact(X, p, f_ref, [](int, int, int) { return real_t(1); });

  real_t f[numNodes];
  QK::computeStiffnessTermSumFactZDeformedScatter(J00, J11, Z, p, f, [](int) { return real_t(1); });

  for (int i = 0; i < numNodes; ++i) EXPECT_NEAR(f[i], f_ref[i], TOL_NUMERICAL) << "node " << i;
}

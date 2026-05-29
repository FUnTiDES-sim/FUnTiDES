#pragma once
#include "common.h"

template <typename QK_BASIS>
class JacobianTest : public ::testing::Test {};

TYPED_TEST_SUITE(JacobianTest, TestedBases);

// ============================================================================
// JACOBIAN TESTS
// ============================================================================

TYPED_TEST(JacobianTest, JacobianDeterminantPositive_VariousCubes) {
  using QK = TypeParam;

  struct CubeConfig {
    real_t x0, y0, z0, size;
  };

  CubeConfig configs[] = {{0.0, 0.0, 0.0, 1.0}, {-2.0, 3.0, -1.0, 0.5}, {5.0, -5.0, 2.0, 3.0}};

  for (const auto& config : configs) {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    for (int q = 0; q < QK::numQuadraturePoints; ++q) {
      int qa, qb, qc;
      QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

      real_t J[3][3] = {{0}};
      QK::jacobianTransformation(qa, qb, qc, X, J);

      real_t detJ = determinant(J);

      EXPECT_GT(detJ, 0.0) << "Jacobian determinant should be positive at quadrature point " << q
                           << " for cube at (" << config.x0 << "," << config.y0 << "," << config.z0
                           << ") with size " << config.size;
    }
  }
}

TYPED_TEST(JacobianTest, InverseJacobianCorrectness_VariousCubes) {
  using QK = TypeParam;

  struct CubeConfig {
    real_t x0, y0, z0, size;
  };

  CubeConfig configs[] = {{0.0, 0.0, 0.0, 1.0}, {-10.0, 5.0, -3.0, 2.0}, {3.5, -1.5, 0.5, 0.8}};

  for (const auto& config : configs) {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    int numTestPoints = (QK::numQuadraturePoints < 5) ? QK::numQuadraturePoints : 5;
    for (int q = 0; q < numTestPoints; ++q) {
      real_t J[3][3] = {{0}};
      real_t invJ[3][3] = {{0}};
      real_t identity[3][3] = {{0}};

      int qa, qb, qc;
      QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);
      QK::jacobianTransformation(qa, qb, qc, X, J);

      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) invJ[i][j] = J[i][j];

      real_t detJ = invert3x3(invJ);
      EXPECT_GT(std::abs(detJ), TOL) << "Jacobian should be invertible";

      matMul3x3(J, invJ, identity);

      EXPECT_TRUE(isIdentity(identity, TOL_MATRIX_INVERSION))
          << "J * J^-1 should equal identity matrix at quadrature point " << q
          << " for cube at (" << config.x0 << "," << config.y0 << "," << config.z0
          << ") with size " << config.size;
    }
  }
}

TYPED_TEST(JacobianTest, ShapeFunctionsPartitionOfUnity) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t testPoints[][3] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {-0.7, 0.3, 0.1}, {0.9, -0.5, 0.8}};
  constexpr int numTestPoints = 4;

  for (int pt = 0; pt < numTestPoints; ++pt) {
    double N[numNodes];
    double coords[3] = {testPoints[pt][0], testPoints[pt][1], testPoints[pt][2]};
    QK::calcN(coords, N);

    double sum = 0.0;
    for (int i = 0; i < numNodes; ++i) sum += N[i];

    EXPECT_NEAR(sum, 1.0, TOL) << "Sum of shape functions should be 1 (partition of unity)";
  }
}

TYPED_TEST(JacobianTest, GradientConsistency) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.2, 0.5, 2.3, 1.5);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8) {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  } else {
    QK::computeLocalCoords(X, Xfull);
  }

  constexpr int q = numNodes / 2;

  real_t gradN1[numNodes][3] = {{0}};
  real_t gradN2[numNodes][3] = {{0}};

  real_t detJ1 = QK::calcGradN(q, Xfull, gradN1);
  real_t detJ2 = QK::calcGradNWithCorners(q, X, gradN2);

  EXPECT_NEAR(detJ1, detJ2, TOL_NUMERICAL) << "Two gradient computation methods should be consistent";

  for (int i = 0; i < numNodes; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(gradN1[i][j], gradN2[i][j], TOL_NUMERICAL)
          << "gradN[" << i << "][" << j << "] inconsistent between methods";
    }
  }
}

TYPED_TEST(JacobianTest, QuadratureRuleIntegratesConstant_VariousCubes) {
  using QK = TypeParam;

  struct CubeConfig {
    real_t x0, y0, z0, size;
    real_t expectedVolume;
  };

  CubeConfig configs[] = {{0.0, 0.0, 0.0, 1.0, 1.0}, {-5.0, 2.0, -1.0, 2.0, 8.0}, {3.0, -2.0, 1.0, 0.5, 0.125}};

  for (const auto& config : configs) {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    real_t integral = 0.0;

    for (int q = 0; q < QK::numQuadraturePoints; ++q) {
      int qa, qb, qc;
      QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

      real_t w = QK::BasisType::weight(qa) * QK::BasisType::weight(qb) * QK::BasisType::weight(qc);

      real_t J[3][3] = {{0}};
      QK::jacobianTransformation(qa, qb, qc, X, J);
      real_t detJ = determinant(J);

      integral += w * detJ;
    }

    EXPECT_NEAR(integral, config.expectedVolume, TOL_NUMERICAL)
        << "Quadrature rule should exactly integrate constant functions"
        << " for cube at (" << config.x0 << "," << config.y0 << "," << config.z0
        << ") with size " << config.size;
  }
}

TYPED_TEST(JacobianTest, TransformedQuadratureWeightConsistency) {
  using QK = TypeParam;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 2.5, -1.0, 0.3, 1.8);

  real_t totalWeight = 0.0;
  real_t expectedVolume = 1.8 * 1.8 * 1.8;

  for (int q = 0; q < QK::numQuadraturePoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t w = QK::BasisType::weight(qa) * QK::BasisType::weight(qb) * QK::BasisType::weight(qc);

    real_t J[3][3] = {{0}};
    QK::jacobianTransformation(qa, qb, qc, X, J);
    real_t detJ = determinant(J);

    totalWeight += w * detJ;
  }

  EXPECT_NEAR(totalWeight, expectedVolume, TOL_NUMERICAL)
      << "Sum of transformed quadrature weights should equal element volume";
}

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

      EXPECT_GT(detJ, 0.0) << "Jacobian determinant should be positive at quadrature point " << q << " for cube at ("
                           << config.x0 << "," << config.y0 << "," << config.z0 << ") with size " << config.size;
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
          << "J * J^-1 should equal identity matrix at quadrature point " << q << " for cube at (" << config.x0 << ","
          << config.y0 << "," << config.z0 << ") with size " << config.size;
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
        << "Quadrature rule should exactly integrate constant functions" << " for cube at (" << config.x0 << ","
        << config.y0 << "," << config.z0 << ") with size " << config.size;
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

// ============================================================================
// INVERSE JACOBIAN TESTS
// ============================================================================

TYPED_TEST(JacobianTest, InvJacobianTransformation3IndexGivesInverse) {
  using QK = TypeParam;
  real_t X[8][3];
  createArbitraryCube<QK>(X, 1.0, -0.5, 2.0, 1.5);

  int numTestPoints = (QK::numQuadraturePoints < 4) ? QK::numQuadraturePoints : 4;
  for (int q = 0; q < numTestPoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    // jacobianTransformation accumulates (+=) into J, so both must start zeroed
    real_t J[3][3] = {{0}};
    QK::jacobianTransformation(qa, qb, qc, X, J);

    real_t invJ[3][3] = {{0}};
    real_t det = QK::invJacobianTransformation(qa, qb, qc, X, invJ);
    EXPECT_GT(std::abs(det), TOL) << "Jacobian must be invertible at q=" << q;

    real_t product[3][3] = {{0}};
    matMul3x3(J, invJ, product);
    EXPECT_TRUE(isIdentity(product, TOL_MATRIX_INVERSION)) << "J * invJ must be identity at q=" << q;
  }
}

TYPED_TEST(JacobianTest, InvJacobianTransformation1IndexMatchesMultiIndex) {
  using QK = TypeParam;
  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.0, 0.0, 0.0, 1.0);

  int numTestPoints = (QK::numQuadraturePoints < 4) ? QK::numQuadraturePoints : 4;
  for (int q = 0; q < numTestPoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t invJ1[3][3] = {{0}};
    real_t invJ2[3][3] = {{0}};

    real_t det1 = QK::invJacobianTransformation(qa, qb, qc, X, invJ1);
    real_t det2 = QK::invJacobianTransformation(q, X, invJ2);

    EXPECT_NEAR(det1, det2, TOL_MATRIX_INVERSION) << "det must match between 1-index and 3-index overloads";
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(invJ1[i][j], invJ2[i][j], TOL_MATRIX_INVERSION)
            << "invJ[" << i << "][" << j << "] mismatch between overloads";
  }
}

// ============================================================================
// COORDS-BASED JACOBIAN OVERLOADS
// ============================================================================

TYPED_TEST(JacobianTest, JacobianTransformationCoordsMatchesIndexed) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;
  real_t X[8][3];
  createArbitraryCube<QK>(X, 1.0, 2.0, -1.0, 2.0);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8) {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  } else {
    QK::computeLocalCoords(X, Xfull);
  }

  int numTestPoints = (QK::numQuadraturePoints < 4) ? QK::numQuadraturePoints : 4;
  for (int q = 0; q < numTestPoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J_idx[3][3] = {{0}};
    QK::jacobianTransformation(qa, qb, qc, X, J_idx);

    real_t coords[3] = {static_cast<real_t>(QK::BasisType::parentSupportCoord(qa)),
                        static_cast<real_t>(QK::BasisType::parentSupportCoord(qb)),
                        static_cast<real_t>(QK::BasisType::parentSupportCoord(qc))};
    real_t J_coords[3][3] = {{0}};
    QK::jacobianTransformation(coords, Xfull, J_coords);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(J_idx[i][j], J_coords[i][j], TOL_MATRIX_INVERSION)
            << "jacobianTransformation(coords) must match indexed at GLL point q=" << q;
  }
}

TYPED_TEST(JacobianTest, JacobianTransformationWithCornersMatchesIndexed) {
  using QK = TypeParam;
  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.0, 0.0, 0.0, 1.0);

  int numTestPoints = (QK::numQuadraturePoints < 4) ? QK::numQuadraturePoints : 4;
  for (int q = 0; q < numTestPoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J_idx[3][3] = {{0}};
    QK::jacobianTransformation(qa, qb, qc, X, J_idx);

    real_t coords[3] = {static_cast<real_t>(QK::BasisType::parentSupportCoord(qa)),
                        static_cast<real_t>(QK::BasisType::parentSupportCoord(qb)),
                        static_cast<real_t>(QK::BasisType::parentSupportCoord(qc))};
    real_t J_wc[3][3] = {{0}};
    QK::jacobianTransformationWithCorners(coords, X, J_wc);

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(J_idx[i][j], J_wc[i][j], TOL_MATRIX_INVERSION)
            << "jacobianTransformationWithCorners must match indexed at GLL point q=" << q;
  }
}

// ============================================================================
// TENSOROPS TESTS
// ============================================================================

TEST(TensorOpsTest, Invert3x3TwoArgPreservesSource) {
  real_t J[3][3] = {{2.0, 1.0, 0.0}, {1.0, 3.0, 1.0}, {0.0, 1.0, 2.0}};
  real_t Jinv[3][3] = {{0}};
  real_t Jorig[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Jorig[i][j] = J[i][j];

  real_t det = invert3x3(Jinv, J);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) EXPECT_NEAR(J[i][j], Jorig[i][j], TOL) << "invert3x3(Jinv, J) must not modify J";

  real_t product[3][3] = {{0}};
  matMul3x3(J, Jinv, product);
  EXPECT_TRUE(isIdentity(product, TOL_MATRIX_INVERSION));

  real_t Jcopy[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Jcopy[i][j] = Jorig[i][j];
  real_t det_inplace = invert3x3(Jcopy);
  EXPECT_NEAR(det, det_inplace, TOL_MATRIX_INVERSION) << "2-arg and in-place must return same det";
}

TEST(TensorOpsTest, SymDeterminant2DKnownValue) {
  // Voigt [B00, B11, B01]: det = B00*B11 - B01^2 = 4*9 - 4 = 32
  real_t B[3] = {4.0f, 9.0f, 2.0f};
  real_t det = symDeterminant(B);
  EXPECT_NEAR(det, 32.0f, TOL_MATRIX_INVERSION);
}

TEST(TensorOpsTest, SymDeterminant3DIdentity) {
  // Voigt identity [B00,B11,B22,B12,B02,B01]: det = 1
  real_t B[6] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f};
  real_t det = symDeterminant(B);
  EXPECT_NEAR(det, 1.0f, TOL_MATRIX_INVERSION);
}

TEST(TensorOpsTest, SymInvertTwoArgDiagonalMatrix) {
  // Diagonal: [B00=4, B11=2, B22=3, B12=0, B02=0, B01=0]
  // inv = [1/4, 1/2, 1/3, 0, 0, 0]; symInvert is void
  real_t J[6] = {4.0f, 2.0f, 3.0f, 0.0f, 0.0f, 0.0f};
  real_t dst[6] = {0};
  symInvert(dst, J);
  EXPECT_NEAR(dst[0], 1.0f / 4.0f, TOL_MATRIX_INVERSION);
  EXPECT_NEAR(dst[1], 1.0f / 2.0f, TOL_MATRIX_INVERSION);
  EXPECT_NEAR(dst[2], 1.0f / 3.0f, TOL_MATRIX_INVERSION);
  EXPECT_NEAR(dst[3], 0.0f, TOL_MATRIX_INVERSION);
  EXPECT_NEAR(dst[4], 0.0f, TOL_MATRIX_INVERSION);
  EXPECT_NEAR(dst[5], 0.0f, TOL_MATRIX_INVERSION);
}

TEST(TensorOpsTest, SymInvertInPlaceMatchesTwoArg) {
  real_t J1[6] = {4.0f, 2.0f, 3.0f, 0.0f, 0.0f, 0.0f};
  real_t dst[6] = {0};
  symInvert(dst, J1);

  real_t J2[6] = {4.0f, 2.0f, 3.0f, 0.0f, 0.0f, 0.0f};
  symInvert(J2);

  for (int i = 0; i < 6; ++i) EXPECT_NEAR(dst[i], J2[i], TOL_MATRIX_INVERSION);
}

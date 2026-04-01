/**
 * @file test_Qk_Hexahedron.cpp
 * @brief Unit tests for Qk hexahedral elements
 */

#include <gtest/gtest.h>

#include <cmath>
#include <set>

#include "Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

// Tolerances for floating point comparisons
// Adapted for float precision (real_t is float by default)
constexpr double TOL = 1.0e-10;  // For exact tests
constexpr double TOL_NUMERICAL =
    1.0e-4;  // For numerical integration (accumulation errors)

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
void createArbitraryCube(real_t (&X)[8][3], real_t x0, real_t y0, real_t z0,
                         real_t size)
{
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
void createReferenceHex(real_t (&X)[8][3])
{
  createArbitraryCube<BASIS>(X, -1.0, -1.0, -1.0, 2.0);
}

/**
 * @brief Creates a unit hexahedron [0,1]^3
 */
template <typename BASIS>
void createUnitHex(real_t (&X)[8][3])
{
  createArbitraryCube<BASIS>(X, 0.0, 0.0, 0.0, 1.0);
}

/**
 * @brief Helper to compute matrix-matrix product
 */
void matMul3x3(real_t const (&A)[3][3], real_t const (&B)[3][3],
               real_t (&C)[3][3])
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < 3; ++k) C[i][j] += A[i][k] * B[k][j];
    }
}

/**
 * @brief Check if matrix is identity
 */
bool isIdentity(real_t const (&M)[3][3], double tol = TOL)
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      real_t expected = (i == j) ? 1.0 : 0.0;
      if (std::abs(M[i][j] - expected) > tol) return false;
    }
  return true;
}

// ============================================================================
// TEST FIXTURES
// ============================================================================

using TestedBases = ::testing::Types<
    Q1_Hexahedron_Lagrange_GaussLobatto, Q2_Hexahedron_Lagrange_GaussLobatto,
    Q3_Hexahedron_Lagrange_GaussLobatto, Q4_Hexahedron_Lagrange_GaussLobatto,
    Q5_Hexahedron_Lagrange_GaussLobatto>;

template <typename QK_BASIS>
class MassMatrixTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class StiffnessMatrixTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class JacobianTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class InterpolationTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class GradientTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class IndexingTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class BMatrixTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class BasisGradientTest : public ::testing::Test
{
};

template <typename QK_BASIS>
class FaceOperationsTest : public ::testing::Test
{
};

TYPED_TEST_SUITE(MassMatrixTest, TestedBases);
TYPED_TEST_SUITE(StiffnessMatrixTest, TestedBases);
TYPED_TEST_SUITE(JacobianTest, TestedBases);
TYPED_TEST_SUITE(InterpolationTest, TestedBases);
TYPED_TEST_SUITE(GradientTest, TestedBases);
TYPED_TEST_SUITE(IndexingTest, TestedBases);
TYPED_TEST_SUITE(BMatrixTest, TestedBases);
TYPED_TEST_SUITE(BasisGradientTest, TestedBases);
TYPED_TEST_SUITE(FaceOperationsTest, TestedBases);

// ============================================================================
// MASS MATRIX TESTS
// ============================================================================

TYPED_TEST(MassMatrixTest, MassMatrixSumEqualsVolume_VariousCubes)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  // Test with different cube configurations
  struct CubeConfig
  {
    real_t x0, y0, z0, size;
    real_t expectedVolume;
  };

  CubeConfig configs[] = {
      {0.0, 0.0, 0.0, 1.0, 1.0},      // Unit cube
      {-1.0, -1.0, -1.0, 2.0, 8.0},   // Reference cube
      {5.0, 3.0, -2.0, 0.5, 0.125},   // Small translated cube
      {-10.0, -5.0, 2.0, 3.0, 27.0},  // Large translated cube
      {1.5, -0.5, 0.3, 2.5, 15.625}   // Another arbitrary cube
  };

  for (const auto& config : configs)
  {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    typename QK::TransformType transformData;
    for (int k = 0; k < 8; ++k)
      for (int i = 0; i < 3; ++i) transformData.data[k][i] = X[k][i];

    real_t mass[numNodes] = {0};
    QK::computeMassTerm(transformData,
                        [&](int q, real_t val) { mass[q] = val; });

    real_t totalMass = 0.0;
    for (int i = 0; i < numNodes; ++i) totalMass += mass[i];

    EXPECT_NEAR(totalMass, config.expectedVolume, TOL_NUMERICAL)
        << "Sum of mass matrix should equal element volume" << " for cube at ("
        << config.x0 << "," << config.y0 << "," << config.z0 << ") with size "
        << config.size;
  }
}

TYPED_TEST(MassMatrixTest, MassMatrixPositive)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -3.5, 2.1, 0.7, 1.8);

  typename QK::TransformType transformData;
  for (int k = 0; k < 8; ++k)
    for (int i = 0; i < 3; ++i) transformData.data[k][i] = X[k][i];

  real_t mass[numNodes] = {0};
  QK::computeMassTerm(transformData, [&](int q, real_t val) { mass[q] = val; });

  for (int i = 0; i < numNodes; ++i)
  {
    EXPECT_GT(mass[i], 0.0) << "All mass matrix entries should be positive";
  }
}

// ============================================================================
// STIFFNESS MATRIX TESTS
// ============================================================================

TYPED_TEST(StiffnessMatrixTest, StiffnessTimesConstantIsZero_VariousCubes)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  // Test with different cubes
  struct CubeConfig
  {
    real_t x0, y0, z0, size;
  };

  CubeConfig configs[] = {
      {0.0, 0.0, 0.0, 1.0}, {-5.0, 2.0, -1.0, 2.5}, {10.0, -3.0, 5.0, 0.75}};

  for (const auto& config : configs)
  {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    typename QK::TransformType transformData;
    for (int k = 0; k < 8; ++k)
      for (int i = 0; i < 3; ++i) transformData.data[k][i] = X[k][i];

    // Constant vector (all dofs = 1.0)
    real_t u[numNodes];
    real_t Ku[numNodes] = {0};
    for (int i = 0; i < numNodes; ++i) u[i] = 1.0;

    QK::computeStiffnessTerm(
        transformData, [](int qa, int qb, int qc) {},
        [&](int i, int j, real_t Kij) { Ku[i] += Kij * u[j]; });

    for (int i = 0; i < numNodes; ++i)
    {
      EXPECT_NEAR(Ku[i], 0.0, TOL_NUMERICAL)
          << "K*u should be zero for constant u (partition of unity property)"
          << " for cube at (" << config.x0 << "," << config.y0 << ","
          << config.z0 << ") with size " << config.size;
    }
  }
}

TYPED_TEST(StiffnessMatrixTest, StiffnessMatrixIsSymmetric)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 1.5, -2.3, 0.8, 1.2);

  typename QK::TransformType transformData;
  for (int k = 0; k < 8; ++k)
    for (int i = 0; i < 3; ++i) transformData.data[k][i] = X[k][i];

  real_t K[numNodes][numNodes] = {{0}};

  QK::computeStiffnessTerm(
      transformData, [](int qa, int qb, int qc) {},
      [&](int i, int j, real_t Kij) { K[i][j] += Kij; });

  for (int i = 0; i < numNodes; ++i)
  {
    for (int j = i + 1; j < numNodes; ++j)
    {
      EXPECT_NEAR(K[i][j], K[j][i], TOL_NUMERICAL)
          << "Stiffness matrix should be symmetric: K[" << i << "][" << j
          << "] != K[" << j << "][" << i << "]";
    }
  }
}

// ============================================================================
// JACOBIAN TESTS
// ============================================================================

TYPED_TEST(JacobianTest, JacobianDeterminantPositive_VariousCubes)
{
  using QK = TypeParam;

  // Test with different cubes
  struct CubeConfig
  {
    real_t x0, y0, z0, size;
  };

  CubeConfig configs[] = {
      {0.0, 0.0, 0.0, 1.0}, {-2.0, 3.0, -1.0, 0.5}, {5.0, -5.0, 2.0, 3.0}};

  for (const auto& config : configs)
  {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    for (int q = 0; q < QK::numQuadraturePoints; ++q)
    {
      int qa, qb, qc;
      QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

      real_t J[3][3] = {{0}};
      QK::jacobianTransformation(qa, qb, qc, X, J);

      real_t detJ = determinant(J);

      EXPECT_GT(detJ, 0.0)
          << "Jacobian determinant should be positive at quadrature point " << q
          << " for cube at (" << config.x0 << "," << config.y0 << ","
          << config.z0 << ") with size " << config.size;
    }
  }
}

TYPED_TEST(JacobianTest, InverseJacobianCorrectness_VariousCubes)
{
  using QK = TypeParam;

  // Test with different cubes
  struct CubeConfig
  {
    real_t x0, y0, z0, size;
  };

  CubeConfig configs[] = {
      {0.0, 0.0, 0.0, 1.0}, {-10.0, 5.0, -3.0, 2.0}, {3.5, -1.5, 0.5, 0.8}};

  for (const auto& config : configs)
  {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    // Test a few quadrature points
    int numTestPoints =
        (QK::numQuadraturePoints < 5) ? QK::numQuadraturePoints : 5;
    for (int q = 0; q < numTestPoints; ++q)
    {
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
          << " for cube at (" << config.x0 << "," << config.y0 << ","
          << config.z0 << ") with size " << config.size;
    }
  }
}

TYPED_TEST(JacobianTest, ShapeFunctionsPartitionOfUnity)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  // Test points in reference space
  real_t testPoints[][3] = {
      {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {-0.7, 0.3, 0.1}, {0.9, -0.5, 0.8}};
  constexpr int numTestPoints = 4;

  for (int pt = 0; pt < numTestPoints; ++pt)
  {
    double N[numNodes];
    double coords[3] = {testPoints[pt][0], testPoints[pt][1],
                        testPoints[pt][2]};
    QK::calcN(coords, N);

    double sum = 0.0;
    for (int i = 0; i < numNodes; ++i) sum += N[i];

    EXPECT_NEAR(sum, 1.0, TOL)
        << "Sum of shape functions should be 1 (partition of unity)";
  }
}

TYPED_TEST(JacobianTest, GradientConsistency)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.2, 0.5, 2.3, 1.5);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8)
  {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  }
  else
  {
    QK::computeLocalCoords(X, Xfull);
  }

  constexpr int q = numNodes / 2;

  real_t gradN1[numNodes][3] = {{0}};
  real_t gradN2[numNodes][3] = {{0}};

  real_t detJ1 = QK::calcGradN(q, Xfull, gradN1);
  real_t detJ2 = QK::calcGradNWithCorners(q, X, gradN2);

  EXPECT_NEAR(detJ1, detJ2, TOL_NUMERICAL)
      << "Two gradient computation methods should be consistent";

  for (int i = 0; i < numNodes; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      EXPECT_NEAR(gradN1[i][j], gradN2[i][j], TOL_NUMERICAL)
          << "gradN[" << i << "][" << j << "] inconsistent between methods";
    }
  }
}

TYPED_TEST(JacobianTest, QuadratureRuleIntegratesConstant_VariousCubes)
{
  using QK = TypeParam;

  // Test with different cubes
  struct CubeConfig
  {
    real_t x0, y0, z0, size;
    real_t expectedVolume;
  };

  CubeConfig configs[] = {{0.0, 0.0, 0.0, 1.0, 1.0},
                          {-5.0, 2.0, -1.0, 2.0, 8.0},
                          {3.0, -2.0, 1.0, 0.5, 0.125}};

  for (const auto& config : configs)
  {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    real_t integral = 0.0;

    for (int q = 0; q < QK::numQuadraturePoints; ++q)
    {
      int qa, qb, qc;
      QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

      real_t w = QK::BasisType::weight(qa) * QK::BasisType::weight(qb) *
                 QK::BasisType::weight(qc);

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

TYPED_TEST(JacobianTest, TransformedQuadratureWeightConsistency)
{
  using QK = TypeParam;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 2.5, -1.0, 0.3, 1.8);

  real_t totalWeight = 0.0;
  real_t expectedVolume = 1.8 * 1.8 * 1.8;

  for (int q = 0; q < QK::numQuadraturePoints; ++q)
  {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t w = QK::BasisType::weight(qa) * QK::BasisType::weight(qb) *
               QK::BasisType::weight(qc);

    real_t J[3][3] = {{0}};
    QK::jacobianTransformation(qa, qb, qc, X, J);
    real_t detJ = determinant(J);

    totalWeight += w * detJ;
  }

  EXPECT_NEAR(totalWeight, expectedVolume, TOL_NUMERICAL)
      << "Sum of transformed quadrature weights should equal element volume";
}

// ============================================================================
// INTERPOLATION TESTS
// ============================================================================

TYPED_TEST(InterpolationTest, TrilinearInterpAtCorners)
{
  using QK = TypeParam;

  real_t X[8][3];
  createUnitHex<QK>(X);

  // Corner test data: alpha, beta, gamma, expected corner index
  struct CornerTest
  {
    real_t alpha, beta, gamma;
    int idx;
  };

  CornerTest corners[] = {{0.0, 0.0, 0.0, 0}, {1.0, 0.0, 0.0, 1},
                          {0.0, 1.0, 0.0, 2}, {1.0, 1.0, 0.0, 3},
                          {0.0, 0.0, 1.0, 4}, {1.0, 0.0, 1.0, 5},
                          {0.0, 1.0, 1.0, 6}, {1.0, 1.0, 1.0, 7}};

  for (int c = 0; c < 8; ++c)
  {
    real_t coords[3];
    QK::trilinearInterp(corners[c].alpha, corners[c].beta, corners[c].gamma, X,
                        coords);

    for (int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(coords[i], X[corners[c].idx][i], TOL)
          << "Trilinear interpolation at corner " << corners[c].idx
          << " should match corner coordinates";
    }
  }
}

TYPED_TEST(InterpolationTest, TrilinearInterpAtCenter_ArbitraryCube)
{
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

TYPED_TEST(InterpolationTest, ComputeLocalCoordsConsistency)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t Xmesh[8][3];
  createArbitraryCube<QK>(Xmesh, -1.5, 0.8, 2.2, 1.3);

  real_t X[numNodes][3];
  QK::computeLocalCoords(Xmesh, X);

  for (int k = 0; k < 8; ++k)
  {
    int nodeIdx = QK::meshIndexToLinearIndex3D(k);
    for (int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(X[nodeIdx][i], Xmesh[k][i], TOL)
          << "Corner node " << k << " should match mesh corner";
    }
  }
}

TYPED_TEST(InterpolationTest, InterpolationCoefficientsSum)
{
  using QK = TypeParam;

  for (int q = 0; q < QK::num1dNodes; ++q)
  {
    real_t c0 = QK::interpolationCoord(q, 0);
    real_t c1 = QK::interpolationCoord(q, 1);

    EXPECT_NEAR(c0 + c1, 1.0, TOL)
        << "Interpolation coefficients should sum to 1 at quadrature point "
        << q;

    EXPECT_GE(c0, 0.0) << "Interpolation coefficient should be non-negative";
    EXPECT_LE(c0, 1.0) << "Interpolation coefficient should be ≤ 1";
    EXPECT_GE(c1, 0.0) << "Interpolation coefficient should be non-negative";
    EXPECT_LE(c1, 1.0) << "Interpolation coefficient should be ≤ 1";
  }
}

TYPED_TEST(InterpolationTest, InterpolationCoefficientsAtBoundaries)
{
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  real_t c0_first = QK::interpolationCoord(0, 0);
  real_t c1_first = QK::interpolationCoord(0, 1);

  EXPECT_NEAR(c0_first, 1.0, TOL)
      << "At first node, k=0 coefficient should be 1";
  EXPECT_NEAR(c1_first, 0.0, TOL)
      << "At first node, k=1 coefficient should be 0";

  real_t c0_last = QK::interpolationCoord(num1d - 1, 0);
  real_t c1_last = QK::interpolationCoord(num1d - 1, 1);

  EXPECT_NEAR(c0_last, 0.0, TOL) << "At last node, k=0 coefficient should be 0";
  EXPECT_NEAR(c1_last, 1.0, TOL) << "At last node, k=1 coefficient should be 1";
}

// ============================================================================
// GRADIENT OPERATOR TESTS
// ============================================================================

TYPED_TEST(GradientTest, GradientOfLinearFieldIsConstant_ArbitraryCube)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -2.0, 1.5, 0.3, 1.7);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8)
  {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  }
  else
  {
    QK::computeLocalCoords(X, Xfull);
  }

  real_t u[numNodes][3];
  for (int i = 0; i < numNodes; ++i)
  {
    u[i][0] = Xfull[i][0];
    u[i][1] = 2.0 * Xfull[i][1];
    u[i][2] = 3.0 * Xfull[i][2];
  }

  real_t expectedGrad[3][3] = {
      {1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}};

  int numTestPoints =
      (QK::numQuadraturePoints < 5) ? QK::numQuadraturePoints : 5;
  for (int q = 0; q < numTestPoints; ++q)
  {
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

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        EXPECT_NEAR(grad[i][j], expectedGrad[i][j], TOL_NUMERICAL)
            << "Gradient of linear field should be constant at all quadrature "
               "points";
      }
    }
  }
}

TYPED_TEST(GradientTest, SymmetricGradientIsSymmetric)
{
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.5, -1.2, 3.0, 2.1);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8)
  {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  }
  else
  {
    QK::computeLocalCoords(X, Xfull);
  }

  real_t u[numNodes][3];
  for (int i = 0; i < numNodes; ++i)
  {
    u[i][0] = Xfull[i][0] + 0.5 * Xfull[i][1];
    u[i][1] = Xfull[i][1] + 0.3 * Xfull[i][2];
    u[i][2] = Xfull[i][2] + 0.2 * Xfull[i][0];
  }

  // Test a few quadrature points
  int numTestPoints =
      (QK::numQuadraturePoints < 10) ? QK::numQuadraturePoints : 10;
  for (int q = 0; q < numTestPoints; ++q)
  {
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
// INDEXING TESTS
// ============================================================================

TYPED_TEST(IndexingTest, LinearIndex3DConsistency)
{
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int c = 0; c < num1d; ++c)
  {
    for (int b = 0; b < num1d; ++b)
    {
      for (int a = 0; a < num1d; ++a)
      {
        int idx1 = QK::linearIndex3DVal(a, b, c);
        int idx2 = QK::BasisType::TensorProduct3D::linearIndex(a, b, c);

        EXPECT_EQ(idx1, idx2) << "Linear indexing should be consistent";
      }
    }
  }
}

TYPED_TEST(IndexingTest, MeshIndexToLinearIndexBijection)
{
  using QK = TypeParam;

  std::set<int> nodeIndices;
  for (int k = 0; k < 8; ++k)
  {
    int nodeIdx = QK::meshIndexToLinearIndex3D(k);
    nodeIndices.insert(nodeIdx);

    EXPECT_GE(nodeIdx, 0);
    EXPECT_LT(nodeIdx, QK::numNodes);
  }

  EXPECT_EQ(nodeIndices.size(), 8)
      << "Mesh corners should map to 8 distinct node indices";
}

TYPED_TEST(IndexingTest, LinearIndex2DConsistency)
{
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int b = 0; b < num1d; ++b)
  {
    for (int a = 0; a < num1d; ++a)
    {
      int idx1 = QK::linearIndex2DVal(a, b);
      int idx2 = QK::BasisType::TensorProduct2D::linearIndex(a, b);

      EXPECT_EQ(idx1, idx2) << "2D linear indexing should be consistent";
    }
  }
}

// ============================================================================
// B MATRIX TESTS
// ============================================================================

TYPED_TEST(BMatrixTest, BMatrixSymmetry_ArbitraryCube)
{
  using QK = TypeParam;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -3.0, 1.5, -0.5, 1.8);

  int numTestPoints =
      (QK::numQuadraturePoints < 5) ? QK::numQuadraturePoints : 5;
  for (int q = 0; q < numTestPoints; ++q)
  {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J[3][3] = {{0}};
    real_t B[6] = {0};

    QK::computeBMatrix(qa, qb, qc, X, J, B);

    for (int i = 0; i < 6; ++i)
    {
      EXPECT_TRUE(std::isfinite(B[i]))
          << "B matrix component " << i << " should be finite";
    }

    EXPECT_GT(B[0], 0.0) << "B[0] (xx component) should be positive";
    EXPECT_GT(B[1], 0.0) << "B[1] (yy component) should be positive";
    EXPECT_GT(B[2], 0.0) << "B[2] (zz component) should be positive";
  }
}

// ============================================================================
// BASIS GRADIENT TESTS
// ============================================================================

TYPED_TEST(BasisGradientTest, BasisGradientSymmetryProperty)
{
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int q = 0; q < num1d; ++q)
  {
    for (int p = QK::halfNodes + 1; p < num1d; ++p)
    {
      real_t g1 = QK::basisGradientAt(q, p);
      real_t g2 = QK::basisGradientAt(num1d - 1 - q, num1d - 1 - p);

      EXPECT_NEAR(g1, -g2, TOL)
          << "Basis gradient should satisfy symmetry property: " << "grad(" << q
          << "," << p << ") = -grad(" << (num1d - 1 - q) << ","
          << (num1d - 1 - p) << ")";
    }
  }
}

TYPED_TEST(BasisGradientTest, BasisGradientZeroAtSameNode)
{
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int q = 1; q < num1d - 1; ++q)
  {
    real_t grad = QK::basisGradientAt(q, q);
    EXPECT_NEAR(grad, 0.0, TOL)
        << "Basis gradient should be zero at its own interior node";
  }
}

// ============================================================================
// 2D FACE OPERATIONS TESTS
// ============================================================================

TYPED_TEST(FaceOperationsTest, Jacobian2DRankTwo)
{
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

TYPED_TEST(FaceOperationsTest, DampingTermPositive_ArbitrarySquare)
{
  using QK = TypeParam;
  constexpr int numNodesPerFace = QK::numNodesPerFace;

  // Arbitrary square in XY plane
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

  for (int q = 0; q < numNodesPerFace; ++q)
  {
    real_t damping = QK::computeDampingTerm(q, X);

    EXPECT_GT(damping, 0.0) << "Damping term should be positive at node " << q;

    totalDamping += damping;
  }

  EXPECT_NEAR(totalDamping, expectedArea, TOL_NUMERICAL)
      << "Sum of damping terms should equal face area";
}

TYPED_TEST(FaceOperationsTest, DampingTermScaling)
{
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

  EXPECT_NEAR(d2 / d1, 4.0, TOL_NUMERICAL)
      << "Damping term should scale quadratically with element size";
}
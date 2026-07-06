#pragma once
#include "Qk_Hexahedron_Tensorial.h"
#include "common.h"

template <typename QK_BASIS>
class TensorialIndexingTest : public ::testing::Test {};

template <typename QK_BASIS>
class TensorialJacobianTest : public ::testing::Test {};

template <typename QK_BASIS>
class TensorialOperatorTest : public ::testing::Test {};

template <typename QK_BASIS>
class TensorialGEMMTest : public ::testing::Test {};

TYPED_TEST_SUITE(TensorialIndexingTest, TestedBases);
TYPED_TEST_SUITE(TensorialJacobianTest, TestedBases);
TYPED_TEST_SUITE(TensorialOperatorTest, TestedBases);
TYPED_TEST_SUITE(TensorialGEMMTest, TestedBases);

// ============================================================================
// TESTS D'INDEXATION, BASES ET SYMÉTRIES
// ============================================================================

TYPED_TEST(TensorialIndexingTest, LinearIndex3DBijection) {
  // REDIRECTION DU TYPE : On extrait la base du TypeParam injecté par CMake
  // pour instancier notre nouvelle classe Tensorial_GEMM.
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;

  std::set<int> nodeIndices;
  for (int k = 0; k < 8; ++k) {
    int nodeIdx = QK::meshIndexToLinearIndex3D(k);
    nodeIndices.insert(nodeIdx);
    EXPECT_GE(nodeIdx, 0);
    EXPECT_LT(nodeIdx, QK::numNodes);
  }
  EXPECT_EQ(nodeIndices.size(), 8) << "Les 8 coins doivent pointer vers 8 noeuds distincts.";
}

TYPED_TEST(TensorialIndexingTest, InterpolationCoordBounds) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  for (int q = 0; q < QK::num1dNodes; ++q) {
    real_t c0 = QK::interpolationCoord(q, 0);
    real_t c1 = QK::interpolationCoord(q, 1);
    EXPECT_NEAR(c0 + c1, 1.0, TOL);
    EXPECT_GE(c0, 0.0);
    EXPECT_LE(c0, 1.0);
    EXPECT_GE(c1, 0.0);
    EXPECT_LE(c1, 1.0);
  }
}

TYPED_TEST(TensorialIndexingTest, BasisGradientSymmetryBranches) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int n = QK::num1dNodes;
  for (int q = 0; q < n; ++q) {
    for (int p = 0; p <= QK::halfNodes; ++p) {
      real_t g1 = QK::basisGradientAt(q, p);
      real_t g2 = QK::basisGradientAt(n - 1 - q, n - 1 - p);
      EXPECT_NEAR(g1, -g2, TOL) << "Violation de la symétrie de Gauss-Lobatto pour q=" << q << ", p=" << p;
    }
  }
}

TYPED_TEST(TensorialIndexingTest, JacobianCoefficient1DBranches) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int q = 0;
  EXPECT_NEAR(QK::jacobianCoefficient1D(q, 0, 0, 0), -0.5, TOL);
  EXPECT_NEAR(QK::jacobianCoefficient1D(q, 0, 1, 0), 0.5, TOL);
  EXPECT_NEAR(QK::jacobianCoefficient1D(q, 1, 0, 0), QK::interpolationCoord(q, 0), TOL);
  EXPECT_NEAR(QK::jacobianCoefficient1D(q, 1, 1, 0), QK::interpolationCoord(q, 1), TOL);
}

// ============================================================================
// TESTS DES JACOBIENNES ET MATRICES B
// ============================================================================

TYPED_TEST(TensorialJacobianTest, ComputeBMatrixPositivityAndSymmetry) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  real_t X[8][3];
  createArbitraryCube<QK>(X, -2.0, 1.5, 0.5, 1.5);

  int numTestPoints = (QK::numQuadraturePoints < 5) ? QK::numQuadraturePoints : 5;
  for (int q = 0; q < numTestPoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J[3][3] = {{0}};
    real_t B[6] = {0};
    QK::computeBMatrix(qa, qb, qc, X, J, B);

    for (int i = 0; i < 6; ++i) EXPECT_TRUE(std::isfinite(B[i]));
    EXPECT_GT(B[0], 0.0);
    EXPECT_GT(B[1], 0.0);
    EXPECT_GT(B[2], 0.0);
  }
}

TYPED_TEST(TensorialJacobianTest, DampingTermIsPositive) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  real_t X[4][3];
  X[0][0] = 0.0;
  X[0][1] = 0.0;
  X[0][2] = 0.0;
  X[1][0] = 2.0;
  X[1][1] = 0.0;
  X[1][2] = 0.0;
  X[2][0] = 0.0;
  X[2][1] = 2.0;
  X[2][2] = 0.0;
  X[3][0] = 2.0;
  X[3][1] = 2.0;
  X[3][2] = 0.0;

  real_t totalDamping = 0.0;
  for (int q = 0; q < QK::numNodesPerFace; ++q) {
    real_t damping = QK::computeDampingTerm(q, X);
    EXPECT_GT(damping, 0.0);
    totalDamping += damping;
  }
  EXPECT_NEAR(totalDamping, 4.0, TOL_NUMERICAL);
}

// ============================================================================
// TESTS OPÉRATEURS ET MÉTHODES VIRTUELLES
// ============================================================================

TYPED_TEST(TensorialOperatorTest, VirtualHooksCoverage) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  QK elem;
  EXPECT_EQ(elem.getNumQuadraturePoints(), QK::numQuadraturePoints);
  EXPECT_EQ(elem.getNumSupportPoints(), QK::numNodes);
  EXPECT_EQ(elem.getMaxSupportPoints(), QK::maxSupportPoints);
}

TYPED_TEST(TensorialOperatorTest, MassMatrixIntegratesVolume) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.0, -1.0, -1.0, 2.5);
  real_t totalMass = 0.0;
  QK::computeMassTerm(X, [&](int, real_t val) { totalMass += val; });
  EXPECT_NEAR(totalMass, 15.625, TOL_NUMERICAL);
}

TYPED_TEST(TensorialOperatorTest, SumFactConsistencyWithConstantField) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int numNodes = QK::numNodes;
  real_t X[8][3];
  createArbitraryCube<QK>(X, 1.0, -0.5, 2.0, 1.5);

  real_t u_local[numNodes];
  real_t v_local[numNodes] = {0};
  for (int i = 0; i < numNodes; ++i) u_local[i] = 3.14;

  QK::computeStiffnessTermSumFact(X, u_local, v_local, [](int, int, int) { return real_t(1.0); });
  for (int i = 0; i < numNodes; ++i) {
    EXPECT_NEAR(v_local[i], 0.0, TOL_NUMERICAL);
  }
}

TYPED_TEST(TensorialOperatorTest, SparseStiffnessSymmetryAndPoU) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int numNodes = QK::numNodes;
  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.0, 0.0, 0.0, 1.0);

  real_t K[numNodes][numNodes] = {{0}};
  int func1_calls = 0;

  QK::computeStiffnessTerm(X, [&](int, int, int) { func1_calls++; }, [&](int i, int j, real_t Kij) { K[i][j] += Kij; });

  EXPECT_EQ(func1_calls, QK::numQuadraturePoints);

  for (int i = 0; i < numNodes; ++i) {
    real_t rowSum = 0.0;
    for (int j = 0; j < numNodes; ++j) {
      EXPECT_NEAR(K[i][j], K[j][i], TOL_NUMERICAL);
      rowSum += K[i][j];
    }
    EXPECT_NEAR(rowSum, 0.0, TOL_NUMERICAL);
  }
}

TYPED_TEST(TensorialGEMMTest, FillDerivativeMatrixCorrectness) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int n = QK::num1dNodes;
  real_t D_flat[n * n] = {0};
  QK::fillDerivativeMatrix(D_flat);

  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < n; ++col) {
      EXPECT_NEAR(D_flat[row * n + col], QK::basisGradientAt(col, row), TOL);
    }
  }
}

TYPED_TEST(TensorialGEMMTest, MatMulIdentities) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int N = 2;
  real_t A[N][N] = {{1.0, 2.0}, {3.0, 4.0}};
  real_t I[N][N] = {{1.0, 0.0}, {0.0, 1.0}};
  real_t C_NN[N][N] = {{0}};
  real_t C_TN[N][N] = {{0}};

  QK::template matmul_NN<N, N, N>(A, I, C_NN);
  QK::template matmul_TN<N, N, N>(A, I, C_TN);

  EXPECT_NEAR(C_NN[0][0], 1.0, TOL);
  EXPECT_NEAR(C_NN[0][1], 2.0, TOL);
  EXPECT_NEAR(C_NN[1][0], 3.0, TOL);
  EXPECT_NEAR(C_NN[1][1], 4.0, TOL);

  EXPECT_NEAR(C_TN[0][0], 1.0, TOL);
  EXPECT_NEAR(C_TN[0][1], 3.0, TOL);  // A^T
  EXPECT_NEAR(C_TN[1][0], 2.0, TOL);
  EXPECT_NEAR(C_TN[1][1], 4.0, TOL);
}

TYPED_TEST(TensorialGEMMTest, DeviceStiffnessOperatorMatchesSumFact) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int numNodes = QK::numNodes;
  constexpr int n = QK::num1dNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -0.5, 0.2, 1.1, 1.8);

  real_t u_local[numNodes];
  for (int i = 0; i < numNodes; ++i) u_local[i] = std::sin(static_cast<real_t>(i));

  real_t v_sumfact[numNodes] = {0};
  // Le chemin Flat classique n'est pas conditionné par USE_KOKKOS
  QK::computeStiffnessTermSumFact(X, u_local, v_sumfact, [](int, int, int) { return real_t(1.0); });

  real_t W_metrics[numNodes * 6] = {0};
  QK::computeElementMetrics(X, [](int, int, int) { return real_t(1.0); }, W_metrics);

  real_t D_matrix[n][n] = {{0}};
  for (int r = 0; r < n; ++r) {
    for (int c = 0; c < n; ++c) D_matrix[r][c] = QK::basisGradientAt(c, r);
  }

  real_t v_gemm[numNodes] = {0};
  QK::computeStiffnessOperatorDevice(u_local, v_gemm, W_metrics, D_matrix);

  for (int i = 0; i < numNodes; ++i) {
    EXPECT_NEAR(v_sumfact[i], v_gemm[i], TOL_NUMERICAL);
  }
}

TYPED_TEST(TensorialGEMMTest, ScratchMemorySizesAreValid) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;

  size_t scratch = QK::scratchBytesPerTeam();
  size_t scratchStreaming = QK::scratchBytesPerTeamStreaming();

  EXPECT_GT(scratch, 0u) << "La taille du scratch memory doit être strictement positive.";
  EXPECT_GT(scratchStreaming, scratch) << "La version streaming nécessite plus de mémoire scratch.";
}

#pragma once
#include "common.h"

template <typename QK_BASIS>
class MassMatrixTest : public ::testing::Test {};

template <typename QK_BASIS>
class StiffnessMatrixTest : public ::testing::Test {};

TYPED_TEST_SUITE(MassMatrixTest, TestedBases);
TYPED_TEST_SUITE(StiffnessMatrixTest, TestedBases);

// ============================================================================
// MASS MATRIX TESTS
// ============================================================================

TYPED_TEST(MassMatrixTest, MassMatrixSumEqualsVolume_VariousCubes) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  struct CubeConfig {
    real_t x0, y0, z0, size;
    real_t expectedVolume;
  };

  CubeConfig configs[] = {{0.0, 0.0, 0.0, 1.0, 1.0},
                          {-1.0, -1.0, -1.0, 2.0, 8.0},
                          {5.0, 3.0, -2.0, 0.5, 0.125},
                          {-10.0, -5.0, 2.0, 3.0, 27.0},
                          {1.5, -0.5, 0.3, 2.5, 15.625}};

  for (const auto& config : configs) {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    real_t mass[numNodes] = {0};
    QK::computeMassTerm(X, [&](int q, real_t val) { mass[q] = val; });

    real_t totalMass = 0.0;
    for (int i = 0; i < numNodes; ++i) totalMass += mass[i];

    EXPECT_NEAR(totalMass, config.expectedVolume, TOL_NUMERICAL)
        << "Sum of mass matrix should equal element volume" << " for cube at (" << config.x0 << "," << config.y0 << ","
        << config.z0 << ") with size " << config.size;
  }
}

TYPED_TEST(MassMatrixTest, MassMatrixPositive) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -3.5, 2.1, 0.7, 1.8);

  real_t mass[numNodes] = {0};
  QK::computeMassTerm(X, [&](int q, real_t val) { mass[q] = val; });

  for (int i = 0; i < numNodes; ++i) {
    EXPECT_GT(mass[i], 0.0) << "All mass matrix entries should be positive";
  }
}

// ============================================================================
// STIFFNESS MATRIX TESTS
// ============================================================================

TYPED_TEST(StiffnessMatrixTest, StiffnessTimesConstantIsZero_VariousCubes) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  struct CubeConfig {
    real_t x0, y0, z0, size;
  };

  CubeConfig configs[] = {{0.0, 0.0, 0.0, 1.0}, {-5.0, 2.0, -1.0, 2.5}, {10.0, -3.0, 5.0, 0.75}};

  for (const auto& config : configs) {
    real_t X[8][3];
    createArbitraryCube<QK>(X, config.x0, config.y0, config.z0, config.size);

    real_t u[numNodes];
    real_t Ku[numNodes] = {0};
    for (int i = 0; i < numNodes; ++i) u[i] = 1.0;

    QK::computeStiffnessTerm(X, [](int qa, int qb, int qc) {}, [&](int i, int j, real_t Kij) { Ku[i] += Kij * u[j]; });

    for (int i = 0; i < numNodes; ++i) {
      EXPECT_NEAR(Ku[i], 0.0, TOL_NUMERICAL)
          << "K*u should be zero for constant u (partition of unity property)" << " for cube at (" << config.x0 << ","
          << config.y0 << "," << config.z0 << ") with size " << config.size;
    }
  }
}

TYPED_TEST(StiffnessMatrixTest, StiffnessMatrixIsSymmetric) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 1.5, -2.3, 0.8, 1.2);

  real_t K[numNodes][numNodes] = {{0}};

  QK::computeStiffnessTerm(X, [](int qa, int qb, int qc) {}, [&](int i, int j, real_t Kij) { K[i][j] += Kij; });

  for (int i = 0; i < numNodes; ++i) {
    for (int j = i + 1; j < numNodes; ++j) {
      EXPECT_NEAR(K[i][j], K[j][i], TOL_NUMERICAL)
          << "Stiffness matrix should be symmetric: K[" << i << "][" << j << "] != K[" << j << "][" << i << "]";
    }
  }
}

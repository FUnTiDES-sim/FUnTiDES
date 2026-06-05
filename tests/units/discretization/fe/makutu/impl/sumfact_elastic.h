#pragma once
#include "common.h"

template <typename QK_BASIS>
class SumFactElasticTest : public ::testing::Test {};

TYPED_TEST_SUITE(SumFactElasticTest, TestedBases);

// ============================================================================
// SUM FACTORIZATION — ELASTIC (computeElasticStiffnessSumFact)
// ============================================================================

TYPED_TEST(SumFactElasticTest, ZeroDisplacementGivesZeroForce) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.0, 0.0, 0.0, 1.0);

  real_t u[3][numNodes] = {{0}};
  real_t f[3][numNodes] = {{0}};

  QK::computeElasticStiffnessSumFact(
      X, u, f, [](int, int, int, real_t const(&)[3][3], real_t const(&grad)[3][3], real_t(&flux)[3][3]) {
        for (int p = 0; p < 3; ++p)
          for (int s = 0; s < 3; ++s) flux[p][s] = grad[p][s];
      });

  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i)
      EXPECT_NEAR(f[c][i], 0.0, TOL_NUMERICAL) << "Zero u must give zero f at comp " << c << " node " << i;
}

TYPED_TEST(SumFactElasticTest, ConstantDisplacementGivesZeroForce) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.0, 0.5, 2.0, 1.5);

  real_t u[3][numNodes];
  real_t f[3][numNodes] = {{0}};
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i) u[c][i] = real_t(1) + static_cast<real_t>(c);

  QK::computeElasticStiffnessSumFact(
      X, u, f, [](int, int, int, real_t const(&)[3][3], real_t const(&grad)[3][3], real_t(&flux)[3][3]) {
        for (int p = 0; p < 3; ++p)
          for (int s = 0; s < 3; ++s) flux[p][s] = grad[p][s];
      });

  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i)
      EXPECT_NEAR(f[c][i], 0.0, TOL_NUMERICAL) << "Constant u must give zero f at comp " << c << " node " << i;
}

TYPED_TEST(SumFactElasticTest, NonTrivialDisplacementOutputIsFinite) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -2.0, 1.0, 0.5, 1.8);

  real_t Xfull[numNodes][3];
  if constexpr (numNodes == 8) {
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 3; ++j) Xfull[i][j] = X[i][j];
  } else {
    QK::computeLocalCoords(X, Xfull);
  }

  real_t u[3][numNodes];
  real_t f[3][numNodes] = {{0}};
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i) u[c][i] = Xfull[i][c];

  QK::computeElasticStiffnessSumFact(
      X, u, f, [](int, int, int, real_t const(&)[3][3], real_t const(&grad)[3][3], real_t(&flux)[3][3]) {
        for (int p = 0; p < 3; ++p)
          for (int s = 0; s < 3; ++s) flux[p][s] = grad[p][s];
      });

  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i)
      EXPECT_TRUE(std::isfinite(f[c][i])) << "Force must be finite at comp " << c << " node " << i;
}

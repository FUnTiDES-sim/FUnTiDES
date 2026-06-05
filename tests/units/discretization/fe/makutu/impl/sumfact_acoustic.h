#pragma once
#include "common.h"

template <typename QK_BASIS>
class SumFactAcousticTest : public ::testing::Test {};

TYPED_TEST_SUITE(SumFactAcousticTest, TestedBases);

// ============================================================================
// SUM FACTORIZATION — ACOUSTIC (computeStiffnessTermSumFact)
// ============================================================================

TYPED_TEST(SumFactAcousticTest, ZeroInputGivesZeroOutput) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 1.0, 2.0, 3.0, 1.5);

  real_t p[numNodes] = {0};
  real_t f[numNodes] = {0};
  QK::computeStiffnessTermSumFact(X, p, f, [](int, int, int) { return real_t(1); });

  for (int i = 0; i < numNodes; ++i) EXPECT_NEAR(f[i], 0.0, TOL_NUMERICAL) << "Zero p must give zero f at node " << i;
}

TYPED_TEST(SumFactAcousticTest, ConstantPressureGivesZeroForce) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.0, 2.5, 0.3, 2.0);

  real_t p[numNodes];
  real_t f[numNodes] = {0};
  for (int i = 0; i < numNodes; ++i) p[i] = real_t(3.7);

  QK::computeStiffnessTermSumFact(X, p, f, [](int, int, int) { return real_t(1); });

  for (int i = 0; i < numNodes; ++i) EXPECT_NEAR(f[i], 0.0, TOL_NUMERICAL) << "K*const should be zero at node " << i;
}

TYPED_TEST(SumFactAcousticTest, ConsistencyWithDirectAssembly) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.5, -1.5, 2.0, 1.3);

  real_t p[numNodes];
  for (int i = 0; i < numNodes; ++i) p[i] = std::sin(static_cast<real_t>(i));

  // Direct assembly: K·p via computeStiffnessTerm with alpha=1
  real_t Ku_direct[numNodes] = {0};
  QK::computeStiffnessTerm(X, [](int, int, int) {}, [&](int i, int j, real_t Kij) { Ku_direct[i] += Kij * p[j]; });

  // Sum-factorized: same physics (alpha=1)
  real_t Ku_sumfact[numNodes] = {0};
  QK::computeStiffnessTermSumFact(X, p, Ku_sumfact, [](int, int, int) { return real_t(1); });

  for (int i = 0; i < numNodes; ++i)
    EXPECT_NEAR(Ku_sumfact[i], Ku_direct[i], TOL_NUMERICAL) << "Sum-fact and direct assembly disagree at node " << i;
}

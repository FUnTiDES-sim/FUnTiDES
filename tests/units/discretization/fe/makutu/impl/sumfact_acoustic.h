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

// ============================================================================
// C-PML (computeStiffnessTermSumFactPML)
// ============================================================================

TYPED_TEST(SumFactAcousticTest, PmlZeroProfileMatchesPlainKernel) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.5, -1.5, 2.0, 1.3);

  real_t p[numNodes];
  for (int i = 0; i < numNodes; ++i) p[i] = std::sin(static_cast<real_t>(i));

  // Identity PML profile: kappa=1, coef0=1, coef1=0 (psi stays zero).
  auto identity_pml = [](int, int, int, real_t(&kappa)[3], real_t(&coef0)[3], real_t(&coef1)[3]) {
    for (int j = 0; j < 3; ++j) {
      kappa[j] = real_t(1);
      coef0[j] = real_t(1);
      coef1[j] = real_t(0);
    }
  };

  real_t f_plain[numNodes] = {0};
  QK::computeStiffnessTermSumFact(X, p, f_plain, [](int, int, int) { return real_t(1); });

  real_t f_pml[numNodes] = {0};
  real_t mem[6][numNodes] = {{0}};
  QK::computeStiffnessTermSumFactPML(X, p, f_pml, mem, [](int, int, int) { return real_t(1); }, identity_pml);

  for (int i = 0; i < numNodes; ++i)
    EXPECT_NEAR(f_pml[i], f_plain[i], TOL_NUMERICAL) << "Zero-profile PML must match plain kernel at node " << i;
  for (int j = 0; j < 6; ++j)
    for (int i = 0; i < numNodes; ++i)
      EXPECT_NEAR(mem[j][i], 0.0, TOL_NUMERICAL)
          << "Identity profile must keep mem zero at comp " << j << " node " << i;
}

TYPED_TEST(SumFactAcousticTest, PmlNonZeroProfileChangesForce) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, 0.5, -1.5, 2.0, 1.3);

  real_t p[numNodes];
  for (int i = 0; i < numNodes; ++i) p[i] = std::sin(static_cast<real_t>(i));

  // Non-trivial profile: stretch x only.
  auto pml = [](int, int, int, real_t(&kappa)[3], real_t(&coef0)[3], real_t(&coef1)[3]) {
    for (int j = 0; j < 3; ++j) {
      kappa[j] = (j == 0) ? real_t(1.5) : real_t(1);
      coef0[j] = (j == 0) ? real_t(0.5) : real_t(1);
      coef1[j] = (j == 0) ? real_t(0.3) : real_t(0);
    }
  };

  real_t f_plain[numNodes] = {0};
  QK::computeStiffnessTermSumFact(X, p, f_plain, [](int, int, int) { return real_t(1); });

  real_t f_pml[numNodes] = {0};
  real_t mem[6][numNodes] = {{0}};
  QK::computeStiffnessTermSumFactPML(X, p, f_pml, mem, [](int, int, int) { return real_t(1); }, pml);

  // The stretched operator must differ from the plain one somewhere.
  bool differs = false;
  for (int i = 0; i < numNodes; ++i)
    if (std::fabs(f_pml[i] - f_plain[i]) > TOL_NUMERICAL) differs = true;
  EXPECT_TRUE(differs) << "Non-zero PML profile must change the force";
  // The memory variables must have been advanced (non-zero somewhere).
  bool mem_moved = false;
  for (int j = 0; j < 6; ++j)
    for (int i = 0; i < numNodes; ++i)
      if (std::fabs(mem[j][i]) > TOL_NUMERICAL) mem_moved = true;
  EXPECT_TRUE(mem_moved) << "Non-zero PML profile must advance the memory variables";
}

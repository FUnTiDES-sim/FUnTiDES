#pragma once
#include "common.h"

template <typename QK_BASIS>
class LagrangeBasisTest : public ::testing::Test {};

TYPED_TEST_SUITE(LagrangeBasisTest, TestedBases);

// ============================================================================
// Legacy 1D basis functions — value0, value1, valueBubble
// ============================================================================

TYPED_TEST(LagrangeBasisTest, Value0AndValue1AreFinite) {
  using BASIS = typename TypeParam::BasisType;
  double testPoints[] = {-1.0, -0.5, 0.0, 0.5, 1.0};
  for (double xi : testPoints) {
    EXPECT_TRUE(std::isfinite(BASIS::value0(xi))) << "value0 finite at xi=" << xi;
    EXPECT_TRUE(std::isfinite(BASIS::value1(xi))) << "value1 finite at xi=" << xi;
  }
}

TYPED_TEST(LagrangeBasisTest, Value0IsOneAtMinusOne) {
  using BASIS = typename TypeParam::BasisType;
  // First GLL node always at -1; last always at +1
  EXPECT_NEAR(BASIS::value0(-1.0), 1.0, TOL);
  EXPECT_NEAR(BASIS::value0(1.0), 0.0, TOL);
}

TYPED_TEST(LagrangeBasisTest, Value1IsZeroAtFirstNode) {
  using BASIS = typename TypeParam::BasisType;
  // Lagrange basis: basis_1 is 0 at node 0 (xi=-1), regardless of order
  EXPECT_NEAR(BASIS::value1(-1.0), 0.0, TOL);
}

// ============================================================================
// Legacy 1D gradient functions — gradient0, gradient1
// ============================================================================

TYPED_TEST(LagrangeBasisTest, Gradient0AndGradient1AreFinite) {
  using BASIS = typename TypeParam::BasisType;
  double testPoints[] = {-1.0, -0.5, 0.0, 0.5, 1.0};
  for (double xi : testPoints) {
    EXPECT_TRUE(std::isfinite(BASIS::gradient0(xi))) << "gradient0 finite at xi=" << xi;
    EXPECT_TRUE(std::isfinite(BASIS::gradient1(xi))) << "gradient1 finite at xi=" << xi;
  }
}

TYPED_TEST(LagrangeBasisTest, Gradient0IsDerivativeOfValue0) {
  using BASIS = typename TypeParam::BasisType;
  // Finite difference check at xi=0
  constexpr double h = 1e-5;
  double fd = (BASIS::value0(h) - BASIS::value0(-h)) / (2.0 * h);
  EXPECT_NEAR(BASIS::gradient0(0.0), fd, 1e-4);
}

TYPED_TEST(LagrangeBasisTest, Gradient1IsDerivativeOfValue1) {
  using BASIS = typename TypeParam::BasisType;
  constexpr double h = 1e-5;
  double fd = (BASIS::value1(h) - BASIS::value1(-h)) / (2.0 * h);
  EXPECT_NEAR(BASIS::gradient1(0.0), fd, 1e-4);
}

// ============================================================================
// TensorProduct2D::value — partition of unity
// ============================================================================

TYPED_TEST(LagrangeBasisTest, TensorProduct2DValuePartitionOfUnity) {
  using BASIS = typename TypeParam::BasisType;
  constexpr int n2D = BASIS::TensorProduct2D::numSupportPoints;

  double testCoords[][2] = {{-0.5, 0.3}, {0.0, 0.0}, {0.7, -0.8}, {-1.0, 1.0}};
  for (auto& c : testCoords) {
    double N[n2D];
    BASIS::TensorProduct2D::value(c, N);
    double sum = 0.0;
    for (int i = 0; i < n2D; ++i) sum += N[i];
    EXPECT_NEAR(sum, 1.0, TOL) << "TensorProduct2D::value: partition of unity at (" << c[0] << "," << c[1] << ")";
  }
}

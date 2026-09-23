#pragma once
#include <type_traits>

#include "common.h"

template <typename QK_BASIS>
class LagrangeBasisTest : public ::testing::Test {};

TYPED_TEST_SUITE(LagrangeBasisTest, TestedBases);

// LagrangeBasis3GL exposes only the unrolled value0..value3 / gradient0..3 and
// has no index-dispatching value(index, xi) / gradient(index, xi). The tests
// below must stay compilable for every order, so they are guarded by these
// detectors instead of being duplicated per basis.
template <typename BASIS, typename = void>
struct HasIndexedValue : std::false_type {};
template <typename BASIS>
struct HasIndexedValue<BASIS, std::void_t<decltype(BASIS::value(0, 0.0))>> : std::true_type {};

template <typename BASIS, typename = void>
struct HasIndexedGradient : std::false_type {};
template <typename BASIS>
struct HasIndexedGradient<BASIS, std::void_t<decltype(BASIS::gradient(0, 0.0))>> : std::true_type {};

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

// ============================================================================
// parentSupportCoord — GLL nodes span [-1, 1] in increasing order
// ============================================================================

TYPED_TEST(LagrangeBasisTest, ParentSupportCoordSpansReferenceInterval) {
  using BASIS = typename TypeParam::BasisType;
  constexpr int n = BASIS::numSupportPoints;

  EXPECT_NEAR(BASIS::parentSupportCoord(0), -1.0, TOL) << "first GLL node must be at -1";
  EXPECT_NEAR(BASIS::parentSupportCoord(n - 1), 1.0, TOL) << "last GLL node must be at +1";

  for (int i = 1; i < n; ++i) {
    double const prev = BASIS::parentSupportCoord(i - 1);
    double const cur = BASIS::parentSupportCoord(i);
    EXPECT_GT(cur, prev) << "GLL nodes must be strictly increasing at index " << i;
    EXPECT_LE(std::fabs(cur), 1.0 + TOL) << "GLL node " << i << " outside [-1,1]";
  }

  // Out-of-range index takes the fallback branch; only finiteness is contractual.
  EXPECT_TRUE(std::isfinite(BASIS::parentSupportCoord(n)));
}

// ============================================================================
// value(index, xi) — cardinality and partition of unity
// ============================================================================

TYPED_TEST(LagrangeBasisTest, IndexedValueIsCardinalAtSupportPoints) {
  using BASIS = typename TypeParam::BasisType;
  if constexpr (HasIndexedValue<BASIS>::value) {
    constexpr int n = BASIS::numSupportPoints;

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        double const expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(BASIS::value(i, BASIS::parentSupportCoord(j)), expected, 1e-9)
            << "value(" << i << ", node " << j << ") breaks the cardinal property";
      }

    double testPoints[] = {-0.83, -0.27, 0.0, 0.41, 0.95};
    for (double xi : testPoints) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) sum += BASIS::value(i, xi);
      EXPECT_NEAR(sum, 1.0, TOL) << "value(index, xi): partition of unity at xi=" << xi;
    }

    EXPECT_TRUE(std::isfinite(BASIS::value(n, 0.3)));
  }
}

// ============================================================================
// gradient(index, xi) — derivative of value(index, xi), sums to zero
// ============================================================================

TYPED_TEST(LagrangeBasisTest, IndexedGradientMatchesValueDerivative) {
  using BASIS = typename TypeParam::BasisType;
  if constexpr (HasIndexedGradient<BASIS>::value && HasIndexedValue<BASIS>::value) {
    constexpr int n = BASIS::numSupportPoints;
    constexpr double h = 1e-5;

    double testPoints[] = {-0.83, -0.27, 0.0, 0.41, 0.95};
    for (double xi : testPoints) {
      // The basis sums to 1 everywhere, so its derivative sums to 0.
      double sum = 0.0;
      for (int i = 0; i < n; ++i) sum += BASIS::gradient(i, xi);
      EXPECT_NEAR(sum, 0.0, TOL) << "gradient(index, xi) must sum to zero at xi=" << xi;

      for (int i = 0; i < n; ++i) {
        double const fd = (BASIS::value(i, xi + h) - BASIS::value(i, xi - h)) / (2.0 * h);
        EXPECT_NEAR(BASIS::gradient(i, xi), fd, 1e-4)
            << "gradient(" << i << ", " << xi << ") disagrees with the derivative of value";
      }
    }

    EXPECT_TRUE(std::isfinite(BASIS::gradient(n, 0.3)));
  }
}

// ============================================================================
// TensorProduct3D — index round-trip and partition of unity
// ============================================================================

TYPED_TEST(LagrangeBasisTest, TensorProduct3DIndexRoundTrip) {
  using BASIS = typename TypeParam::BasisType;
  using TP3D = typename BASIS::TensorProduct3D;
  constexpr int n1D = BASIS::numSupportPoints;

  for (int k = 0; k < n1D; ++k)
    for (int j = 0; j < n1D; ++j)
      for (int i = 0; i < n1D; ++i) {
        int const lin = TP3D::linearIndex(i, j, k);
        ASSERT_GE(lin, 0);
        ASSERT_LT(lin, TP3D::numSupportPoints);

        int i0 = -1, i1 = -1, i2 = -1;
        TP3D::multiIndex(lin, i0, i1, i2);
        EXPECT_EQ(i0, i) << "multiIndex(linearIndex(" << i << "," << j << "," << k << ")) broke i";
        EXPECT_EQ(i1, j) << "multiIndex(linearIndex(" << i << "," << j << "," << k << ")) broke j";
        EXPECT_EQ(i2, k) << "multiIndex(linearIndex(" << i << "," << j << "," << k << ")) broke k";
      }
}

TYPED_TEST(LagrangeBasisTest, TensorProduct3DValuePartitionOfUnity) {
  using BASIS = typename TypeParam::BasisType;
  using TP3D = typename BASIS::TensorProduct3D;
  constexpr int n3D = TP3D::numSupportPoints;

  double testCoords[][3] = {{-0.5, 0.3, 0.9}, {0.0, 0.0, 0.0}, {0.7, -0.8, 0.2}, {-1.0, 1.0, -1.0}};
  for (auto& c : testCoords) {
    double N[n3D];
    TP3D::value(c, N);
    double sum = 0.0;
    for (int i = 0; i < n3D; ++i) sum += N[i];
    // Looser than the 2D case: the 3D sum is the cube of the 1D one, so it
    // accumulates roughly three times the round-off at order 9.
    EXPECT_NEAR(sum, 1.0, 1e-9) << "TensorProduct3D::value: partition of unity at (" << c[0] << "," << c[1] << ","
                                << c[2] << ")";
  }
}

/**
 * @file test_pml_coefficients.cc
 * @brief Unit tests for the C-PML coefficient builder (fillPmlCoefficients).
 *
 * The builder implements the Komatitsch & Martin (2007) stretching profiles
 * and the Wang, Lee & Teixeira (2006) exact-convolution recursion documented
 * in pml_coefficients.h. The oracle is that documentation: hand-computed
 * values for a canonical configuration, plus invariants any correct
 * implementation must satisfy (identity in the interior, coef1 = 1 - coef0
 * for a homogeneous layer, coef0 = exp(-a dt)).
 */

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "pml_coefficients.h"

namespace solver {
namespace test {

namespace {

// Canonical configuration for the hand-computed reference values below:
//   domain 1000 m, PML thickness 100 m, vp 3000 m/s, R = 1e-3, N = 2, dt = 1e-3
//   d_max = -(N+1) vp / (2 L) ln(R) = 310.8489875541961
constexpr float kD = 1000.0f;
constexpr float kL = 100.0f;
constexpr float kVp = 3000.0f;
constexpr float kR = 1e-3f;
constexpr float kDt = 1e-3f;
constexpr float kDMax = 310.8489875541961f;

const float kDomain[3] = {kD, kD, kD};
const float kPml[3] = {kL, kL, kL};

// Fill with the canonical configuration, overriding the profile parameters.
PmlCoefficients fill(float x, float y, float z, float profile = 2.0f, float reflection = kR, float alphaMax = 0.0f,
                     float kappaMax = 1.0f, float vp = kVp, const float pml[3] = kPml) {
  PmlCoefficients pc;
  fillPmlCoefficients(x, y, z, kDomain, pml, kDt, vp, profile, reflection, alphaMax, kappaMax, pc);
  return pc;
}

}  // namespace

// ----------------------------------------------------------------------
// Interior: identity coefficients
// ----------------------------------------------------------------------

TEST(PmlCoefficientsTest, InteriorNodeIsIdentity) {
  PmlCoefficients pc = fill(500.0f, 500.0f, 500.0f);
  EXPECT_FALSE(pc.isPml);
  for (int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(pc.d[i], 0.0f);
    EXPECT_FLOAT_EQ(pc.kappa[i], 1.0f);
    EXPECT_FLOAT_EQ(pc.alpha[i], 0.0f);
    EXPECT_FLOAT_EQ(pc.coef0[i], 1.0f);
    EXPECT_FLOAT_EQ(pc.coef1[i], 0.0f);
    EXPECT_FLOAT_EQ(pc.coef2[i], 0.0f);
  }
}

TEST(PmlCoefficientsTest, ZeroPmlSizeIsIdentityEverywhere) {
  const float none[3] = {0.0f, 0.0f, 0.0f};
  PmlCoefficients pc = fill(0.0f, 0.0f, 0.0f, 2.0f, kR, 0.0f, 1.0f, kVp, none);
  EXPECT_FALSE(pc.isPml);
  for (int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(pc.d[i], 0.0f);
    EXPECT_FLOAT_EQ(pc.coef0[i], 1.0f);
    EXPECT_FLOAT_EQ(pc.coef1[i], 0.0f);
  }
}

// The layer is c < L (strict): a node exactly on the inner edge is interior.
TEST(PmlCoefficientsTest, InnerEdgeIsInterior) {
  PmlCoefficients low = fill(kL, 500.0f, 500.0f);
  EXPECT_FALSE(low.isPml);
  EXPECT_FLOAT_EQ(low.d[0], 0.0f);
  EXPECT_FLOAT_EQ(low.coef0[0], 1.0f);

  PmlCoefficients high = fill(kD - kL, 500.0f, 500.0f);
  EXPECT_FALSE(high.isPml);
  EXPECT_FLOAT_EQ(high.d[0], 0.0f);
  EXPECT_FLOAT_EQ(high.coef0[0], 1.0f);
}

// ----------------------------------------------------------------------
// Hand-computed reference values
// ----------------------------------------------------------------------

TEST(PmlCoefficientsTest, OuterBoundaryMatchesHandComputedReference) {
  // c = 0: delta = L, r = 1, rN = 1 -> d = d_max, kappa = 1, alpha = 0.
  PmlCoefficients pc = fill(0.0f, 500.0f, 500.0f);
  EXPECT_TRUE(pc.isPml);
  EXPECT_NEAR(pc.d[0], kDMax, 1e-3f);
  EXPECT_FLOAT_EQ(pc.kappa[0], 1.0f);
  EXPECT_FLOAT_EQ(pc.alpha[0], 0.0f);
  // a = d/kappa = 310.8489875541961; coef0 = exp(-a dt) = 0.732824533138904
  EXPECT_NEAR(pc.coef0[0], 0.732824533138904f, 1e-6f);
  EXPECT_NEAR(pc.coef1[0], 0.267175466861096f, 1e-6f);
  EXPECT_NEAR(pc.coef2[0], pc.coef1[0], 1e-7f);
  // Directions y and z are interior.
  EXPECT_FLOAT_EQ(pc.d[1], 0.0f);
  EXPECT_FLOAT_EQ(pc.coef0[1], 1.0f);
}

TEST(PmlCoefficientsTest, MidLayerMatchesHandComputedReference) {
  // c = 50: delta = 50, r = 0.5, rN = 0.25 -> d = 0.25 d_max = 77.71224688854903
  PmlCoefficients pc = fill(50.0f, 500.0f, 500.0f);
  EXPECT_TRUE(pc.isPml);
  EXPECT_NEAR(pc.d[0], 77.71224688854903f, 1e-3f);
  EXPECT_NEAR(pc.coef0[0], 0.9252306262333163f, 1e-6f);
  EXPECT_NEAR(pc.coef1[0], 1.0f - 0.9252306262333163f, 1e-6f);
}

TEST(PmlCoefficientsTest, KappaAndAlphaProfiles) {
  // r = 0.5, N = 2: kappa = 1 + (kappaMax-1) r^N = 1 + 2*0.25 = 1.5
  //                 alpha = alphaMax (1-r)^N = 2 * 0.25 = 0.5
  PmlCoefficients pc = fill(50.0f, 500.0f, 500.0f, 2.0f, kR, 2.0f, 3.0f);
  EXPECT_NEAR(pc.kappa[0], 1.5f, 1e-6f);
  EXPECT_NEAR(pc.alpha[0], 0.5f, 1e-6f);
  // a = alpha + d/kappa = 0.5 + 77.71224688854903/1.5 = 52.30816459236602
  EXPECT_NEAR(pc.coef0[0], 0.949036362371693f, 1e-6f);
  EXPECT_NEAR(pc.coef1[0], 0.05047648960824684f, 1e-6f);
}

// The two boundaries are symmetric: delta = L at both outer edges.
TEST(PmlCoefficientsTest, HighBoundaryMirrorsLowBoundary) {
  PmlCoefficients lo = fill(0.0f, 500.0f, 500.0f);
  PmlCoefficients hi = fill(kD, 500.0f, 500.0f);
  EXPECT_TRUE(hi.isPml);
  EXPECT_FLOAT_EQ(hi.d[0], lo.d[0]);
  EXPECT_FLOAT_EQ(hi.kappa[0], lo.kappa[0]);
  EXPECT_FLOAT_EQ(hi.coef0[0], lo.coef0[0]);
  EXPECT_FLOAT_EQ(hi.coef1[0], lo.coef1[0]);
}

// ----------------------------------------------------------------------
// Invariants
// ----------------------------------------------------------------------

// With alphaMax = 0 and kappaMax = 1, a = d/kappa = d, so
// coef1 = (d/kappa)/a (1 - coef0) = 1 - coef0 exactly.
TEST(PmlCoefficientsTest, HomogeneousLayerCoef1IsOneMinusCoef0) {
  for (float c = 0.0f; c < kL; c += 5.0f) {
    PmlCoefficients pc = fill(c, 500.0f, 500.0f);
    ASSERT_TRUE(pc.isPml) << "c=" << c;
    EXPECT_NEAR(pc.coef1[0], 1.0f - pc.coef0[0], 1e-6f) << "c=" << c;
    EXPECT_NEAR(pc.coef2[0], pc.coef1[0], 1e-7f) << "c=" << c;
  }
}

// Independent oracle for the recursion: coef0 = exp(-(alpha + d/kappa) dt),
// evaluated from the profile values the builder itself reports.
TEST(PmlCoefficientsTest, Coef0IsExponentialOfRelaxationRate) {
  for (float c = 0.0f; c < kL; c += 7.0f) {
    PmlCoefficients pc = fill(c, 500.0f, 500.0f, 2.0f, kR, 1.5f, 2.5f);
    ASSERT_TRUE(pc.isPml) << "c=" << c;
    float const a = pc.alpha[0] + pc.d[0] / pc.kappa[0];
    EXPECT_NEAR(pc.coef0[0], std::exp(-a * kDt), 1e-6f) << "c=" << c;
  }
}

// Damping is maximal at the outer boundary (c = 0) and decays toward the
// inner edge (c = L).
TEST(PmlCoefficientsTest, DampingDecreasesAwayFromOuterBoundary) {
  float prev = std::numeric_limits<float>::infinity();
  for (float c = 0.0f; c < kL; c += 10.0f) {
    PmlCoefficients pc = fill(c, 500.0f, 500.0f);
    ASSERT_TRUE(pc.isPml) << "c=" << c;
    EXPECT_LT(pc.d[0], prev) << "c=" << c;
    prev = pc.d[0];
  }
}

// A physical layer must be dissipative: 0 < coef0 <= 1 and coef1 >= 0.
TEST(PmlCoefficientsTest, CoefficientsAreDissipative) {
  for (float c = 0.0f; c < kL; c += 3.0f) {
    PmlCoefficients pc = fill(c, 500.0f, 500.0f, 2.0f, kR, 1.0f, 2.0f);
    ASSERT_TRUE(pc.isPml) << "c=" << c;
    EXPECT_GT(pc.coef0[0], 0.0f) << "c=" << c;
    EXPECT_LE(pc.coef0[0], 1.0f) << "c=" << c;
    EXPECT_GE(pc.coef1[0], 0.0f) << "c=" << c;
  }
}

// ----------------------------------------------------------------------
// Degenerate parameter branches
// ----------------------------------------------------------------------

// vp <= 0 -> d_max = 0; with alphaMax = 0 the rate a is 0 and the recursion
// degenerates to the identity.
TEST(PmlCoefficientsTest, ZeroVelocityGivesZeroDMax) {
  PmlCoefficients pc = fill(0.0f, 500.0f, 500.0f, 2.0f, kR, 0.0f, 1.0f, 0.0f);
  EXPECT_TRUE(pc.isPml);  // still inside the layer geometrically
  EXPECT_FLOAT_EQ(pc.d[0], 0.0f);
  EXPECT_FLOAT_EQ(pc.coef0[0], 1.0f);
  EXPECT_FLOAT_EQ(pc.coef1[0], 0.0f);
  EXPECT_FLOAT_EQ(pc.coef2[0], 0.0f);
}

// reflection <= 0 -> d_max = 0.
TEST(PmlCoefficientsTest, NonPositiveReflectionGivesZeroDMax) {
  PmlCoefficients pc = fill(0.0f, 500.0f, 500.0f, 2.0f, 0.0f);
  EXPECT_TRUE(pc.isPml);
  EXPECT_FLOAT_EQ(pc.d[0], 0.0f);
  EXPECT_FLOAT_EQ(pc.coef0[0], 1.0f);
  EXPECT_FLOAT_EQ(pc.coef1[0], 0.0f);
}

// profile <= 0 falls back to the quadratic exponent N = 2.
TEST(PmlCoefficientsTest, NonPositiveProfileFallsBackToQuadratic) {
  PmlCoefficients zero = fill(50.0f, 500.0f, 500.0f, 0.0f);
  PmlCoefficients neg = fill(50.0f, 500.0f, 500.0f, -1.0f);
  PmlCoefficients two = fill(50.0f, 500.0f, 500.0f, 2.0f);
  EXPECT_FLOAT_EQ(zero.d[0], two.d[0]);
  EXPECT_FLOAT_EQ(neg.d[0], two.d[0]);
  EXPECT_FLOAT_EQ(zero.coef0[0], two.coef0[0]);
  EXPECT_FLOAT_EQ(neg.coef0[0], two.coef0[0]);
}

// ----------------------------------------------------------------------
// Per-direction activation
// ----------------------------------------------------------------------

TEST(PmlCoefficientsTest, IsPmlTrueWhenOnlyOneDirectionIsActive) {
  const float pmlX[3] = {kL, 0.0f, 0.0f};
  PmlCoefficients pc = fill(0.0f, 500.0f, 500.0f, 2.0f, kR, 0.0f, 1.0f, kVp, pmlX);
  EXPECT_TRUE(pc.isPml);
  EXPECT_GT(pc.d[0], 0.0f);
  EXPECT_FLOAT_EQ(pc.d[1], 0.0f);
  EXPECT_FLOAT_EQ(pc.d[2], 0.0f);
  EXPECT_FLOAT_EQ(pc.coef0[1], 1.0f);
  EXPECT_FLOAT_EQ(pc.coef0[2], 1.0f);
}

TEST(PmlCoefficientsTest, CornerNodeActivatesAllThreeDirections) {
  PmlCoefficients pc = fill(0.0f, 0.0f, 0.0f);
  EXPECT_TRUE(pc.isPml);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(pc.d[i], kDMax, 1e-3f);
    EXPECT_NEAR(pc.coef0[i], 0.732824533138904f, 1e-6f);
  }
}

}  // namespace test
}  // namespace solver

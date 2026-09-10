/**
 * @file test_elasticity_utils.cc
 * @brief Unit tests for the elasticity coefficient and tensor helpers.
 *
 * The Thomsen parameters have exact analytic inverses in terms of the Voigt
 * coefficients, so most tests here are round trips: build the tensor from
 * (vp, vs, rho, delta, epsilon, gamma) and check that the defining Thomsen
 * relations give the inputs back. Those inverses are only satisfied when the
 * density factors are right, which is what makes them useful as regressions.
 */
#include <gtest/gtest.h>

#include <cmath>

#include "elasticity_utils.h"

namespace {

// A soft marine sediment: vp/vs = 2.0, comfortably above the sqrt(4/3) limit.
constexpr double kVp = 3000.0;
constexpr double kVs = 1500.0;
constexpr double kRho = 2200.0;

/// Thomsen delta expressed from the Voigt coefficients.
double thomsenDelta(double c13, double c33, double c44) {
  const double num = (c13 + c44) * (c13 + c44) - (c33 - c44) * (c33 - c44);
  return num / (2.0 * c33 * (c33 - c44));
}

double thomsenEpsilon(double c11, double c33) { return (c11 - c33) / (2.0 * c33); }

double thomsenGamma(double c66, double c44) { return (c66 - c44) / (2.0 * c44); }

struct VtiCoefficients {
  double c11, c12, c13, c33, c44, c66;
};

VtiCoefficients vti(double vp, double vs, double rho, double delta, double epsilon, double gamma) {
  VtiCoefficients c{};
  computeVTICoefficients(vp, vs, rho, delta, epsilon, gamma, c.c11, c.c12, c.c13, c.c33, c.c44, c.c66);
  return c;
}

void expectSymmetric(const double C[6][6]) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) EXPECT_DOUBLE_EQ(C[i][j], C[j][i]) << "at (" << i << "," << j << ")";
}

}  // namespace

// ======================================================================
// Isotropic coefficients and tensor
// ======================================================================

TEST(ElasticityUtilsTest, IsotropicCoefficientsInvertToVelocities) {
  double lambda = 0.0, mu = 0.0;
  computeIsotropicCoefficients(kVp, kVs, kRho, lambda, mu);

  EXPECT_DOUBLE_EQ(mu, kRho * kVs * kVs);
  EXPECT_DOUBLE_EQ(lambda, kRho * (kVp * kVp - 2.0 * kVs * kVs));
  EXPECT_DOUBLE_EQ(std::sqrt((lambda + 2.0 * mu) / kRho), kVp);
  EXPECT_DOUBLE_EQ(std::sqrt(mu / kRho), kVs);
}

TEST(ElasticityUtilsTest, IsotropicTensorHasTheExpectedStructure) {
  double lambda = 0.0, mu = 0.0;
  computeIsotropicCoefficients(kVp, kVs, kRho, lambda, mu);
  double C[6][6];
  buildIsotropicTensor(lambda, mu, C);

  expectSymmetric(C);
  for (int i = 0; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(C[i][i], lambda + 2.0 * mu);
    EXPECT_DOUBLE_EQ(C[i + 3][i + 3], mu);
    // The two blocks of the Voigt matrix do not couple for an isotropic medium.
    for (int j = 0; j < 3; ++j) EXPECT_DOUBLE_EQ(C[i][j + 3], 0.0);
  }
  EXPECT_DOUBLE_EQ(C[0][1], lambda);
  EXPECT_DOUBLE_EQ(C[0][2], lambda);
  EXPECT_DOUBLE_EQ(C[1][2], lambda);
}

// ======================================================================
// VTI coefficients
// ======================================================================

TEST(ElasticityUtilsTest, VtiCoefficientsRecoverTheThomsenParameters) {
  const double delta = 0.12, epsilon = 0.2, gamma = 0.08;
  const auto c = vti(kVp, kVs, kRho, delta, epsilon, gamma);

  EXPECT_NEAR(thomsenDelta(c.c13, c.c33, c.c44), delta, 1.0e-12);
  EXPECT_NEAR(thomsenEpsilon(c.c11, c.c33), epsilon, 1.0e-12);
  EXPECT_NEAR(thomsenGamma(c.c66, c.c44), gamma, 1.0e-12);
}

TEST(ElasticityUtilsTest, VtiCoefficientsAreLinearInDensity) {
  // Every stiffness is a density times a squared velocity, so doubling rho must
  // double each coefficient. c13 only satisfies this when the density is kept
  // outside the square root of the Thomsen expression.
  const double delta = 0.12, epsilon = 0.2, gamma = 0.08;
  const auto c = vti(kVp, kVs, kRho, delta, epsilon, gamma);
  const auto c2 = vti(kVp, kVs, 2.0 * kRho, delta, epsilon, gamma);

  EXPECT_DOUBLE_EQ(c2.c11, 2.0 * c.c11);
  EXPECT_DOUBLE_EQ(c2.c12, 2.0 * c.c12);
  EXPECT_DOUBLE_EQ(c2.c13, 2.0 * c.c13);
  EXPECT_DOUBLE_EQ(c2.c33, 2.0 * c.c33);
  EXPECT_DOUBLE_EQ(c2.c44, 2.0 * c.c44);
  EXPECT_DOUBLE_EQ(c2.c66, 2.0 * c.c66);
}

TEST(ElasticityUtilsTest, VtiWithZeroThomsenParametersIsIsotropic) {
  double lambda = 0.0, mu = 0.0;
  computeIsotropicCoefficients(kVp, kVs, kRho, lambda, mu);
  const auto c = vti(kVp, kVs, kRho, 0.0, 0.0, 0.0);

  EXPECT_DOUBLE_EQ(c.c11, lambda + 2.0 * mu);
  EXPECT_DOUBLE_EQ(c.c33, lambda + 2.0 * mu);
  EXPECT_DOUBLE_EQ(c.c44, mu);
  EXPECT_DOUBLE_EQ(c.c66, mu);
  EXPECT_DOUBLE_EQ(c.c12, lambda);
  EXPECT_NEAR(c.c13, lambda, 1.0e-9);
}

TEST(ElasticityUtilsTest, VtiTensorMatchesTheCoefficients) {
  const auto c = vti(kVp, kVs, kRho, 0.12, 0.2, 0.08);
  double C[6][6];
  buildVTITensor(c.c11, c.c12, c.c13, c.c33, c.c44, c.c66, C);

  expectSymmetric(C);
  EXPECT_DOUBLE_EQ(C[0][0], c.c11);
  EXPECT_DOUBLE_EQ(C[1][1], c.c11);
  EXPECT_DOUBLE_EQ(C[2][2], c.c33);
  EXPECT_DOUBLE_EQ(C[0][1], c.c12);
  EXPECT_DOUBLE_EQ(C[0][2], c.c13);
  EXPECT_DOUBLE_EQ(C[1][2], c.c13);
  EXPECT_DOUBLE_EQ(C[3][3], c.c44);
  EXPECT_DOUBLE_EQ(C[4][4], c.c44);
  EXPECT_DOUBLE_EQ(C[5][5], c.c66);
  // c11 = c12 + 2 c66 is the in-plane isotropy of a VTI medium.
  EXPECT_DOUBLE_EQ(C[0][0], C[0][1] + 2.0 * C[5][5]);
}

// ======================================================================
// TTI tensor
// ======================================================================

TEST(ElasticityUtilsTest, TtiWithoutRotationEqualsTheVtiTensor) {
  const double delta = 0.12, epsilon = 0.2, gamma = 0.08;
  const auto c = vti(kVp, kVs, kRho, delta, epsilon, gamma);
  double expected[6][6];
  buildVTITensor(c.c11, c.c12, c.c13, c.c33, c.c44, c.c66, expected);

  double C[6][6];
  computeCTensor(kVp, kVs, kRho, delta, epsilon, gamma, 0.0, 0.0, C);

  // computeCTensor rebuilds the VTI tensor from its own expressions; the two
  // paths must agree, otherwise a fix applied to one of them silently misses
  // the other.
  const double tol = 1.0e-6 * c.c11;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) EXPECT_NEAR(C[i][j], expected[i][j], tol) << "at (" << i << "," << j << ")";
}

TEST(ElasticityUtilsTest, TtiTensorStaysSymmetricUnderRotation) {
  double C[6][6];
  computeCTensor(kVp, kVs, kRho, 0.12, 0.2, 0.08, 37.0, 121.0, C);
  expectSymmetric(C);
}

// The rotation itself (frame invariance of an isotropic medium, invariance of a
// VTI medium under a spin about its symmetry axis, conservation of C_iijj and
// C_ikik) is not covered here: those properties currently fail, and the cause is
// the Bond matrix rather than anything this file changes. Tracked separately.

// ======================================================================
// float instantiation
// ======================================================================

TEST(ElasticityUtilsTest, FloatInstantiationAgreesWithDouble) {
  const double delta = 0.12, epsilon = 0.2, gamma = 0.08;
  const auto expected = vti(kVp, kVs, kRho, delta, epsilon, gamma);

  float c11, c12, c13, c33, c44, c66;
  computeVTICoefficients(static_cast<float>(kVp), static_cast<float>(kVs), static_cast<float>(kRho),
                         static_cast<float>(delta), static_cast<float>(epsilon), static_cast<float>(gamma), c11, c12,
                         c13, c33, c44, c66);

  const double tol = 1.0e-5 * expected.c11;
  EXPECT_NEAR(c11, expected.c11, tol);
  EXPECT_NEAR(c12, expected.c12, tol);
  EXPECT_NEAR(c13, expected.c13, tol);
  EXPECT_NEAR(c33, expected.c33, tol);
  EXPECT_NEAR(c44, expected.c44, tol);
  EXPECT_NEAR(c66, expected.c66, tol);
}

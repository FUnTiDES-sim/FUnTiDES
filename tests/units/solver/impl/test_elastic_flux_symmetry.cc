/**
 * @file test_elastic_flux_symmetry.cc
 * @brief Major-symmetry regression test for the elastic stiffness contraction.
 *
 * The constitutive callbacks of computeElementContributions_{Iso,Vti,Tti} map a
 * reference displacement gradient grad_u_ref (3x3) to a reference-element flux
 * (3x3). That map is linear, i.e. a 9x9 matrix M. The stiffness operator has
 * major symmetry iff M == M^T. A previous implementation rebuilt a rank-4 tensor
 * per (p,r) pair and reused the same off-diagonal coefficients for the (f,l) and
 * (l,f) contractions, which breaks the symmetry as soon as J_inv is not
 * diagonal (the residual is proportional to (lambda - mu) times the
 * antisymmetric part of the gradient). This test locks the invariant so the
 * pattern cannot silently come back.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "elastic_flux.h"

namespace solver {
namespace fe {
namespace test {

namespace {

// Deliberately non-diagonal inverse Jacobian: a diagonal J_inv hides the bug
// (see DiagonalJacobianStaysSymmetric below).
constexpr float kJinv[3][3] = {{1.3f, 0.2f, -0.4f}, {0.1f, 0.9f, 0.3f}, {-0.2f, 0.5f, 1.1f}};

constexpr float kDiagJinv[3][3] = {{1.3f, 0.0f, 0.0f}, {0.0f, 0.9f, 0.0f}, {0.0f, 0.0f, 1.1f}};

/**
 * @brief Assemble the 9x9 matrix of the linear map grad_u_ref -> flux and check
 *   M(i,j) == M(j,i). Column c is obtained by feeding a unit gradient.
 * @tparam FLUX Callable (float const (&grad)[3][3], float (&flux)[3][3]).
 */
template <typename FLUX>
void ExpectMajorSymmetry(FLUX&& flux_fn, float scale) {
  float M[9][9];
  float max_abs = 0.0f;
  for (int c = 0; c < 9; ++c) {
    float grad[3][3] = {};
    grad[c / 3][c % 3] = 1.0f;
    float flux[3][3];
    flux_fn(grad, flux);
    for (int r = 0; r < 9; ++r) {
      M[r][c] = flux[r / 3][r % 3];
      max_abs = std::fmax(max_abs, std::fabs(M[r][c]));
    }
  }

  EXPECT_GT(max_abs, 0.0f) << "flux map is identically zero";

  float const tol = 1e-4f * std::fmax(max_abs, scale);
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 9; ++j)
      EXPECT_NEAR(M[i][j], M[j][i], tol) << "major symmetry broken at (" << i << "," << j << ")";
}

}  // namespace

// ======================================================================
// Isotropic
// ======================================================================
TEST(ElasticFluxSymmetry, Isotropic) {
  // Vp = 3000, Vs = 1500, rho = 2000 -> lambda != mu (lambda == mu only at
  // Vp/Vs = sqrt(3), which would mask the residual).
  float const mu = 2000.0f * 1500.0f * 1500.0f;
  float const lambda = 2000.0f * (3000.0f * 3000.0f - 2.0f * 1500.0f * 1500.0f);
  ExpectMajorSymmetry([&](float const(&g)[3][3], float(&f)[3][3]) { flux::elasticFluxIso(kJinv, mu, lambda, g, f); },
                      mu + lambda);
}

// ======================================================================
// VTI — also guards the second, independent VTI bug (c66<->c12, c44<->c13
// swapped in the old off-diagonal coefficients).
// ======================================================================
TEST(ElasticFluxSymmetry, Vti) {
  float const rho_vp2 = 2000.0f * 3200.0f * 3200.0f;
  float const rho_vs2 = 2000.0f * 1700.0f * 1700.0f;
  float const c33 = rho_vp2;
  float const c44 = rho_vs2;
  float const c11 = rho_vp2 * (1.0f + 2.0f * 0.15f);  // epsilon = 0.15
  float const c66 = rho_vs2 * (1.0f + 2.0f * 0.08f);  // gamma = 0.08
  float const c13 = 0.6f * rho_vp2;                   // c13 != c44
  float const c12 = c11 - 2.0f * c66;                 // c12 != c66
  ExpectMajorSymmetry(
      [&](float const(&g)[3][3], float(&f)[3][3]) { flux::elasticFluxVti(kJinv, c11, c12, c13, c33, c44, c66, g, f); },
      c11);
}

// ======================================================================
// TTI — arbitrary symmetric 6x6 stiffness matrix.
// ======================================================================
TEST(ElasticFluxSymmetry, Tti) {
  float C[6][6];
  float v = 1.0f;
  for (int a = 0; a < 6; ++a)
    for (int b = a; b < 6; ++b) {
      float const val = 1.0e9f * (1.0f + 0.31f * v);
      C[a][b] = val;
      C[b][a] = val;
      v += 1.0f;
    }
  ExpectMajorSymmetry([&](float const(&g)[3][3], float(&f)[3][3]) { flux::elasticFluxTti(kJinv, C, g, f); }, C[0][0]);
}

// ======================================================================
// A non-symmetric CTTI must produce a non-symmetric map: proves the test
// is actually sensitive to the constitutive tensor, not trivially passing.
// ======================================================================
TEST(ElasticFluxSymmetry, TtiAsymmetricTensorIsDetected) {
  float C[6][6] = {};
  for (int a = 0; a < 6; ++a) C[a][a] = 1.0e9f;
  C[0][3] = 5.0e8f;  // no matching C[3][0]
  float M[9][9];
  for (int c = 0; c < 9; ++c) {
    float grad[3][3] = {};
    grad[c / 3][c % 3] = 1.0f;
    float f[3][3];
    flux::elasticFluxTti(kJinv, C, grad, f);
    for (int r = 0; r < 9; ++r) M[r][c] = f[r / 3][r % 3];
  }
  float max_asym = 0.0f;
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 9; ++j) max_asym = std::fmax(max_asym, std::fabs(M[i][j] - M[j][i]));
  EXPECT_GT(max_asym, 1.0f) << "test insensitive to tensor asymmetry";
}

// ======================================================================
// Documents why the homogeneous cartesian validation could never see the
// bug: with a diagonal J_inv the map is symmetric regardless.
// ======================================================================
TEST(ElasticFluxSymmetry, DiagonalJacobianStaysSymmetric) {
  float const mu = 2000.0f * 1500.0f * 1500.0f;
  float const lambda = 2000.0f * (3000.0f * 3000.0f - 2.0f * 1500.0f * 1500.0f);
  ExpectMajorSymmetry(
      [&](float const(&g)[3][3], float(&f)[3][3]) { flux::elasticFluxIso(kDiagJinv, mu, lambda, g, f); }, mu + lambda);
}

}  // namespace test
}  // namespace fe
}  // namespace solver

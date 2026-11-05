#pragma once
// Tests of the 2-D element Qk_Quadrilateral_Lagrange_GaussLobatto.
// TestedBases must be defined by the per-order driver before this header is included.

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Qk_Quadrilateral_Lagrange_GaussLobatto.h"

namespace quad_test {

#ifdef USE_DOUBLE
constexpr double kTol = 1.0e-10;
#else
constexpr double kTol = 2.0e-4;
#endif

/// Quadrilateral from its 4 vertices in counter-clockwise order (v0, v1, v2, v3);
/// the element numbering is k % 2 along x, k / 2 along z: X = {v0, v1, v3, v2}.
inline void makeQuad(real_t (&X)[4][2], real_t const (&v)[4][2]) {
  int const order[4] = {0, 1, 3, 2};
  for (int k = 0; k < 4; ++k) {
    X[k][0] = v[order[k]][0];
    X[k][1] = v[order[k]][1];
  }
}

inline void makeRectangle(real_t (&X)[4][2], real_t x0, real_t z0, real_t lx, real_t lz) {
  real_t const v[4][2] = {{x0, z0}, {x0 + lx, z0}, {x0 + lx, z0 + lz}, {x0, z0 + lz}};
  makeQuad(X, v);
}

/// Shoelace area of the counter-clockwise polygon v.
inline real_t area(real_t const (&v)[4][2]) {
  real_t a = 0;
  for (int k = 0; k < 4; ++k) {
    int const n = (k + 1) % 4;
    a += v[k][0] * v[n][1] - v[n][0] * v[k][1];
  }
  return 0.5 * a;
}

/// A rectangle, a parallelogram and two general (non-affine) convex quadrilaterals.
inline std::vector<std::vector<real_t>> shapes() {
  return {{0, 0, 1, 0, 1, 1, 0, 1},
          {2, -1, 5, -1, 5, 0.5, 2, 0.5},
          {0, 0, 2, 0.5, 2.5, 2, 0.5, 1.5},
          {-1, -1, 1.2, -0.8, 0.9, 1.1, -0.7, 0.8},
          {10, 20, 13, 21, 12.5, 24, 9.5, 23}};
}

inline void shapeVertices(std::vector<real_t> const& s, real_t (&v)[4][2]) {
  for (int k = 0; k < 4; ++k) {
    v[k][0] = s[2 * k];
    v[k][1] = s[2 * k + 1];
  }
}

/// Physical coordinates of the element nodes (bilinear map).
template <typename QK>
void nodeCoords(real_t const (&X)[4][2], real_t (&x)[QK::numNodes], real_t (&z)[QK::numNodes]) {
  for (int qb = 0; qb < QK::num1dNodes; ++qb)
    for (int qa = 0; qa < QK::num1dNodes; ++qa) {
      int const q = QK::linearIndex2D(qa, qb);
      x[q] = 0;
      z[q] = 0;
      for (int k = 0; k < 4; ++k) {
        real_t const n = QK::interpolationCoord(qa, k % 2) * QK::interpolationCoord(qb, k / 2);
        x[q] += n * X[k][0];
        z[q] += n * X[k][1];
      }
    }
}

/// Element stiffness matrix (unit coefficient), column by column.
template <typename QK>
void stiffnessMatrix(real_t const (&X)[4][2], std::vector<real_t>& K) {
  constexpr int n = QK::numNodes;
  K.assign(n * n, 0);
  for (int j = 0; j < n; ++j) {
    real_t u[n] = {0};
    real_t v[n] = {0};
    u[j] = 1;
    QK::computeStiffnessTermSumFact(X, u, v, [](int, int) { return real_t(1); });
    for (int i = 0; i < n; ++i) K[i * n + j] = v[i];
  }
}

}  // namespace quad_test

template <typename QK>
class QuadrilateralTest : public ::testing::Test {};
TYPED_TEST_SUITE(QuadrilateralTest, TestedBases);

// ---------------------------------------------------------------------------
// Numbering
// ---------------------------------------------------------------------------

TYPED_TEST(QuadrilateralTest, VertexAndEdgeNumbering) {
  using QK = TypeParam;
  constexpr int n = QK::num1dNodes;
  EXPECT_EQ(QK::numNodes, n * n);
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(0), 0);
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(1), n - 1);
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(2), (n - 1) * n);
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(3), n * n - 1);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(QK::edgeNodeIndex(0, i), i);
    EXPECT_EQ(QK::edgeNodeIndex(1, i), (n - 1) + i * n);
    EXPECT_EQ(QK::edgeNodeIndex(2, i), i + (n - 1) * n);
    EXPECT_EQ(QK::edgeNodeIndex(3, i), i * n);
  }
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

TYPED_TEST(QuadrilateralTest, RectangleJacobianIsDiagonalHalfSizes) {
  using QK = TypeParam;
  real_t X[4][2];
  quad_test::makeRectangle(X, 3.0, -2.0, 4.0, 1.5);
  for (int qb = 0; qb < QK::num1dNodes; ++qb)
    for (int qa = 0; qa < QK::num1dNodes; ++qa) {
      real_t J[2][2];
      QK::jacobianTransformation(qa, qb, X, J);
      EXPECT_NEAR(J[0][0], 2.0, quad_test::kTol);
      EXPECT_NEAR(J[1][1], 0.75, quad_test::kTol);
      EXPECT_NEAR(J[0][1], 0.0, quad_test::kTol);
      EXPECT_NEAR(J[1][0], 0.0, quad_test::kTol);
    }
}

TYPED_TEST(QuadrilateralTest, BMatrixIsDetTimesInverseMetric) {
  using QK = TypeParam;
  for (auto const& s : quad_test::shapes()) {
    real_t v[4][2], X[4][2];
    quad_test::shapeVertices(s, v);
    quad_test::makeQuad(X, v);
    for (int qb = 0; qb < QK::num1dNodes; ++qb)
      for (int qa = 0; qa < QK::num1dNodes; ++qa) {
        real_t J[2][2], B[3];
        real_t const det = QK::computeBMatrix(qa, qb, X, J, B);
        EXPECT_GT(det, 0.0);
        // B J^T J = |det J| I
        real_t const M[2][2] = {{J[0][0] * J[0][0] + J[1][0] * J[1][0], J[0][0] * J[0][1] + J[1][0] * J[1][1]},
                                {J[0][0] * J[0][1] + J[1][0] * J[1][1], J[0][1] * J[0][1] + J[1][1] * J[1][1]}};
        real_t const Bm[2][2] = {{B[0], B[2]}, {B[2], B[1]}};
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j) {
            real_t p = 0;
            for (int k = 0; k < 2; ++k) p += Bm[i][k] * M[k][j];
            EXPECT_NEAR(p, i == j ? std::abs(det) : 0.0, quad_test::kTol * (1 + std::abs(det)));
          }
      }
  }
}

// ---------------------------------------------------------------------------
// Mass
// ---------------------------------------------------------------------------

TYPED_TEST(QuadrilateralTest, MassSumEqualsAreaAndIsPositive) {
  using QK = TypeParam;
  for (auto const& s : quad_test::shapes()) {
    real_t v[4][2], X[4][2];
    quad_test::shapeVertices(s, v);
    quad_test::makeQuad(X, v);
    real_t mass[QK::numNodes] = {0};
    QK::computeMassTerm(X, [&](int q, real_t val) { mass[q] = val; });
    real_t total = 0;
    for (int q = 0; q < QK::numNodes; ++q) {
      EXPECT_GT(mass[q], 0.0);
      total += mass[q];
    }
    real_t const a = quad_test::area(v);
    EXPECT_NEAR(total, a, quad_test::kTol * a);
  }
}

// ---------------------------------------------------------------------------
// Stiffness
// ---------------------------------------------------------------------------

TYPED_TEST(QuadrilateralTest, StiffnessIsSymmetricWithConstantsInKernel) {
  using QK = TypeParam;
  constexpr int n = QK::numNodes;
  for (auto const& s : quad_test::shapes()) {
    real_t v[4][2], X[4][2];
    quad_test::shapeVertices(s, v);
    quad_test::makeQuad(X, v);
    std::vector<real_t> K;
    quad_test::stiffnessMatrix<QK>(X, K);
    real_t kmax = 0;
    for (real_t k : K) kmax = std::max(kmax, std::abs(k));
    for (int i = 0; i < n; ++i) {
      real_t rowSum = 0;
      for (int j = 0; j < n; ++j) {
        EXPECT_NEAR(K[i * n + j], K[j * n + i], quad_test::kTol * kmax);
        rowSum += K[i * n + j];
      }
      EXPECT_NEAR(rowSum, 0.0, quad_test::kTol * kmax);
      EXPECT_GT(K[i * n + i], 0.0);
    }
  }
}

TYPED_TEST(QuadrilateralTest, StiffnessEnergyOfLinearFieldIsExact) {
  // u = a x + b z is in the element space for any bilinear map, and
  // u^T K u = |grad u|^2 * area exactly with GLL quadrature.
  using QK = TypeParam;
  constexpr int n = QK::numNodes;
  real_t const a = 0.7, b = -1.3;
  for (auto const& s : quad_test::shapes()) {
    real_t v[4][2], X[4][2];
    quad_test::shapeVertices(s, v);
    quad_test::makeQuad(X, v);
    real_t x[n], z[n], u[n], Ku[n] = {0};
    quad_test::nodeCoords<QK>(X, x, z);
    for (int q = 0; q < n; ++q) u[q] = a * x[q] + b * z[q];
    QK::computeStiffnessTermSumFact(X, u, Ku, [](int, int) { return real_t(1); });
    real_t energy = 0;
    for (int q = 0; q < n; ++q) energy += u[q] * Ku[q];
    real_t const expected = (a * a + b * b) * quad_test::area(v);
    EXPECT_NEAR(energy, expected, 10 * quad_test::kTol * (1 + expected));
  }
}

TYPED_TEST(QuadrilateralTest, StiffnessScalesWithCoefficient) {
  using QK = TypeParam;
  constexpr int n = QK::numNodes;
  real_t v[4][2], X[4][2];
  quad_test::shapeVertices(quad_test::shapes()[3], v);
  quad_test::makeQuad(X, v);
  real_t u[n], K1[n] = {0}, K3[n] = {0};
  for (int q = 0; q < n; ++q) u[q] = std::sin(0.7 * q) + 0.1 * q;
  QK::computeStiffnessTermSumFact(X, u, K1, [](int, int) { return real_t(1); });
  QK::computeStiffnessTermSumFact(X, u, K3, [](int, int) { return real_t(3); });
  for (int q = 0; q < n; ++q) EXPECT_NEAR(K3[q], 3 * K1[q], quad_test::kTol * (1 + std::abs(K3[q])));
}

// ---------------------------------------------------------------------------
// Edges and shape functions
// ---------------------------------------------------------------------------

TYPED_TEST(QuadrilateralTest, DampingTermSumsToEdgeLength) {
  using QK = TypeParam;
  real_t const E[2][2] = {{1.0, 2.0}, {4.0, 6.0}};  // length 5
  real_t total = 0;
  for (int i = 0; i < QK::numNodesPerEdge; ++i) {
    real_t const d = QK::computeDampingTerm(i, E);
    EXPECT_GT(d, 0.0);
    total += d;
  }
  EXPECT_NEAR(total, 5.0, quad_test::kTol * 5);
}

TYPED_TEST(QuadrilateralTest, ShapeFunctionsInterpolate) {
  using QK = TypeParam;
  // Kronecker property at the nodes (node coordinates are real_t: float rounding)
  constexpr double kNodeTol = 1e-6;
  for (int qb = 0; qb < QK::num1dNodes; ++qb)
    for (int qa = 0; qa < QK::num1dNodes; ++qa) {
      double const c[2] = {QK::parentCoord(qa), QK::parentCoord(qb)};
      double N[QK::numNodes];
      QK::calcN(c, N);
      for (int q = 0; q < QK::numNodes; ++q) EXPECT_NEAR(N[q], q == QK::linearIndex2D(qa, qb) ? 1.0 : 0.0, kNodeTol);
    }
  // Partition of unity and exact bilinear interpolation inside the element
  double const pts[3][2] = {{0.1, -0.4}, {-0.77, 0.9}, {0.33, 0.25}};
  for (auto const& p : pts) {
    double const c[2] = {p[0], p[1]};
    double N[QK::numNodes];
    QK::calcN(c, N);
    double sum = 0, f = 0;
    for (int qb = 0; qb < QK::num1dNodes; ++qb)
      for (int qa = 0; qa < QK::num1dNodes; ++qa) {
        int const q = QK::linearIndex2D(qa, qb);
        sum += N[q];
        f += N[q] * (1 + 2 * QK::parentCoord(qa) - QK::parentCoord(qb) + 0.5 * QK::parentCoord(qa) * QK::parentCoord(qb));
      }
    EXPECT_NEAR(sum, 1.0, 1e-10);
    EXPECT_NEAR(f, 1 + 2 * p[0] - p[1] + 0.5 * p[0] * p[1], kNodeTol);
  }
}

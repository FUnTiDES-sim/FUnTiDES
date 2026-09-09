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

// ============================================================================
// SUM FACTORIZATION — ELASTIC, TEAM VARIANTS
// (computeElasticStiffnessSumFactTeam)
// ============================================================================

/**
 * @brief Pass-through constitutive callback, callable from a device kernel.
 */
struct IdentityFlux {
  KOKKOS_INLINE_FUNCTION void operator()(int, int, int, real_t const (&)[3][3], real_t const (&grad)[3][3],
                                         real_t (&flux)[3][3]) const {
    for (int p = 0; p < 3; ++p)
      for (int s = 0; s < 3; ++s) flux[p][s] = grad[p][s];
  }
};

/**
 * @brief Runs one team over a single element, general (per-point Jacobian)
 *   overload.
 */
template <typename QK>
inline void launchElasticTeam(Kokkos::View<float[8][3]> X_dev, Kokkos::View<real_t*> u_dev, Kokkos::View<real_t*> f_dev,
                              Kokkos::View<real_t*> F_dev) {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamMember = TeamPolicy::member_type;

  Kokkos::parallel_for(
      TeamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamMember& team) {
        float X_local[8][3];
        for (int i = 0; i < 8; ++i)
          for (int j = 0; j < 3; ++j) X_local[i][j] = X_dev(i, j);

        QK::computeElasticStiffnessSumFactTeam(team, X_local, u_dev.data(), f_dev.data(), F_dev.data(), IdentityFlux{});
      });
}

/**
 * @brief Same, constant-Jacobian overload: the caller evaluates the geometry
 *   once and passes it in.
 */
template <typename QK>
inline void launchElasticTeamConstJac(Kokkos::View<real_t*> geom_dev, Kokkos::View<real_t*> u_dev,
                                      Kokkos::View<real_t*> f_dev, Kokkos::View<real_t*> F_dev) {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamMember = TeamPolicy::member_type;

  Kokkos::parallel_for(
      TeamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamMember& team) {
        QK::computeElasticStiffnessSumFactTeam(team, geom_dev.data(), u_dev.data(), f_dev.data(), F_dev.data(),
                                               IdentityFlux{});
      });
}

TYPED_TEST(SumFactElasticTest, TeamMatchesSerial) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.0, 0.5, 0.0, 2.0);

  real_t u[3][numNodes];
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i) u[c][i] = std::cos(static_cast<real_t>(c * numNodes + i));

  real_t f_serial[3][numNodes] = {{0}};
  QK::computeElasticStiffnessSumFact(X, u, f_serial, IdentityFlux{});

  Kokkos::View<float[8][3]> X_dev("X_dev");
  Kokkos::View<real_t*> u_dev("u_dev", 3 * numNodes);
  Kokkos::View<real_t*> f_dev("f_dev", 3 * numNodes);
  Kokkos::View<real_t*> F_dev("F_dev", 9 * numNodes);

  auto X_host = Kokkos::create_mirror_view(X_dev);
  auto u_host = Kokkos::create_mirror_view(u_dev);
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 3; ++j) X_host(i, j) = static_cast<float>(X[i][j]);
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i) u_host(c * numNodes + i) = u[c][i];
  Kokkos::deep_copy(X_dev, X_host);
  Kokkos::deep_copy(u_dev, u_host);

  launchElasticTeam<QK>(X_dev, u_dev, f_dev, F_dev);
  Kokkos::fence();

  auto f_host = Kokkos::create_mirror_view(f_dev);
  Kokkos::deep_copy(f_host, f_dev);
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i)
      EXPECT_NEAR(f_host(c * numNodes + i), f_serial[c][i], TOL_NUMERICAL) << "comp " << c << " node " << i;
}

// The cube is affine, so its Jacobian is the same at every quadrature point:
// the constant-Jacobian overload must reproduce the general one exactly. This is
// the property the solver relies on to skip the per-point jacobianTransformation
// on structured meshes.
TYPED_TEST(SumFactElasticTest, TeamConstantJacobianMatchesSerial) {
  using QK = TypeParam;
  constexpr int numNodes = QK::numNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.0, 0.5, 0.0, 2.0);

  real_t u[3][numNodes];
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i) u[c][i] = std::cos(static_cast<real_t>(c * numNodes + i));

  real_t f_serial[3][numNodes] = {{0}};
  QK::computeElasticStiffnessSumFact(X, u, f_serial, IdentityFlux{});

  // invJacobianTransformation accumulates into J, so it must start at zero.
  real_t J_inv[3][3] = {{0}};
  real_t const detJ = QK::invJacobianTransformation(0, 0, 0, X, J_inv);

  Kokkos::View<real_t*> geom_dev("geom_dev", 10);
  Kokkos::View<real_t*> u_dev("u_dev", 3 * numNodes);
  Kokkos::View<real_t*> f_dev("f_dev", 3 * numNodes);
  Kokkos::View<real_t*> F_dev("F_dev", 9 * numNodes);

  auto geom_host = Kokkos::create_mirror_view(geom_dev);
  auto u_host = Kokkos::create_mirror_view(u_dev);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) geom_host(a * 3 + b) = J_inv[a][b];
  geom_host(9) = detJ;
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i) u_host(c * numNodes + i) = u[c][i];
  Kokkos::deep_copy(geom_dev, geom_host);
  Kokkos::deep_copy(u_dev, u_host);

  launchElasticTeamConstJac<QK>(geom_dev, u_dev, f_dev, F_dev);
  Kokkos::fence();

  auto f_host = Kokkos::create_mirror_view(f_dev);
  Kokkos::deep_copy(f_host, f_dev);
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes; ++i)
      EXPECT_NEAR(f_host(c * numNodes + i), f_serial[c][i], TOL_NUMERICAL) << "comp " << c << " node " << i;
}

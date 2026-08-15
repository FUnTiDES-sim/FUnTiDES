#pragma once
#include "Qk_Hexahedron_Tensorial.h"
#include "common.h"

struct AlphaOne {
  KOKKOS_INLINE_FUNCTION real_t operator()(int, int, int) const { return real_t(1.0); }
};

template <typename QK>
inline void launchTeamVectorStreaming(Kokkos::View<float[8][3]> X_dev, Kokkos::View<real_t*> u_dev,
                                      Kokkos::View<real_t*> v_gemm_dev, Kokkos::View<real_t*> D_dev) {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamMember = TeamPolicy::member_type;

  TeamPolicy policy(1, Kokkos::AUTO);
  size_t scratch_bytes = QK::scratchBytesPerTeamStreaming();
  policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const TeamMember& team) {
        float X_local[8][3];
        for (int i = 0; i < 8; ++i)
          for (int j = 0; j < 3; ++j) X_local[i][j] = X_dev(i, j);

        QK::computeStiffnessOperatorTeamVectorStreaming(team, u_dev.data(), v_gemm_dev.data(), X_local, D_dev.data(),
                                                        AlphaOne{});
      });
}

template <typename QK>
inline void launchTeamVector(Kokkos::View<real_t*> u_dev, Kokkos::View<real_t*> v_dev, Kokkos::View<real_t*> W_dev,
                             Kokkos::View<real_t*> D_dev) {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamMember = TeamPolicy::member_type;

  TeamPolicy policy(1, Kokkos::AUTO);
  policy.set_scratch_size(0, Kokkos::PerTeam(QK::scratchBytesPerTeam()));

  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const TeamMember& team) {
        QK::computeStiffnessOperatorTeamVector(team, u_dev.data(), v_dev.data(), W_dev.data(), D_dev.data());
      });
}

template <typename QK_BASIS>
class TensorialGEMMKokkosTest : public ::testing::Test {};

TYPED_TEST_SUITE(TensorialGEMMKokkosTest, TestedBases);

TYPED_TEST(TensorialGEMMKokkosTest, TeamVectorStreamingMatchesSerial) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int numNodes = QK::numNodes;
  constexpr int n = QK::num1dNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -1.0, 0.5, 0.0, 2.0);

  real_t u_local[numNodes];
  for (int i = 0; i < numNodes; ++i) u_local[i] = std::cos(static_cast<real_t>(i));

  real_t v_serial[numNodes] = {0};
  QK::computeStiffnessTermSumFact(X, u_local, v_serial, [](int, int, int) { return real_t(1.0); });

  real_t D_flat[n * n] = {0};
  QK::fillDerivativeMatrix(D_flat);

  Kokkos::View<float[8][3]> X_dev("X_dev");
  auto X_host = Kokkos::create_mirror_view(X_dev);
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 3; ++j) X_host(i, j) = static_cast<float>(X[i][j]);
  Kokkos::deep_copy(X_dev, X_host);

  Kokkos::View<real_t*> u_dev("u_dev", numNodes);
  Kokkos::View<real_t*> v_gemm_dev("v_gemm_dev", numNodes);
  Kokkos::View<real_t*> D_dev("D_dev", n * n);

  auto u_host = Kokkos::create_mirror_view(u_dev);
  auto D_host = Kokkos::create_mirror_view(D_dev);
  for (int i = 0; i < numNodes; ++i) u_host(i) = u_local[i];
  for (int i = 0; i < n * n; ++i) D_host(i) = D_flat[i];
  Kokkos::deep_copy(u_dev, u_host);
  Kokkos::deep_copy(D_dev, D_host);

  launchTeamVectorStreaming<QK>(X_dev, u_dev, v_gemm_dev, D_dev);
  Kokkos::fence();

  auto v_gemm_host = Kokkos::create_mirror_view(v_gemm_dev);
  Kokkos::deep_copy(v_gemm_host, v_gemm_dev);

  for (int i = 0; i < numNodes; ++i) {
    EXPECT_NEAR(v_gemm_host(i), v_serial[i], TOL_NUMERICAL);
  }
}

// Production GEMM path: computeElementMetrics (precompute W) + non-streaming
// computeStiffnessOperatorTeamVector. This is the operator actually run by the
// solver's _Gemm path (the streaming variant above precomputes W in-kernel).
TYPED_TEST(TensorialGEMMKokkosTest, TeamVectorMatchesSerial) {
  using QK = Qk_Hexahedron_Tensorial_GEMM<typename TypeParam::BasisType>;
  constexpr int numNodes = QK::numNodes;
  constexpr int n = QK::num1dNodes;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -0.7, 0.3, 0.4, 1.6);

  real_t u_local[numNodes];
  for (int i = 0; i < numNodes; ++i) u_local[i] = std::sin(static_cast<real_t>(i) * real_t(0.5));

  // Reference: serial sum-factorization (alpha = 1).
  real_t v_serial[numNodes] = {0};
  QK::computeStiffnessTermSumFact(X, u_local, v_serial, [](int, int, int) { return real_t(1.0); });

  // Same operator, decomposed: precompute W = w*alpha*B (alpha = 1) and D on host.
  real_t W_flat[numNodes * 6] = {0};
  QK::computeElementMetrics(X, AlphaOne{}, W_flat);
  real_t D_flat[n * n] = {0};
  QK::fillDerivativeMatrix(D_flat);

  Kokkos::View<real_t*> u_dev("u_dev", numNodes);
  Kokkos::View<real_t*> v_dev("v_dev", numNodes);
  Kokkos::View<real_t*> W_dev("W_dev", numNodes * 6);
  Kokkos::View<real_t*> D_dev("D_dev", n * n);

  auto u_host = Kokkos::create_mirror_view(u_dev);
  auto W_host = Kokkos::create_mirror_view(W_dev);
  auto D_host = Kokkos::create_mirror_view(D_dev);
  for (int i = 0; i < numNodes; ++i) u_host(i) = u_local[i];
  for (int i = 0; i < numNodes * 6; ++i) W_host(i) = W_flat[i];
  for (int i = 0; i < n * n; ++i) D_host(i) = D_flat[i];
  Kokkos::deep_copy(u_dev, u_host);
  Kokkos::deep_copy(W_dev, W_host);
  Kokkos::deep_copy(D_dev, D_host);

  launchTeamVector<QK>(u_dev, v_dev, W_dev, D_dev);
  Kokkos::fence();

  auto v_host = Kokkos::create_mirror_view(v_dev);
  Kokkos::deep_copy(v_host, v_dev);

  for (int i = 0; i < numNodes; ++i) {
    EXPECT_NEAR(v_host(i), v_serial[i], TOL_NUMERICAL);
  }
}

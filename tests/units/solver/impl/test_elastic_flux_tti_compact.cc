/**
 * @file test_elastic_flux_tti_compact.cc
 * @brief The compact TTI flux must match the flux of the rotated 6x6 tensor.
 *
 * SEMsolver::computeTtiCompact() and SEMsolver::computeCMatrix() describe the same medium, so
 * flux::elasticFluxTtiCompact() and flux::elasticFluxTti() must give the same flux for any
 * gradient. The acousto-elastic validation case has all Thomsen parameters at zero, so this
 * test is what covers the anisotropic terms.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "Integrals.h"
#include "elastic_flux.h"
#include "model_unstruct.h"
#include "sem_solver.h"
#include "sem_solver_impl.h"

namespace solver {
namespace fe {
namespace test {

namespace {

using Integral = typename IntegralTypeSelector<1, IntegralType::MAKUTU>::type;
using ElasticSolver = SEMsolver<1, Integral, model::ModelUnstruct<float, int>, true, utils::enums::physicType::kElastic>;

constexpr float kJinv[3][3] = {{1.3f, 0.2f, -0.4f}, {0.1f, 0.9f, 0.3f}, {-0.2f, 0.5f, 1.1f}};

/// Deterministic values in [lo, hi).
struct Lcg {
  unsigned state = 12345u;
  float Uniform(float lo, float hi) {
    state = state * 1664525u + 1013904223u;
    return lo + (hi - lo) * static_cast<float>(state >> 8) / 16777216.0f;
  }
};

}  // namespace

TEST(ElasticFluxTtiCompact, MatchesRotatedTensor) {
  Lcg rng;
  for (int trial = 0; trial < 200; ++trial) {
    float const vp = rng.Uniform(2000.0f, 4000.0f);
    float const vs = rng.Uniform(800.0f, 1800.0f);
    float const rho = rng.Uniform(1800.0f, 2600.0f);
    float const delta = rng.Uniform(-0.1f, 0.2f);
    float const epsilon = rng.Uniform(0.0f, 0.3f);
    float const gamma = rng.Uniform(0.0f, 0.3f);
    float const phi = rng.Uniform(0.0f, 360.0f);
    float const theta = rng.Uniform(0.0f, 180.0f);

    float C[6][6];
    ElasticSolver::computeCMatrix(vp, vs, rho, delta, epsilon, gamma, phi, theta, C);
    float p[flux::kTtiCompactSize];
    ElasticSolver::computeTtiCompact(vp, vs, rho, delta, epsilon, gamma, phi, theta, p);

    float grad[3][3];
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) grad[a][b] = rng.Uniform(-1.0f, 1.0f);

    float ref[3][3];
    float compact[3][3];
    flux::elasticFluxTti(kJinv, C, grad, ref);
    flux::elasticFluxTtiCompact(kJinv, p, grad, compact);

    float scale = 0.0f;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) scale = std::fmax(scale, std::fabs(ref[a][b]));
    ASSERT_GT(scale, 0.0f);
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        EXPECT_NEAR(compact[a][b], ref[a][b], 1e-5f * scale) << "trial " << trial << ", entry (" << a << "," << b << ")";
  }
}

}  // namespace test
}  // namespace fe
}  // namespace solver

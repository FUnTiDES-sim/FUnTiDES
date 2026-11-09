#include <gtest/gtest.h>

#include "three_layer_case.h"

TEST(SEM2DLTS, ThreeLayersAdaptToVelocityAndMatchSmallGlobalStep) {
  using C = solver::validation::ThreeLayerResult;
  auto const result = solver::validation::runThreeLayerValidation();
  ASSERT_EQ(result.plan.coarse_elements.size(), 600u);
  ASSERT_EQ(result.plan.fine_elements.size(), 1200u);
  std::vector<bool> fine(C::ex * C::ez, false);
  for (int e : result.plan.fine_elements) fine[e] = true;
  for (int e = 0; e < C::ex * C::ez; ++e) {
    EXPECT_EQ(fine[e], result.vp[e] > 1500);
    EXPECT_NEAR(result.plan.element_steps[e].estimate, C::cfl * 50 / result.vp[e], 1e-7);
  }
  ASSERT_EQ(result.snapshots.size(), 5u);
  ASSERT_EQ(result.history.size(), 151u);
  double const initial_energy = result.history.front().energy;
  ASSERT_GT(initial_energy, 0);
  double max_drift = 0;
  for (auto const& sample : result.history) {
    ASSERT_TRUE(std::isfinite(sample.energy));
    ASSERT_TRUE(std::isfinite(sample.relative_error));
    max_drift = std::max(max_drift, std::abs(sample.energy / initial_energy - 1));
  }
  EXPECT_LT(max_drift, 1e-4);
  EXPECT_LT(result.history.back().relative_error, 0.03);
  std::cout << "[three layers] final relative error=" << result.history.back().relative_error
            << ", max relative energy drift=" << max_drift << '\n';
}

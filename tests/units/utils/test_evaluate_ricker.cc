#include <gtest/gtest.h>

#include <cmath>

#include "utils.h"

namespace utils
{
namespace test
{

class EvaluateRickerTest : public ::testing::Test
{
 protected:
  SolverUtils ricker;
  static constexpr float f0 = 5.0f;
  static constexpr float tpeak = 1.0f / f0;
  static constexpr float pi = static_cast<float>(M_PI);
  static constexpr float lam = (f0 * pi) * (f0 * pi);
};

// --- Order 0: g(t) = -(t - tpeak) * exp(-2*lam*(t-tpeak)^2) ---
// At t = tpeak => 0
TEST_F(EvaluateRickerTest, Order0_AtPeak_IsZero)
{
  EXPECT_NEAR(ricker.evaluateRicker(tpeak, f0, 0), 0.0f, 1e-7f);
}

// --- Order 1 (first derivative of Gaussian): -2*lam*(t-tpeak)*exp(...)  ---
// At t = tpeak => 0  (odd function around tpeak)
TEST_F(EvaluateRickerTest, Order1_AtPeak_IsZero)
{
  EXPECT_NEAR(ricker.evaluateRicker(tpeak, f0, 1), 0.0f, 1e-7f);
}

// Order 1 is antisymmetric: f(tpeak + dt) = -f(tpeak - dt)
TEST_F(EvaluateRickerTest, Order1_IsAntisymmetric)
{
  float dt = 0.02f;
  float vp = ricker.evaluateRicker(tpeak + dt, f0, 1);
  float vm = ricker.evaluateRicker(tpeak - dt, f0, 1);
  EXPECT_NEAR(vp, -vm, 1e-6f);
}

// --- Order 2 (Ricker wavelet / Mexican hat): peak at tpeak ---
// At t = tpeak: 2*lam*(0 - 1)*exp(0) = -2*lam
TEST_F(EvaluateRickerTest, Order2_AtPeak)
{
  float expected = -2.0f * lam;
  EXPECT_NEAR(ricker.evaluateRicker(tpeak, f0, 2), expected, 1e-3f);
}

// Order 2 is symmetric: f(tpeak + dt) = f(tpeak - dt)
TEST_F(EvaluateRickerTest, Order2_IsSymmetric)
{
  float dt = 0.03f;
  float vp = ricker.evaluateRicker(tpeak + dt, f0, 2);
  float vm = ricker.evaluateRicker(tpeak - dt, f0, 2);
  EXPECT_NEAR(vp, vm, 1e-6f);
}

// --- Order 3 (third derivative): antisymmetric ---
// At t = tpeak: (time_n - o_tpeak) = 0 => pulse = 0
TEST_F(EvaluateRickerTest, Order3_AtPeak_IsZero)
{
  EXPECT_NEAR(ricker.evaluateRicker(tpeak, f0, 3), 0.0f, 1e-7f);
}

// Order 3 is antisymmetric
TEST_F(EvaluateRickerTest, Order3_IsAntisymmetric)
{
  float dt = 0.02f;
  float vp = ricker.evaluateRicker(tpeak + dt, f0, 3);
  float vm = ricker.evaluateRicker(tpeak - dt, f0, 3);
  EXPECT_NEAR(vp, -vm, 1e-5f);
}

// --- Order 4 (fourth derivative): symmetric ---
// At t = tpeak: 4*lam^2 * (3 - 0 + 0) * exp(0) = 12*lam^2
TEST_F(EvaluateRickerTest, Order4_AtPeak)
{
  float expected = 12.0f * lam * lam;
  EXPECT_NEAR(ricker.evaluateRicker(tpeak, f0, 4), expected, 1e-1f);
}

// Order 4 is symmetric
TEST_F(EvaluateRickerTest, Order4_IsSymmetric)
{
  float dt = 0.02f;
  float vp = ricker.evaluateRicker(tpeak + dt, f0, 4);
  float vm = ricker.evaluateRicker(tpeak - dt, f0, 4);
  EXPECT_NEAR(vp, vm, 1e-4f);
}

// --- Window clipping: outside [-0.9*tpeak, 2.9*tpeak] returns 0 ---
TEST_F(EvaluateRickerTest, OutsideWindow_ReturnsZero)
{
  for (int ord = 0; ord <= 4; ++ord)
  {
    EXPECT_FLOAT_EQ(ricker.evaluateRicker(-0.9f * tpeak - 0.001f, f0, ord),
                    0.0f)
        << "order=" << ord << " before window";
    EXPECT_FLOAT_EQ(ricker.evaluateRicker(2.9f * tpeak + 0.001f, f0, ord),
                    0.0f)
        << "order=" << ord << " after window";
  }
}

// --- Derivative consistency: finite-diff of order n-1 ≈ order n ---
// Verifies that d/dt[ricker(t, order)] ≈ ricker(t, order+1) numerically.
TEST_F(EvaluateRickerTest, DerivativeConsistency_Order1vs2)
{
  float t = tpeak + 0.04f;
  float h = 1e-4f;
  float fd = (ricker.evaluateRicker(t + h, f0, 1) -
              ricker.evaluateRicker(t - h, f0, 1)) /
             (2.0f * h);
  float analytical = ricker.evaluateRicker(t, f0, 2);
  EXPECT_NEAR(fd, analytical, std::abs(analytical) * 0.01f);
}

TEST_F(EvaluateRickerTest, DerivativeConsistency_Order2vs3)
{
  float t = tpeak + 0.04f;
  float h = 1e-4f;
  float fd = (ricker.evaluateRicker(t + h, f0, 2) -
              ricker.evaluateRicker(t - h, f0, 2)) /
             (2.0f * h);
  float analytical = ricker.evaluateRicker(t, f0, 3);
  EXPECT_NEAR(fd, analytical, std::abs(analytical) * 0.01f);
}

TEST_F(EvaluateRickerTest, DerivativeConsistency_Order3vs4)
{
  float t = tpeak + 0.04f;
  float h = 1e-4f;
  float fd = (ricker.evaluateRicker(t + h, f0, 3) -
              ricker.evaluateRicker(t - h, f0, 3)) /
             (2.0f * h);
  float analytical = ricker.evaluateRicker(t, f0, 4);
  EXPECT_NEAR(fd, analytical, std::abs(analytical) * 0.01f);
}

}  // namespace test
}  // namespace utils

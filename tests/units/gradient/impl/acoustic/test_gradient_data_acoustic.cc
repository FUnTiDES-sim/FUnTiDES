#include <gtest/gtest.h>

#include "data_type.h"
#include "device_view_access.h"
#include "differentiator_data_acoustic.h"

namespace gradient {
namespace test {

using test_utils::fillFromHost;
using test_utils::toHost;

class GradientDataAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size = 40;

    pn = allocateVector<vectorReal>(size, "pn");
    qn = allocateVector<vectorReal>(size, "qn");
    qnPrev = allocateVector<vectorReal>(size, "qnPrev");
    qnPrevPrev = allocateVector<vectorReal>(size, "qnPrevPrev");

    gradKappa = allocateVector<vectorReal>(size, "gradKappa");
    gradBuoyancy = allocateVector<vectorReal>(size, "gradBuoyancy");

    fillFromHost(pn, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.1f;
    });
    fillFromHost(qn, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.2f;
    });
    fillFromHost(qnPrev, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.3f;
    });
    fillFromHost(qnPrevPrev, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.4f;
    });
    fillFromHost(gradKappa, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 1.0f;
    });
    fillFromHost(gradBuoyancy, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 2.0f;
    });
  }

  int size;
  vectorReal pn, qn, qnPrev, qnPrevPrev;
  vectorReal gradKappa, gradBuoyancy;
};

TEST_F(GradientDataAcousticTest, Construction) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic gradient(gradKappa, gradBuoyancy);

  EXPECT_NO_THROW((GradientDataAcoustic{fwd, bwd, gradient}));
}

TEST_F(GradientDataAcousticTest, GetForwardFieldReturnsPn) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = toHost(data.getForwardField(0));
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.1f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldZeroReturnsQn) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = toHost(data.getBackwardField(0));
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.2f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldOneReturnsQnPrev) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = toHost(data.getBackwardField(1));
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.3f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldTwoReturnsQnPrevPrev) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = toHost(data.getBackwardField(2));
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.4f);
}

TEST_F(GradientDataAcousticTest, GetGradientZeroReturnsKappa) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto g = toHost(data.getGradient(0));
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 1.0f);
}

TEST_F(GradientDataAcousticTest, GetGradientOneReturnsBuoyancy) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto g = toHost(data.getGradient(1));
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 2.0f);
}

TEST_F(GradientDataAcousticTest, PrintDoesNotThrow) {
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  EXPECT_NO_THROW(data.print());
}

}  // namespace test
}  // namespace gradient

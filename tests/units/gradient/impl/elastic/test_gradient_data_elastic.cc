#include <gtest/gtest.h>

#include "data_type.h"
#include "device_view_access.h"
#include "differentiator_data_elastic.h"

namespace gradient {
namespace test {

using test_utils::fillFromHost;
using test_utils::toHost;

class GradientDataElasticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size = 40;

    ux_n = allocateVector<vectorReal>(size, "ux_n");
    uy_n = allocateVector<vectorReal>(size, "uy_n");
    uz_n = allocateVector<vectorReal>(size, "uz_n");
    ux_dt2 = allocateVector<vectorReal>(size, "ux_dt2");
    uy_dt2 = allocateVector<vectorReal>(size, "uy_dt2");
    uz_dt2 = allocateVector<vectorReal>(size, "uz_dt2");

    gradRho = allocateVector<vectorReal>(size, "gradRho");
    gradLambda = allocateVector<vectorReal>(size, "gradLambda");
    gradMu = allocateVector<vectorReal>(size, "gradMu");

    fillFromHost(ux_n, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.1f;
    });
    fillFromHost(uy_n, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.2f;
    });
    fillFromHost(uz_n, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.3f;
    });
    fillFromHost(ux_dt2, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.4f;
    });
    fillFromHost(uy_dt2, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.5f;
    });
    fillFromHost(uz_dt2, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 0.6f;
    });
    fillFromHost(gradRho, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 1.0f;
    });
    fillFromHost(gradLambda, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 2.0f;
    });
    fillFromHost(gradMu, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 3.0f;
    });
  }

  int size;
  vectorReal ux_n, uy_n, uz_n;
  vectorReal ux_dt2, uy_dt2, uz_dt2;
  vectorReal gradRho, gradLambda, gradMu;
};

TEST_F(GradientDataElasticTest, Construction) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);

  EXPECT_NO_THROW((GradientDataElastic{fwd, bwd, grad}));
}

TEST_F(GradientDataElasticTest, GetForwardFieldReturnsUxN) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  auto field = toHost(data.getForwardField(0));
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.1f);
}

TEST_F(GradientDataElasticTest, GetForwardFieldReturnsUyN) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  auto field = toHost(data.getForwardField(1));
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.2f);
}

TEST_F(GradientDataElasticTest, GetBackwardFieldReturnsUxDt2) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  auto field = toHost(data.getBackwardField(3));
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.4f);
}

TEST_F(GradientDataElasticTest, GetGradientZeroReturnsRho) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  auto g = toHost(data.getGradient(0));
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 1.0f);
}

TEST_F(GradientDataElasticTest, GetGradientOneReturnsLambda) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  auto g = toHost(data.getGradient(1));
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 2.0f);
}

TEST_F(GradientDataElasticTest, GetGradientTwoReturnsMu) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  auto g = toHost(data.getGradient(2));
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 3.0f);
}

TEST_F(GradientDataElasticTest, PrintDoesNotThrow) {
  WavefieldViewForwardElastic fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic grad(gradRho, gradLambda, gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  EXPECT_NO_THROW(data.print());
}

}  // namespace test
}  // namespace gradient

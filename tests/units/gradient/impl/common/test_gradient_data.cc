#include <gtest/gtest.h>

#include "data_type.h"
#include "gradient_data.h"

namespace gradient
{
namespace test
{

// =============================================================================
// GradientDataAcoustic  (GradientData<physicType::kAcoustic>)
// =============================================================================

class GradientDataAcousticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    size = 40;

    pn   = allocateVector<VECTOR_REAL_VIEW>(size, "pn");
    qn   = allocateVector<VECTOR_REAL_VIEW>(size, "qn");
    qdt2 = allocateVector<VECTOR_REAL_VIEW>(size, "qdt2");

    gradKappa    = allocateVector<VECTOR_REAL_VIEW>(size, "gradKappa");
    gradBuoyancy = allocateVector<VECTOR_REAL_VIEW>(size, "gradBuoyancy");

    for (int i = 0; i < size; ++i)
    {
      pn(i)          = static_cast<float>(i) * 0.1f;
      qn(i)          = static_cast<float>(i) * 0.2f;
      qdt2(i)        = static_cast<float>(i) * 0.3f;
      gradKappa(i)   = static_cast<float>(i) * 1.0f;
      gradBuoyancy(i) = static_cast<float>(i) * 2.0f;
    }
  }

  int size;
  VECTOR_REAL_VIEW pn, qn, qdt2;
  VECTOR_REAL_VIEW gradKappa, gradBuoyancy;
};

TEST_F(GradientDataAcousticTest, Construction)
{
  WavefieldViewForwardAcoustic  fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qdt2);
  GradientAcoustic              gradient(gradKappa, gradBuoyancy);

  EXPECT_NO_THROW(
      (GradientDataAcoustic{fwd, bwd, gradient}));
}

TEST_F(GradientDataAcousticTest, GetForwardFieldReturnsPn)
{
  WavefieldViewForwardAcoustic  fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qdt2);
  GradientAcoustic              grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic          data(fwd, bwd, grad);

  auto field = data.getForwardField(0);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.1f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldZeroReturnsQn)
{
  WavefieldViewForwardAcoustic  fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qdt2);
  GradientAcoustic              grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic          data(fwd, bwd, grad);

  auto field = data.getBackwardField(0);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.2f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldOneReturnsQdt2)
{
  WavefieldViewForwardAcoustic  fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qdt2);
  GradientAcoustic              grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic          data(fwd, bwd, grad);

  auto field = data.getBackwardField(1);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.3f);
}

TEST_F(GradientDataAcousticTest, GetGradientZeroReturnsKappa)
{
  WavefieldViewForwardAcoustic  fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qdt2);
  GradientAcoustic              grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic          data(fwd, bwd, grad);

  auto g = data.getGradient(0);
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 1.0f);
}

TEST_F(GradientDataAcousticTest, GetGradientOneReturnsBuoyancy)
{
  WavefieldViewForwardAcoustic  fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qdt2);
  GradientAcoustic              grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic          data(fwd, bwd, grad);

  auto g = data.getGradient(1);
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 2.0f);
}

TEST_F(GradientDataAcousticTest, NumGradientsConstant)
{
  EXPECT_EQ(GradientDataAcoustic::kNumGradients, 2);
}

TEST_F(GradientDataAcousticTest, PrintDoesNotThrow)
{
  WavefieldViewForwardAcoustic  fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qdt2);
  GradientAcoustic              grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic          data(fwd, bwd, grad);

  EXPECT_NO_THROW(data.print());
}

// =============================================================================
// GradientDataElastic  (GradientData<physicType::kElastic>)
// =============================================================================

class GradientDataElasticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    size = 40;

    ux_n   = allocateVector<VECTOR_REAL_VIEW>(size, "ux_n");
    uy_n   = allocateVector<VECTOR_REAL_VIEW>(size, "uy_n");
    uz_n   = allocateVector<VECTOR_REAL_VIEW>(size, "uz_n");
    ux_dt2 = allocateVector<VECTOR_REAL_VIEW>(size, "ux_dt2");
    uy_dt2 = allocateVector<VECTOR_REAL_VIEW>(size, "uy_dt2");
    uz_dt2 = allocateVector<VECTOR_REAL_VIEW>(size, "uz_dt2");

    gradRho    = allocateVector<VECTOR_REAL_VIEW>(size, "gradRho");
    gradLambda = allocateVector<VECTOR_REAL_VIEW>(size, "gradLambda");
    gradMu     = allocateVector<VECTOR_REAL_VIEW>(size, "gradMu");

    for (int i = 0; i < size; ++i)
    {
      ux_n(i)    = static_cast<float>(i) * 0.1f;
      uy_n(i)    = static_cast<float>(i) * 0.2f;
      uz_n(i)    = static_cast<float>(i) * 0.3f;
      ux_dt2(i)  = static_cast<float>(i) * 0.4f;
      uy_dt2(i)  = static_cast<float>(i) * 0.5f;
      uz_dt2(i)  = static_cast<float>(i) * 0.6f;
      gradRho(i)    = static_cast<float>(i) * 1.0f;
      gradLambda(i) = static_cast<float>(i) * 2.0f;
      gradMu(i)     = static_cast<float>(i) * 3.0f;
    }
  }

  int size;
  VECTOR_REAL_VIEW ux_n, uy_n, uz_n;
  VECTOR_REAL_VIEW ux_dt2, uy_dt2, uz_dt2;
  VECTOR_REAL_VIEW gradRho, gradLambda, gradMu;
};

TEST_F(GradientDataElasticTest, Construction)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);

  EXPECT_NO_THROW(
      (GradientDataElastic{fwd, bwd, grad}));
}

TEST_F(GradientDataElasticTest, GetForwardFieldReturnsUxN)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);
  GradientDataElastic          data(fwd, bwd, grad);

  auto field = data.getForwardField(0);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.1f);
}

TEST_F(GradientDataElasticTest, GetForwardFieldReturnsUyN)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);
  GradientDataElastic          data(fwd, bwd, grad);

  auto field = data.getForwardField(1);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.2f);
}

TEST_F(GradientDataElasticTest, GetBackwardFieldReturnsUxDt2)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);
  GradientDataElastic          data(fwd, bwd, grad);

  auto field = data.getBackwardField(3);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.4f);
}

TEST_F(GradientDataElasticTest, GetGradientZeroReturnsRho)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);
  GradientDataElastic          data(fwd, bwd, grad);

  auto g = data.getGradient(0);
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 1.0f);
}

TEST_F(GradientDataElasticTest, GetGradientOneReturnsLambda)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);
  GradientDataElastic          data(fwd, bwd, grad);

  auto g = data.getGradient(1);
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 2.0f);
}

TEST_F(GradientDataElasticTest, GetGradientTwoReturnsMu)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);
  GradientDataElastic          data(fwd, bwd, grad);

  auto g = data.getGradient(2);
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 3.0f);
}

TEST_F(GradientDataElasticTest, NumGradientsConstant)
{
  EXPECT_EQ(GradientDataElastic::kNumGradients, 3);
}

TEST_F(GradientDataElasticTest, PrintDoesNotThrow)
{
  WavefieldViewForwardElastic  fwd(ux_n, uy_n, uz_n);
  WavefieldViewBackwardElastic bwd(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  GradientElastic              grad(gradRho, gradLambda, gradMu);
  GradientDataElastic          data(fwd, bwd, grad);

  EXPECT_NO_THROW(data.print());
}

}  // namespace test
}  // namespace gradient

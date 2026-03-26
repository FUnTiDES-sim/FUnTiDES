#include <gtest/gtest.h>

#include "data_type.h"
#include "differentiator_data_acoustic.h"

namespace gradient
{
namespace test
{

class GradientDataAcousticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    size = 40;

    pn = allocateVector<VECTOR_REAL_VIEW>(size, "pn");
    qn = allocateVector<VECTOR_REAL_VIEW>(size, "qn");
    qnPrev = allocateVector<VECTOR_REAL_VIEW>(size, "qnPrev");
    qnPrevPrev = allocateVector<VECTOR_REAL_VIEW>(size, "qnPrevPrev");

    gradKappa = allocateVector<VECTOR_REAL_VIEW>(size, "gradKappa");
    gradBuoyancy = allocateVector<VECTOR_REAL_VIEW>(size, "gradBuoyancy");

    for (int i = 0; i < size; ++i)
    {
      pn(i) = static_cast<float>(i) * 0.1f;
      qn(i) = static_cast<float>(i) * 0.2f;
      qnPrev(i) = static_cast<float>(i) * 0.3f;
      qnPrevPrev(i) = static_cast<float>(i) * 0.4f;
      gradKappa(i) = static_cast<float>(i) * 1.0f;
      gradBuoyancy(i) = static_cast<float>(i) * 2.0f;
    }
  }

  int size;
  VECTOR_REAL_VIEW pn, qn, qnPrev, qnPrevPrev;
  VECTOR_REAL_VIEW gradKappa, gradBuoyancy;
};

TEST_F(GradientDataAcousticTest, Construction)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic gradient(gradKappa, gradBuoyancy);

  EXPECT_NO_THROW((GradientDataAcoustic{fwd, bwd, gradient}));
}

TEST_F(GradientDataAcousticTest, GetForwardFieldReturnsPn)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = data.getForwardField(0);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.1f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldZeroReturnsQn)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = data.getBackwardField(0);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.2f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldOneReturnsQnPrev)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = data.getBackwardField(1);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.3f);
}

TEST_F(GradientDataAcousticTest, GetBackwardFieldTwoReturnsQnPrevPrev)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto field = data.getBackwardField(2);
  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.4f);
}

TEST_F(GradientDataAcousticTest, GetGradientZeroReturnsKappa)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto g = data.getGradient(0);
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 1.0f);
}

TEST_F(GradientDataAcousticTest, GetGradientOneReturnsBuoyancy)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  auto g = data.getGradient(1);
  ASSERT_EQ(g.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(g(i), static_cast<float>(i) * 2.0f);
}

TEST_F(GradientDataAcousticTest, PrintDoesNotThrow)
{
  WavefieldViewForwardAcoustic fwd(pn);
  WavefieldViewBackwardAcoustic bwd(qn, qnPrev, qnPrevPrev);
  GradientAcoustic grad(gradKappa, gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  EXPECT_NO_THROW(data.print());
}

}  // namespace test
}  // namespace gradient

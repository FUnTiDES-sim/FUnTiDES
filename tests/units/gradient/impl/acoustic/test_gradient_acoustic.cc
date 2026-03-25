#include <gtest/gtest.h>

#include "data_type.h"
#include "gradient_acoustic.h"

namespace gradient
{
namespace test
{

class GradientAcousticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    numElements = 8;
    numNodes    = 27;

    gradKappaElem    = allocateVector<VECTOR_REAL_VIEW>(numElements, "gradKappaElem");
    gradBuoyancyElem = allocateVector<VECTOR_REAL_VIEW>(numElements, "gradBuoyancyElem");
    gradKappaNode    = allocateVector<VECTOR_REAL_VIEW>(numNodes, "gradKappaNode");
    gradBuoyancyNode = allocateVector<VECTOR_REAL_VIEW>(numNodes, "gradBuoyancyNode");

    for (int i = 0; i < numElements; ++i)
    {
      gradKappaElem(i)    = static_cast<float>(i) * 1.5f;
      gradBuoyancyElem(i) = static_cast<float>(i) * 2.5f;
    }
    for (int i = 0; i < numNodes; ++i)
    {
      gradKappaNode(i)    = static_cast<float>(i) * 0.5f;
      gradBuoyancyNode(i) = static_cast<float>(i) * 1.0f;
    }
  }

  int numElements;
  int numNodes;
  VECTOR_REAL_VIEW gradKappaElem;
  VECTOR_REAL_VIEW gradBuoyancyElem;
  VECTOR_REAL_VIEW gradKappaNode;
  VECTOR_REAL_VIEW gradBuoyancyNode;
};

// --- Static constants ---

TEST_F(GradientAcousticTest, NumGradsConstant)
{
  EXPECT_EQ(GradientAcoustic::kNumGrads, 2);
}

TEST_F(GradientAcousticTest, GradsNamesConstants)
{
  EXPECT_STREQ(GradientAcoustic::kGradsNames[0], "gradKappa");
  EXPECT_STREQ(GradientAcoustic::kGradsNames[1], "gradBuoyancy");
}

TEST_F(GradientAcousticTest, ConstructorStoresViews)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);

  EXPECT_EQ(grad.m_gradKappa.extent(0),    static_cast<size_t>(numElements));
  EXPECT_EQ(grad.m_gradBuoyancy.extent(0), static_cast<size_t>(numElements));
}

TEST_F(GradientAcousticTest, GetNumGradientsReturns2)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);
  EXPECT_EQ(grad.getNumGradients(), 2);
}

TEST_F(GradientAcousticTest, GetGradientNamesCorrect)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);
  const char* const* names = grad.getGradientNames();
  EXPECT_STREQ(names[0], "gradKappa");
  EXPECT_STREQ(names[1], "gradBuoyancy");
}

TEST_F(GradientAcousticTest, GetGradientKappaElementBased)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);
  auto kappa = grad.getGradient(0);

  ASSERT_EQ(kappa.extent(0), static_cast<size_t>(numElements));
  for (int i = 0; i < numElements; ++i)
    EXPECT_FLOAT_EQ(kappa(i), static_cast<float>(i) * 1.5f);
}

TEST_F(GradientAcousticTest, GetGradientBuoyancyElementBased)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);
  auto buoyancy = grad.getGradient(1);

  ASSERT_EQ(buoyancy.extent(0), static_cast<size_t>(numElements));
  for (int i = 0; i < numElements; ++i)
    EXPECT_FLOAT_EQ(buoyancy(i), static_cast<float>(i) * 2.5f);
}

TEST_F(GradientAcousticTest, GetGradientKappaNodeBased)
{
  GradientAcoustic grad(gradKappaNode, gradBuoyancyNode);
  auto kappa = grad.getGradient(0);

  ASSERT_EQ(kappa.extent(0), static_cast<size_t>(numNodes));
  for (int i = 0; i < numNodes; ++i)
    EXPECT_FLOAT_EQ(kappa(i), static_cast<float>(i) * 0.5f);
}

TEST_F(GradientAcousticTest, GetGradientBuoyancyNodeBased)
{
  GradientAcoustic grad(gradKappaNode, gradBuoyancyNode);
  auto buoyancy = grad.getGradient(1);

  ASSERT_EQ(buoyancy.extent(0), static_cast<size_t>(numNodes));
  for (int i = 0; i < numNodes; ++i)
    EXPECT_FLOAT_EQ(buoyancy(i), static_cast<float>(i) * 1.0f);
}

TEST_F(GradientAcousticTest, GetGradientOutOfBoundsFallsBackToKappa)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);
  auto fallback = grad.getGradient(99);

  // Default case returns m_gradKappa
  EXPECT_EQ(fallback.extent(0), static_cast<size_t>(numElements));
  EXPECT_FLOAT_EQ(fallback(0), gradKappaElem(0));
}

TEST_F(GradientAcousticTest, ViewsAreShallowCopies)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);

  // Modify through original view — change should be visible through gradient
  gradKappaElem(0) = 999.0f;
  EXPECT_FLOAT_EQ(grad.m_gradKappa(0), 999.0f);
}

TEST_F(GradientAcousticTest, PrintDoesNotThrow)
{
  GradientAcoustic grad(gradKappaElem, gradBuoyancyElem);
  EXPECT_NO_THROW(grad.print());
}

TEST_F(GradientAcousticTest, PolymorphicInterface)
{
  std::unique_ptr<Gradient> grad =
      std::make_unique<GradientAcoustic>(gradKappaElem, gradBuoyancyElem);

  EXPECT_EQ(grad->getNumGradients(), 2);
  EXPECT_STREQ(grad->getGradientNames()[0], "gradKappa");
  EXPECT_EQ(grad->getGradient(0).extent(0), static_cast<size_t>(numElements));
}

}  // namespace test
}  // namespace gradient

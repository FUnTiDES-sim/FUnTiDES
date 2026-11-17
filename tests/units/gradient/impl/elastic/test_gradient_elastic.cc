#include <gtest/gtest.h>

#include "data_type.h"
#include "device_view_access.h"
#include "gradient_elastic.h"

namespace gradient {
namespace test {

using test_utils::fillFromHost;
using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

class GradientElasticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    numElements = 8;
    numNodes = 27;

    gradRhoElem = allocateVector<vectorReal>(numElements, "gradRhoElem");
    gradLambdaElem = allocateVector<vectorReal>(numElements, "gradLambdaElem");
    gradMuElem = allocateVector<vectorReal>(numElements, "gradMuElem");

    gradRhoNode = allocateVector<vectorReal>(numNodes, "gradRhoNode");
    gradLambdaNode = allocateVector<vectorReal>(numNodes, "gradLambdaNode");
    gradMuNode = allocateVector<vectorReal>(numNodes, "gradMuNode");

    fillFromHost(gradRhoElem, [&](const auto& h) {
      for (int i = 0; i < numElements; ++i) h(i) = static_cast<float>(i) * 1.0f;
    });
    fillFromHost(gradLambdaElem, [&](const auto& h) {
      for (int i = 0; i < numElements; ++i) h(i) = static_cast<float>(i) * 2.0f;
    });
    fillFromHost(gradMuElem, [&](const auto& h) {
      for (int i = 0; i < numElements; ++i) h(i) = static_cast<float>(i) * 3.0f;
    });
    fillFromHost(gradRhoNode, [&](const auto& h) {
      for (int i = 0; i < numNodes; ++i) h(i) = static_cast<float>(i) * 0.1f;
    });
    fillFromHost(gradLambdaNode, [&](const auto& h) {
      for (int i = 0; i < numNodes; ++i) h(i) = static_cast<float>(i) * 0.2f;
    });
    fillFromHost(gradMuNode, [&](const auto& h) {
      for (int i = 0; i < numNodes; ++i) h(i) = static_cast<float>(i) * 0.3f;
    });
  }

  int numElements;
  int numNodes;
  vectorReal gradRhoElem;
  vectorReal gradLambdaElem;
  vectorReal gradMuElem;
  vectorReal gradRhoNode;
  vectorReal gradLambdaNode;
  vectorReal gradMuNode;
};

// --- Static constants ---

TEST_F(GradientElasticTest, NumGradsConstant) { EXPECT_EQ(GradientElastic::kNumGrads, 3); }

TEST_F(GradientElasticTest, GetGradientNamesCorrect) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  EXPECT_EQ(grad.getGradientName(0), "gradRho");
  EXPECT_EQ(grad.getGradientName(1), "gradLambda");
  EXPECT_EQ(grad.getGradientName(2), "gradMu");
}

// --- Constructor ---

TEST_F(GradientElasticTest, ConstructorStoresViews) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);

  // Verify via public interface
  EXPECT_EQ(grad.getGradient(0).extent(0), static_cast<size_t>(numElements));
  EXPECT_EQ(grad.getGradient(1).extent(0), static_cast<size_t>(numElements));
  EXPECT_EQ(grad.getGradient(2).extent(0), static_cast<size_t>(numElements));
}

// --- Virtual interface: getNumGradients ---

TEST_F(GradientElasticTest, GetNumGradientsReturns3) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  EXPECT_EQ(grad.getNumGradients(), 3);
}

// --- Virtual interface: getGradientNames ---

TEST_F(GradientElasticTest, GetGradientNameByIndex) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  EXPECT_EQ(grad.getGradientName(0), "gradRho");
  EXPECT_EQ(grad.getGradientName(1), "gradLambda");
  EXPECT_EQ(grad.getGradientName(2), "gradMu");
}

// --- Virtual interface: getGradient ---

TEST_F(GradientElasticTest, GetGradientRhoElementBased) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  auto rho = toHost(grad.getGradient(0));

  ASSERT_EQ(rho.extent(0), static_cast<size_t>(numElements));
  for (int i = 0; i < numElements; ++i) EXPECT_FLOAT_EQ(rho(i), static_cast<float>(i) * 1.0f);
}

TEST_F(GradientElasticTest, GetGradientLambdaElementBased) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  auto lambda = toHost(grad.getGradient(1));

  ASSERT_EQ(lambda.extent(0), static_cast<size_t>(numElements));
  for (int i = 0; i < numElements; ++i) EXPECT_FLOAT_EQ(lambda(i), static_cast<float>(i) * 2.0f);
}

TEST_F(GradientElasticTest, GetGradientMuElementBased) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  auto mu = toHost(grad.getGradient(2));

  ASSERT_EQ(mu.extent(0), static_cast<size_t>(numElements));
  for (int i = 0; i < numElements; ++i) EXPECT_FLOAT_EQ(mu(i), static_cast<float>(i) * 3.0f);
}

TEST_F(GradientElasticTest, GetGradientRhoNodeBased) {
  GradientElastic grad(gradRhoNode, gradLambdaNode, gradMuNode);
  auto rho = toHost(grad.getGradient(0));

  ASSERT_EQ(rho.extent(0), static_cast<size_t>(numNodes));
  for (int i = 0; i < numNodes; ++i) EXPECT_FLOAT_EQ(rho(i), static_cast<float>(i) * 0.1f);
}

TEST_F(GradientElasticTest, GetGradientLambdaNodeBased) {
  GradientElastic grad(gradRhoNode, gradLambdaNode, gradMuNode);
  auto lambda = toHost(grad.getGradient(1));

  ASSERT_EQ(lambda.extent(0), static_cast<size_t>(numNodes));
  for (int i = 0; i < numNodes; ++i) EXPECT_FLOAT_EQ(lambda(i), static_cast<float>(i) * 0.2f);
}

TEST_F(GradientElasticTest, GetGradientMuNodeBased) {
  GradientElastic grad(gradRhoNode, gradLambdaNode, gradMuNode);
  auto mu = toHost(grad.getGradient(2));

  ASSERT_EQ(mu.extent(0), static_cast<size_t>(numNodes));
  for (int i = 0; i < numNodes; ++i) EXPECT_FLOAT_EQ(mu(i), static_cast<float>(i) * 0.3f);
}

TEST_F(GradientElasticTest, GetGradientOutOfBoundsFallsBackToRho) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  auto fallback = toHost(grad.getGradient(99));

  // Default case returns m_gradRho
  EXPECT_EQ(fallback.extent(0), static_cast<size_t>(numElements));
  EXPECT_FLOAT_EQ(fallback(0), getAt(gradRhoElem, 0));
}

// --- Shared-view (shallow copy) semantics ---

TEST_F(GradientElasticTest, RhoViewIsShallowCopy) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  setAt(gradRhoElem, 777.0f, 0);
  EXPECT_FLOAT_EQ(getAt(grad.getGradient(0), 0), 777.0f);
}

TEST_F(GradientElasticTest, LambdaViewIsShallowCopy) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  setAt(gradLambdaElem, 888.0f, 2);
  EXPECT_FLOAT_EQ(getAt(grad.getGradient(1), 2), 888.0f);
}

// --- print() ---

TEST_F(GradientElasticTest, PrintDoesNotThrow) {
  GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
  EXPECT_NO_THROW(grad.print());
}

// --- Polymorphic usage via base pointer ---

TEST_F(GradientElasticTest, PolymorphicInterface) {
  std::unique_ptr<Gradient> grad = std::make_unique<GradientElastic>(gradRhoElem, gradLambdaElem, gradMuElem);

  EXPECT_EQ(grad->getNumGradients(), 3);
  EXPECT_EQ(grad->getGradientName(0), "gradRho");
  EXPECT_EQ(grad->getGradient(0).extent(0), static_cast<size_t>(numElements));
}

}  // namespace test
}  // namespace gradient

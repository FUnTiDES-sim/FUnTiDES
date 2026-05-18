#include <gtest/gtest.h>

#include <stdexcept>

#include "cartesian_unstruct_builder.h"

namespace model {
namespace test {

static CartesianParams<float, int> makeParams(int order, int ex, int ey, int ez, float lx, float ly, float lz,
                                              bool onNodes, bool elastic) {
  return CartesianParams<float, int>(order, ex, ey, ez, lx, ly, lz, onNodes, elastic);
}

// ============================================================================
// CartesianParams
// ============================================================================

TEST(CartesianParamsTest, DefaultConstructorSetsAEFalse) {
  CartesianParams<float, int> p;
  EXPECT_FALSE(p.isAcoustoElastic);
  EXPECT_FLOAT_EQ(p.acoustoElasticBoundaryZ, 0.0f);
}

TEST(CartesianParamsTest, ParameterizedConstructorSetsFields) {
  auto p = makeParams(2, 4, 5, 6, 100.0f, 200.0f, 300.0f, true, false);
  EXPECT_EQ(p.order, 2);
  EXPECT_EQ(p.ex, 4);
  EXPECT_EQ(p.ey, 5);
  EXPECT_EQ(p.ez, 6);
  EXPECT_FLOAT_EQ(p.lx, 100.0f);
  EXPECT_FLOAT_EQ(p.ly, 200.0f);
  EXPECT_FLOAT_EQ(p.lz, 300.0f);
  EXPECT_TRUE(p.isModelOnNodes);
  EXPECT_FALSE(p.isElastic);
  EXPECT_FALSE(p.isAcoustoElastic);
  EXPECT_EQ(p.model_file, "");
}

// ============================================================================
// CartesianUnstructBuilder — invalid order
// ============================================================================

TEST(CartesianUnstructBuilderTest, InvalidOrderThrows) {
  auto p = makeParams(10, 1, 1, 1, 10.0f, 10.0f, 10.0f, false, false);
  EXPECT_THROW((CartesianUnstructBuilder<float, int>(p)), std::runtime_error);
}

// ============================================================================
// CartesianUnstructBuilder — acoustic elements
// ============================================================================

TEST(CartesianUnstructBuilderTest, AcousticElementsBuildsValidModel) {
  auto p = makeParams(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, false, false);
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_EQ(model->getNumberOfElements(), 8);
  EXPECT_EQ(model->getOrder(), 1);
  EXPECT_FALSE(model->isModelOnNodes());
  EXPECT_FALSE(model->isElastic());
  EXPECT_FLOAT_EQ(model->getModelVpOnElement(0), 1500.0f);
  EXPECT_FLOAT_EQ(model->getModelRhoOnElement(0), 1.0f);
}

// ============================================================================
// CartesianUnstructBuilder — acoustic nodes
// ============================================================================

TEST(CartesianUnstructBuilderTest, AcousticNodesBuildsValidModel) {
  auto p = makeParams(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, true, false);
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_TRUE(model->isModelOnNodes());
  EXPECT_FALSE(model->isElastic());
  EXPECT_FLOAT_EQ(model->getModelVpOnNodes(0), 1500.0f);
  EXPECT_FLOAT_EQ(model->getModelRhoOnNodes(0), 1.0f);
}

// ============================================================================
// CartesianUnstructBuilder — elastic elements
// ============================================================================

TEST(CartesianUnstructBuilderTest, ElasticElementsBuildsValidModel) {
  auto p = makeParams(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, false, true);
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_TRUE(model->isElastic());
  EXPECT_FLOAT_EQ(model->getModelVsOnElement(0), 755.0f);
  EXPECT_FLOAT_EQ(model->getModelDeltaOnElement(0), 0.0f);
  EXPECT_FLOAT_EQ(model->getModelEpsilonOnElement(0), 0.0f);
  EXPECT_FLOAT_EQ(model->getModelGammaOnElement(0), 0.0f);
  EXPECT_EQ(model->getModelThetaOnElement(0), 0);
  EXPECT_EQ(model->getModelPhiOnElement(0), 0);
}

// ============================================================================
// CartesianUnstructBuilder — elastic nodes
// ============================================================================

TEST(CartesianUnstructBuilderTest, ElasticNodesBuildsValidModel) {
  auto p = makeParams(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, true, true);
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_TRUE(model->isModelOnNodes());
  EXPECT_TRUE(model->isElastic());
  EXPECT_FLOAT_EQ(model->getModelVsOnNodes(0), 755.0f);
  EXPECT_FLOAT_EQ(model->getModelDeltaOnNodes(0), 0.0f);
  EXPECT_FLOAT_EQ(model->getModelEpsilonOnNodes(0), 0.0f);
  EXPECT_FLOAT_EQ(model->getModelGammaOnNodes(0), 0.0f);
  EXPECT_EQ(model->getModelThetaOnNodes(0), 0);
  EXPECT_EQ(model->getModelPhiOnNodes(0), 0);
}

// ============================================================================
// CartesianUnstructBuilder — AcoustoElastic elements
// ============================================================================

TEST(CartesianUnstructBuilderTest, AcoustoElasticElementsBuildsValidModel) {
  CartesianParams<float, int> p(1, 2, 2, 4, 100.0f, 100.0f, 200.0f, false, true);
  p.isAcoustoElastic = true;
  p.acoustoElasticBoundaryZ = 100.0f;
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_EQ(model->getNumberOfElements(), 2 * 2 * 4);
  EXPECT_FLOAT_EQ(model->domainSize(2), 200.0f);
}

// ============================================================================
// CartesianUnstructBuilder — AcoustoElastic nodes
// ============================================================================

TEST(CartesianUnstructBuilderTest, AcoustoElasticNodesBuildsValidModel) {
  CartesianParams<float, int> p(1, 2, 2, 4, 100.0f, 100.0f, 200.0f, true, true);
  p.isAcoustoElastic = true;
  p.acoustoElasticBoundaryZ = 100.0f;
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_TRUE(model->isModelOnNodes());
}

// ============================================================================
// CartesianUnstructBuilder — face connectivity built automatically
// ============================================================================

TEST(CartesianUnstructBuilderTest, GetModelHasFaceConnectivity) {
  auto p = makeParams(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, false, false);
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_GT(model->getNumberOfFaces(), 0);
}

// ============================================================================
// CartesianUnstructBuilder — freeSurface flag
// ============================================================================

TEST(CartesianUnstructBuilderTest, FreeSurfaceFlagProducesValidModel) {
  auto p = makeParams(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, false, false);
  CartesianUnstructBuilder<float, int> b(p);
  auto model_fs = b.getModel(true);
  auto model_no_fs = b.getModel(false);
  ASSERT_NE(model_fs, nullptr);
  ASSERT_NE(model_no_fs, nullptr);
  EXPECT_EQ(model_fs->getNumberOfElements(), model_no_fs->getNumberOfElements());
}

// ============================================================================
// CartesianUnstructBuilder — higher order
// ============================================================================

TEST(CartesianUnstructBuilderTest, Order3BuildsCorrectNodeCount) {
  auto p = makeParams(3, 2, 2, 2, 100.0f, 100.0f, 100.0f, false, false);
  CartesianUnstructBuilder<float, int> b(p);
  auto model = b.getModel(true);
  ASSERT_NE(model, nullptr);
  EXPECT_EQ(model->getOrder(), 3);
  int nx = 2 * 3 + 1;  // 7
  EXPECT_EQ(model->getNumberOfNodes(), nx * nx * nx);
}

// ============================================================================
// CartesianUnstructBuilder — multiple getModel calls return different instances
// ============================================================================

TEST(CartesianUnstructBuilderTest, MultipleGetModelCallsDifferentInstances) {
  auto p = makeParams(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, false, false);
  CartesianUnstructBuilder<float, int> b(p);
  auto m1 = b.getModel(true);
  auto m2 = b.getModel(true);
  ASSERT_NE(m1.get(), m2.get());
}

// ============================================================================
// CartesianUnstructBuilder — model_file + isModelOnNodes=true throws
// ============================================================================

TEST(CartesianUnstructBuilderTest, ModelFileWithOnNodesModeThrows) {
  CartesianParams<float, int> p(1, 2, 2, 2, 100.0f, 100.0f, 100.0f, true, false);
  p.model_file = "dummy.txt";
  EXPECT_THROW((CartesianUnstructBuilder<float, int>(p)), std::runtime_error);
}

}  // namespace test
}  // namespace model

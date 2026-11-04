#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "differentiator.h"
#include "differentiator_factory.h"

namespace gradient {
namespace testing {

// Helper function to provide a valid baseline configuration
std::unique_ptr<Differentiator> createDefaultTestDifferentiator(
    utils::enums::physicType physic = utils::enums::physicType::kAcoustic,
    utils::enums::meshType mesh = utils::enums::meshType::kStruct,
    utils::enums::modelLocationType loc = utils::enums::modelLocationType::kOnNodes, int order = 1) {
  return createDifferentiator(utils::enums::implemType::kMakutu, mesh, loc, physic, order);
}

TEST(DifferentiatorFactoryTest, CreatesAcousticStructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(utils::enums::physicType::kAcoustic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnNodes, 1);

  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticUnstructuredOnElements) {
  // Casting to an unknown value triggers the 'else' branch for isModelOnNodes
  auto diff = createDefaultTestDifferentiator(utils::enums::physicType::kElastic, utils::enums::meshType::kUnstruct,
                                              static_cast<utils::enums::modelLocationType>(99), 1);

  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, ThrowsOnUnsupportedPolynomialOrder) {
  int unsupported_order = 999;

  EXPECT_THROW(
      {
        createDefaultTestDifferentiator(utils::enums::physicType::kAcoustic, utils::enums::meshType::kStruct,
                                        utils::enums::modelLocationType::kOnNodes, unsupported_order);
      },
      std::runtime_error);
}

TEST(DifferentiatorFactoryTest, ThrowsOnUnsupportedImplementationType) {
  auto invalid_implem = static_cast<utils::enums::implemType>(999);

  EXPECT_THROW(
      {
        createDifferentiator(invalid_implem, utils::enums::meshType::kStruct, utils::enums::modelLocationType::kOnNodes,
                             utils::enums::physicType::kAcoustic, 1);
      },
      std::runtime_error);
}

TEST(DifferentiatorFactoryTest, ThrowsOnUnknownPhysicsType) {
  auto invalid_physic = static_cast<utils::enums::physicType>(999);

  EXPECT_THROW(
      {
        createDefaultTestDifferentiator(invalid_physic, utils::enums::meshType::kStruct,
                                        utils::enums::modelLocationType::kOnNodes, 1);
      },
      std::runtime_error);
}
TEST(DifferentiatorFactoryTest, CreatesAcousticStructuredOnElements) {
  auto diff = createDefaultTestDifferentiator(
      utils::enums::physicType::kAcoustic, utils::enums::meshType::kStruct,
      static_cast<utils::enums::modelLocationType>(99),  // Passe dans le 'else' de isModelOnNodes
      1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticStructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(utils::enums::physicType::kElastic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnNodes, 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticStructuredOnElements) {
  auto diff = createDefaultTestDifferentiator(utils::enums::physicType::kElastic, utils::enums::meshType::kStruct,
                                              static_cast<utils::enums::modelLocationType>(99), 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesAcousticUnstructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(utils::enums::physicType::kAcoustic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnNodes, 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesAcousticUnstructuredOnElements) {
  auto diff = createDefaultTestDifferentiator(utils::enums::physicType::kAcoustic, utils::enums::meshType::kUnstruct,
                                              static_cast<utils::enums::modelLocationType>(99), 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticUnstructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(utils::enums::physicType::kElastic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnNodes, 1);
  EXPECT_NE(diff, nullptr);
}
}  // namespace testing
}  // namespace gradient

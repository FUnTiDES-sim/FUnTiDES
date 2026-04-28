#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "differentiator.h"
#include "differentiator_factory.h"

namespace gradient {
namespace testing {

namespace feenum = utils::enums;

// Helper function to provide a valid baseline configuration
std::unique_ptr<Differentiator> createDefaultTestDifferentiator(
    feenum::physicType physic = feenum::physicType::kAcoustic, feenum::meshType mesh = feenum::meshType::kStruct,
    feenum::modelLocationType loc = feenum::modelLocationType::kOnNodes, int order = 1) {
  return createDifferentiator(feenum::implemType::kMakutu, mesh, loc, physic, order);
}

TEST(DifferentiatorFactoryTest, CreatesAcousticStructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(feenum::physicType::kAcoustic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnNodes, 1);

  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticUnstructuredOnElements) {
  // Casting to an unknown value triggers the 'else' branch for isModelOnNodes
  auto diff = createDefaultTestDifferentiator(feenum::physicType::kElastic, feenum::meshType::kUnstruct,
                                              static_cast<feenum::modelLocationType>(99), 1);

  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, ThrowsOnUnsupportedPolynomialOrder) {
  int unsupported_order = 999;

  EXPECT_THROW(
      {
        createDefaultTestDifferentiator(feenum::physicType::kAcoustic, feenum::meshType::kStruct,
                                        feenum::modelLocationType::kOnNodes, unsupported_order);
      },
      std::runtime_error);
}

TEST(DifferentiatorFactoryTest, ThrowsOnUnsupportedImplementationType) {
  auto invalid_implem = static_cast<feenum::implemType>(999);

  EXPECT_THROW(
      {
        createDifferentiator(invalid_implem, feenum::meshType::kStruct, feenum::modelLocationType::kOnNodes,
                             feenum::physicType::kAcoustic, 1);
      },
      std::runtime_error);
}

TEST(DifferentiatorFactoryTest, ThrowsOnUnknownPhysicsType) {
  auto invalid_physic = static_cast<feenum::physicType>(999);

  EXPECT_THROW(
      {
        createDefaultTestDifferentiator(invalid_physic, feenum::meshType::kStruct, feenum::modelLocationType::kOnNodes,
                                        1);
      },
      std::runtime_error);
}
TEST(DifferentiatorFactoryTest, CreatesAcousticStructuredOnElements) {
  auto diff = createDefaultTestDifferentiator(
      feenum::physicType::kAcoustic, feenum::meshType::kStruct,
      static_cast<feenum::modelLocationType>(99),  // Passe dans le 'else' de isModelOnNodes
      1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticStructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(feenum::physicType::kElastic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnNodes, 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticStructuredOnElements) {
  auto diff = createDefaultTestDifferentiator(feenum::physicType::kElastic, feenum::meshType::kStruct,
                                              static_cast<feenum::modelLocationType>(99), 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesAcousticUnstructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(feenum::physicType::kAcoustic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnNodes, 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesAcousticUnstructuredOnElements) {
  auto diff = createDefaultTestDifferentiator(feenum::physicType::kAcoustic, feenum::meshType::kUnstruct,
                                              static_cast<feenum::modelLocationType>(99), 1);
  EXPECT_NE(diff, nullptr);
}

TEST(DifferentiatorFactoryTest, CreatesElasticUnstructuredOnNodes) {
  auto diff = createDefaultTestDifferentiator(feenum::physicType::kElastic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnNodes, 1);
  EXPECT_NE(diff, nullptr);
}
}  // namespace testing
}  // namespace gradient

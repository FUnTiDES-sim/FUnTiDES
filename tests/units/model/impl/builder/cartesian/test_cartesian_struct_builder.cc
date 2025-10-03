#include <gtest/gtest.h>
#include <memory>
#include "cartesian_struct_builder.h"

namespace model {
namespace test {

// Parameterized test fixture for different orders and isModelOnNodes values
class CartesianStructBuilderFixture : public ::testing::TestWithParam<std::tuple<int, bool>> {
protected:
  static constexpr int ex = 10;
  static constexpr float hx = 0.1f;
  static constexpr int ey = 20;
  static constexpr float hy = 0.2f;
  static constexpr int ez = 30;
  static constexpr float hz = 0.3f;

  std::unique_ptr<ModelBuilderBase<float, int>> createBuilder(int order, bool isModelOnNodes) {
    switch(order) {
      case 1:
        return std::make_unique<CartesianStructBuilder<float, int, 1>>(ex, hx, ey, hy, ez, hz, isModelOnNodes);
      case 2:
        return std::make_unique<CartesianStructBuilder<float, int, 2>>(ex, hx, ey, hy, ez, hz, isModelOnNodes);
      case 3:
        return std::make_unique<CartesianStructBuilder<float, int, 3>>(ex, hx, ey, hy, ez, hz, isModelOnNodes);
      case 4:
        return std::make_unique<CartesianStructBuilder<float, int, 4>>(ex, hx, ey, hy, ez, hz, isModelOnNodes);
      default:
        return nullptr;
    }
  }
};

// Instantiate the parameterized tests with different orders and isModelOnNodes values
INSTANTIATE_TEST_SUITE_P(
    CartesianStructBuilderTests,
    CartesianStructBuilderFixture,
    ::testing::Combine(
        ::testing::Values(1, 2, 3, 4),
        ::testing::Bool()
    )
);


// Test constructor and getModel for all orders and isModelOnNodes values, and check resulting model values
TEST_P(CartesianStructBuilderFixture, GetModelReturnsValidModel) {
  // Prepare
  int order = std::get<0>(GetParam());
  bool isModelOnNodes = std::get<1>(GetParam());

  // Act
  auto builder = createBuilder(order, isModelOnNodes);
  ASSERT_NE(builder, nullptr);
  auto model = builder->getModel();
  ASSERT_NE(model, nullptr);

  // Assert
  EXPECT_EQ(model->getOrder(), order);
  EXPECT_EQ(model->isModelOnNodes(), isModelOnNodes);
  EXPECT_EQ(model->getNumberOfElements(), 10 * 20 * 30);
  EXPECT_EQ(model->getNumberOfNodes(), (10 * order + 1) * (20 * order + 1) * (30 * order + 1));
  EXPECT_FLOAT_EQ(model->domainSize(0), 10 * 0.1f);
  EXPECT_FLOAT_EQ(model->domainSize(1), 20 * 0.2f);
  EXPECT_FLOAT_EQ(model->domainSize(2), 30 * 0.3f);
}

// Test multiple calls return different instances
TEST_P(CartesianStructBuilderFixture, MultipleCallsReturnDifferentInstances) {
    // Prepare
    int order = std::get<0>(GetParam());
    bool isModelOnNodes = std::get<1>(GetParam());

    // Act
    auto builder = createBuilder(order, isModelOnNodes);
    ASSERT_NE(builder, nullptr);

    // Assert
    auto model1 = builder->getModel();
    auto model2 = builder->getModel();
    ASSERT_NE(model1, nullptr);
    ASSERT_NE(model2, nullptr);
    EXPECT_NE(model1.get(), model2.get());
}

// Test polymorphic behavior
TEST_P(CartesianStructBuilderFixture, PolymorphicBehavior) {
    int order = std::get<0>(GetParam());
    bool isModelOnNodes = std::get<1>(GetParam());
    auto builder = createBuilder(order, isModelOnNodes);

    ASSERT_NE(builder, nullptr);
    ModelBuilderBase<float, int>* base_ptr = builder.get();
    auto model = base_ptr->getModel();
    
    ASSERT_NE(model, nullptr);
}

} // namespace test
} // namespace model
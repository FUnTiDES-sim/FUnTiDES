#include <gtest/gtest.h>

#include <memory>

#include "cartesian_struct_builder.h"

namespace model
{
namespace test
{

template <typename T>
class CartesianStructInputs : public ::testing::Test
{
 protected:
  static constexpr int ex = 10;
  static constexpr float hx = 0.1f;
  static constexpr int ey = 20;
  static constexpr float hy = 0.2f;
  static constexpr int ez = 30;
  static constexpr float hz = 0.3f;
};

template <int Order, bool IsModelOnNodes>
struct BuilderConfig
{
  using Type = CartesianStructBuilder<float, int, Order>;
  static constexpr int order = Order;
  static constexpr bool isModelOnNodes = IsModelOnNodes;
};

// Define all combinations of Order (1-4) and isModelOnNodes (true/false)
using BuilderTypes =
    ::testing::Types<BuilderConfig<1, true>, BuilderConfig<1, false>,
                     BuilderConfig<2, true>, BuilderConfig<2, false>,
                     BuilderConfig<3, true>, BuilderConfig<3, false>,
                     BuilderConfig<4, true>, BuilderConfig<4, false>>;

TYPED_TEST_SUITE(CartesianStructInputs, BuilderTypes);

// Test constructor and getModel for all orders and isModelOnNodes values, and
// check resulting model values
TYPED_TEST(CartesianStructInputs, GetModelReturnsValidModel)
{
  // Prepare
  constexpr int order = TypeParam::order;
  constexpr bool isModelOnNodes = TypeParam::isModelOnNodes;

  // Act
  typename TypeParam::Type builder(this->ex, this->hx, this->ey, this->hy,
                                   this->ez, this->hz, isModelOnNodes);
  auto model = builder.getModel();

  // Assert
  ASSERT_NE(model, nullptr);
  EXPECT_EQ(model->getOrder(), order);
  EXPECT_EQ(model->isModelOnNodes(), isModelOnNodes);
  EXPECT_EQ(model->getNumberOfElements(), 10 * 20 * 30);
  EXPECT_EQ(model->getNumberOfNodes(),
            (10 * order + 1) * (20 * order + 1) * (30 * order + 1));
  EXPECT_FLOAT_EQ(model->domainSize(0), 10 * 0.1f);
  EXPECT_FLOAT_EQ(model->domainSize(1), 20 * 0.2f);
  EXPECT_FLOAT_EQ(model->domainSize(2), 30 * 0.3f);
}

// Test multiple calls return different instances
TYPED_TEST(CartesianStructInputs, MultipleCallsReturnDifferentInstances)
{
  // Prepare
  constexpr bool isModelOnNodes = TypeParam::isModelOnNodes;

  // Act
  typename TypeParam::Type builder(this->ex, this->hx, this->ey, this->hy,
                                   this->ez, this->hz, isModelOnNodes);

  // Assert
  auto model1 = builder.getModel();
  auto model2 = builder.getModel();
  ASSERT_NE(model1, nullptr);
  ASSERT_NE(model2, nullptr);
  EXPECT_NE(model1.get(), model2.get());
}

// Test polymorphic behavior
TYPED_TEST(CartesianStructInputs, PolymorphicBehavior)
{
  // Prepare
  constexpr bool isModelOnNodes = TypeParam::isModelOnNodes;

  // Act
  typename TypeParam::Type builder(this->ex, this->hx, this->ey, this->hy,
                                   this->ez, this->hz, isModelOnNodes);
  ModelBuilderBase<float, int>* base_ptr = &builder;
  auto model = base_ptr->getModel();

  // Assert
  ASSERT_NE(model, nullptr);
}

}  // namespace test
}  // namespace model
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
  static constexpr float lx = 1;
  static constexpr int ey = 20;
  static constexpr float ly = 4;
  static constexpr int ez = 30;
  static constexpr float lz = 9;
};

template <int Order, bool IsModelOnNodes, bool IsElastic>
struct BuilderConfig
{
  using Type = CartesianStructBuilder<float, int, Order>;
  static constexpr int order = Order;
  static constexpr bool isModelOnNodes = IsModelOnNodes;
  static constexpr bool isElastic = IsElastic;
};

// Define all combinations of Order (1-4) and isModelOnNodes (true/false)
using BuilderTypes =
    ::testing::Types<BuilderConfig<1, true, true>, BuilderConfig<1, false, false>,
                     BuilderConfig<1, true, false>, BuilderConfig<1, false, true>,
                     BuilderConfig<2, true, true>, BuilderConfig<2, false, false>,
                     BuilderConfig<2, true, false>, BuilderConfig<2, false, true>,
                     BuilderConfig<3, true, true>, BuilderConfig<3, false, false>,
                     BuilderConfig<3, true, false>, BuilderConfig<3, false, true>,
                     BuilderConfig<4, true, true>, BuilderConfig<4, false, false>,
                     BuilderConfig<4, true, false>, BuilderConfig<4, false, true>>;

TYPED_TEST_SUITE(CartesianStructInputs, BuilderTypes);

// Test constructor and getModel for all orders and isModelOnNodes values, and
// check resulting model values
TYPED_TEST(CartesianStructInputs, GetModelReturnsValidModel)
{
  // Prepare
  constexpr int order = TypeParam::order;
  constexpr bool isModelOnNodes = TypeParam::isModelOnNodes;
  constexpr bool isElastic = TypeParam::isElastic;


  // Act
  typename TypeParam::Type builder(this->ex, this->lx, this->ey, this->ly,
                                   this->ez, this->lz, isModelOnNodes, isElastic);
  auto model = builder.getModel();

  // Assert
  ASSERT_NE(model, nullptr);
  EXPECT_EQ(model->getOrder(), order);
  EXPECT_EQ(model->isModelOnNodes(), isModelOnNodes);
  EXPECT_EQ(model->isElastic(), isElastic);
  EXPECT_EQ(model->getNumberOfElements(), 10 * 20 * 30);
  EXPECT_EQ(model->getNumberOfNodes(),
            (10 * order + 1) * (20 * order + 1) * (30 * order + 1));
  EXPECT_FLOAT_EQ(model->domainSize(0), 1);
  EXPECT_FLOAT_EQ(model->domainSize(1), 4);
  EXPECT_FLOAT_EQ(model->domainSize(2), 9);
}

// Test multiple calls return different instances
TYPED_TEST(CartesianStructInputs, MultipleCallsReturnDifferentInstances)
{
  // Prepare
  constexpr bool isModelOnNodes = TypeParam::isModelOnNodes;
  constexpr bool isElastic = TypeParam::isElastic;

  // Act
  typename TypeParam::Type builder(this->ex, this->lx, this->ey, this->ly,
                                   this->ez, this->lz, isModelOnNodes, isElastic);

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
  constexpr bool isElastic = TypeParam::isElastic;

  // Act
  typename TypeParam::Type builder(this->ex, this->lx, this->ey, this->ly,
                                   this->ez, this->lz, isModelOnNodes, isElastic);
  ModelBuilderBase<float, int>* base_ptr = &builder;
  auto model = base_ptr->getModel();

  // Assert
  ASSERT_NE(model, nullptr);
}

}  // namespace test
}  // namespace model
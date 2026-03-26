#include <gtest/gtest.h>

#include "data_type.h"
#include "wavefield_view_forward_elastic.h"

namespace gradient
{
namespace test
{

// =============================================================================
// WavefieldViewForwardElastic
// =============================================================================

class WavefieldViewForwardElasticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    size = 45;
    ux_n = allocateVector<VECTOR_REAL_VIEW>(size, "ux_n");
    uy_n = allocateVector<VECTOR_REAL_VIEW>(size, "uy_n");
    uz_n = allocateVector<VECTOR_REAL_VIEW>(size, "uz_n");

    for (int i = 0; i < size; ++i)
    {
      ux_n(i) = static_cast<float>(i) * 1.0f;
      uy_n(i) = static_cast<float>(i) * 2.0f;
      uz_n(i) = static_cast<float>(i) * 3.0f;
    }
  }

  int size;
  VECTOR_REAL_VIEW ux_n;
  VECTOR_REAL_VIEW uy_n;
  VECTOR_REAL_VIEW uz_n;
};

TEST_F(WavefieldViewForwardElasticTest, NumFieldsConstant)
{
  EXPECT_EQ(WavefieldViewForwardElastic::kNumFields, 3);
}

TEST_F(WavefieldViewForwardElasticTest, GetFieldNamesCorrect)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  EXPECT_EQ(view.getFieldName(0), "ux_n");
  EXPECT_EQ(view.getFieldName(1), "uy_n");
  EXPECT_EQ(view.getFieldName(2), "uz_n");
}

TEST_F(WavefieldViewForwardElasticTest, ConstructorStoresFields)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  EXPECT_EQ(view.getField(0).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(1).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(2).extent(0), static_cast<size_t>(size));
}

TEST_F(WavefieldViewForwardElasticTest, GetNumFieldsReturns3)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  EXPECT_EQ(view.getNumFields(), 3);
}

TEST_F(WavefieldViewForwardElasticTest, GetFieldNameByIndex)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  EXPECT_EQ(view.getFieldName(0), "ux_n");
  EXPECT_EQ(view.getFieldName(1), "uy_n");
  EXPECT_EQ(view.getFieldName(2), "uz_n");
}

TEST_F(WavefieldViewForwardElasticTest, GetFieldZeroReturnsUxN)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  auto field = view.getField(0);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 1.0f);
}

TEST_F(WavefieldViewForwardElasticTest, GetFieldOneReturnsUyN)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  auto field = view.getField(1);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 2.0f);
}

TEST_F(WavefieldViewForwardElasticTest, GetFieldTwoReturnsUzN)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  auto field = view.getField(2);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 3.0f);
}

TEST_F(WavefieldViewForwardElasticTest, GetFieldOutOfBoundsFallsBackToUxN)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  auto fallback = view.getField(99);
  EXPECT_EQ(fallback.extent(0), static_cast<size_t>(size));
  EXPECT_FLOAT_EQ(fallback(0), ux_n(0));
}

TEST_F(WavefieldViewForwardElasticTest, ViewsAreShallowCopies)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  ux_n(0) = 100.0f;
  uy_n(0) = 200.0f;
  uz_n(0) = 300.0f;
  EXPECT_FLOAT_EQ(view.getField(0)(0), 100.0f);
  EXPECT_FLOAT_EQ(view.getField(1)(0), 200.0f);
  EXPECT_FLOAT_EQ(view.getField(2)(0), 300.0f);
}

TEST_F(WavefieldViewForwardElasticTest, PrintDoesNotThrow)
{
  WavefieldViewForwardElastic view(ux_n, uy_n, uz_n);
  EXPECT_NO_THROW(view.print());
}

TEST_F(WavefieldViewForwardElasticTest, PolymorphicInterface)
{
  std::unique_ptr<WavefieldView> view =
      std::make_unique<WavefieldViewForwardElastic>(ux_n, uy_n, uz_n);

  EXPECT_EQ(view->getNumFields(), 3);
  EXPECT_EQ(view->getFieldName(1), "uy_n");
  EXPECT_EQ(view->getField(2).extent(0), static_cast<size_t>(size));
}

}  // namespace test
}  // namespace gradient

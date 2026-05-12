#include <gtest/gtest.h>

#include "data_type.h"
#include "wavefield_view_backward_elastic.h"

namespace gradient {
namespace test {

// =============================================================================
// WavefieldViewBackwardElastic
// =============================================================================

class WavefieldViewBackwardElasticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size = 45;
    ux_n = allocateVector<vectorReal>(size, "ux_n");
    uy_n = allocateVector<vectorReal>(size, "uy_n");
    uz_n = allocateVector<vectorReal>(size, "uz_n");
    ux_dt2 = allocateVector<vectorReal>(size, "ux_dt2");
    uy_dt2 = allocateVector<vectorReal>(size, "uy_dt2");
    uz_dt2 = allocateVector<vectorReal>(size, "uz_dt2");

    for (int i = 0; i < size; ++i) {
      ux_n(i) = static_cast<float>(i) * 1.0f;
      uy_n(i) = static_cast<float>(i) * 2.0f;
      uz_n(i) = static_cast<float>(i) * 3.0f;
      ux_dt2(i) = static_cast<float>(i) * 4.0f;
      uy_dt2(i) = static_cast<float>(i) * 5.0f;
      uz_dt2(i) = static_cast<float>(i) * 6.0f;
    }
  }

  int size;
  vectorReal ux_n, uy_n, uz_n;
  vectorReal ux_dt2, uy_dt2, uz_dt2;
};

TEST_F(WavefieldViewBackwardElasticTest, NumFieldsConstant) { EXPECT_EQ(WavefieldViewBackwardElastic::kNumFields, 6); }

TEST_F(WavefieldViewBackwardElasticTest, GetFieldNamesCorrect) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  EXPECT_EQ(view.getFieldName(0), "ux_n");
  EXPECT_EQ(view.getFieldName(1), "uy_n");
  EXPECT_EQ(view.getFieldName(2), "uz_n");
  EXPECT_EQ(view.getFieldName(3), "ux_dt2");
  EXPECT_EQ(view.getFieldName(4), "uy_dt2");
  EXPECT_EQ(view.getFieldName(5), "uz_dt2");
}

TEST_F(WavefieldViewBackwardElasticTest, ConstructorStoresAllFields) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  EXPECT_EQ(view.getField(0).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(1).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(2).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(3).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(4).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(5).extent(0), static_cast<size_t>(size));
}

TEST_F(WavefieldViewBackwardElasticTest, GetNumFieldsReturns6) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  EXPECT_EQ(view.getNumFields(), 6);
}

TEST_F(WavefieldViewBackwardElasticTest, GetFieldNameByIndex) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  EXPECT_EQ(view.getFieldName(0), "ux_n");
  EXPECT_EQ(view.getFieldName(1), "uy_n");
  EXPECT_EQ(view.getFieldName(2), "uz_n");
  EXPECT_EQ(view.getFieldName(3), "ux_dt2");
  EXPECT_EQ(view.getFieldName(4), "uy_dt2");
  EXPECT_EQ(view.getFieldName(5), "uz_dt2");
}

TEST_F(WavefieldViewBackwardElasticTest, GetFieldReturnsCorrectViews) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);

  auto f0 = view.getField(0);
  auto f1 = view.getField(1);
  auto f2 = view.getField(2);
  auto f3 = view.getField(3);
  auto f4 = view.getField(4);
  auto f5 = view.getField(5);

  for (int i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(f0(i), static_cast<float>(i) * 1.0f);
    EXPECT_FLOAT_EQ(f1(i), static_cast<float>(i) * 2.0f);
    EXPECT_FLOAT_EQ(f2(i), static_cast<float>(i) * 3.0f);
    EXPECT_FLOAT_EQ(f3(i), static_cast<float>(i) * 4.0f);
    EXPECT_FLOAT_EQ(f4(i), static_cast<float>(i) * 5.0f);
    EXPECT_FLOAT_EQ(f5(i), static_cast<float>(i) * 6.0f);
  }
}

TEST_F(WavefieldViewBackwardElasticTest, GetFieldOutOfBoundsFallsBackToUxN) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  auto fallback = view.getField(99);
  EXPECT_EQ(fallback.extent(0), static_cast<size_t>(size));
  EXPECT_FLOAT_EQ(fallback(0), ux_n(0));
}

TEST_F(WavefieldViewBackwardElasticTest, ViewsAreShallowCopies) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  ux_dt2(5) = 999.0f;
  EXPECT_FLOAT_EQ(view.getField(3)(5), 999.0f);
}

TEST_F(WavefieldViewBackwardElasticTest, PrintDoesNotThrow) {
  WavefieldViewBackwardElastic view(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);
  EXPECT_NO_THROW(view.print());
}

TEST_F(WavefieldViewBackwardElasticTest, PolymorphicInterface) {
  std::unique_ptr<WavefieldView> view =
      std::make_unique<WavefieldViewBackwardElastic>(ux_n, uy_n, uz_n, ux_dt2, uy_dt2, uz_dt2);

  EXPECT_EQ(view->getNumFields(), 6);
  EXPECT_EQ(view->getFieldName(3), "ux_dt2");
  EXPECT_EQ(view->getField(5).extent(0), static_cast<size_t>(size));
}

}  // namespace test
}  // namespace gradient

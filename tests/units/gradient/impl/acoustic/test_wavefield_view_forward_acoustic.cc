#include <gtest/gtest.h>

#include "data_type.h"
#include "wavefield_view_forward_acoustic.h"

namespace gradient {
namespace test {

class WavefieldViewForwardAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size = 50;
    pn = allocateVector<vectorReal>(size, "pn");

    for (int i = 0; i < size; ++i) pn(i) = static_cast<float>(i) * 0.1f;
  }

  int size;
  vectorReal pn;
};

TEST_F(WavefieldViewForwardAcousticTest, NumFieldsConstant) { EXPECT_EQ(WavefieldViewForwardAcoustic::kNumFields, 1); }

TEST_F(WavefieldViewForwardAcousticTest, GetFieldNameCorrect) {
  WavefieldViewForwardAcoustic view(pn);
  EXPECT_EQ(view.getFieldName(0), "pn");
}

TEST_F(WavefieldViewForwardAcousticTest, ConstructorStoresPn) {
  WavefieldViewForwardAcoustic view(pn);
  EXPECT_EQ(view.getField(0).extent(0), static_cast<size_t>(size));
}

TEST_F(WavefieldViewForwardAcousticTest, GetNumFieldsReturns1) {
  WavefieldViewForwardAcoustic view(pn);
  EXPECT_EQ(view.getNumFields(), 1);
}

TEST_F(WavefieldViewForwardAcousticTest, GetFieldNameReturnsCorrectString) {
  WavefieldViewForwardAcoustic view(pn);
  EXPECT_EQ(view.getFieldName(0), "pn");
}

TEST_F(WavefieldViewForwardAcousticTest, GetFieldZeroReturnsPn) {
  WavefieldViewForwardAcoustic view(pn);
  auto field = view.getField(0);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 0.1f);
}

TEST_F(WavefieldViewForwardAcousticTest, GetFieldOutOfBoundsReturnsPn) {
  // All indices return m_pn (only one field)
  WavefieldViewForwardAcoustic view(pn);
  auto field = view.getField(5);
  EXPECT_EQ(field.extent(0), static_cast<size_t>(size));
}

TEST_F(WavefieldViewForwardAcousticTest, ViewIsShallowCopy) {
  WavefieldViewForwardAcoustic view(pn);
  pn(0) = 123.0f;
  EXPECT_FLOAT_EQ(view.getField(0)(0), 123.0f);
}

TEST_F(WavefieldViewForwardAcousticTest, PrintDoesNotThrow) {
  WavefieldViewForwardAcoustic view(pn);
  EXPECT_NO_THROW(view.print());
}

TEST_F(WavefieldViewForwardAcousticTest, PolymorphicInterface) {
  std::unique_ptr<WavefieldView> view = std::make_unique<WavefieldViewForwardAcoustic>(pn);

  EXPECT_EQ(view->getNumFields(), 1);
  EXPECT_EQ(view->getFieldName(0), "pn");
  EXPECT_EQ(view->getField(0).extent(0), static_cast<size_t>(size));
}

}  // namespace test
}  // namespace gradient

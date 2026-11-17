#include <gtest/gtest.h>

#include "data_type.h"
#include "device_view_access.h"
#include "wavefield_view_backward_acoustic.h"

namespace gradient {
namespace test {

using test_utils::fillFromHost;
using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

class WavefieldViewBackwardAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size = 50;
    qn = allocateVector<vectorReal>(size, "qn");
    qnPrev = allocateVector<vectorReal>(size, "qnPrev");
    qnPrevPrev = allocateVector<vectorReal>(size, "qnPrevPrev");
    dt = 0.001f;

    fillFromHost(qn, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 2.0f;
    });
    fillFromHost(qnPrev, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 3.0f;
    });
    fillFromHost(qnPrevPrev, [&](const auto& h) {
      for (int i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 4.0f;
    });
  }

  int size;
  float dt;
  vectorReal qn;
  vectorReal qnPrev;
  vectorReal qnPrevPrev;
};

TEST_F(WavefieldViewBackwardAcousticTest, NumFieldsConstant) {
  EXPECT_EQ(WavefieldViewBackwardAcoustic::kNumFields, 3);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldNamesCorrect) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getFieldName(0), "qn");
  EXPECT_EQ(view.getFieldName(1), "qnPrev");
  EXPECT_EQ(view.getFieldName(2), "qnPrevPrev");
}

TEST_F(WavefieldViewBackwardAcousticTest, ConstructorStoresFields) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getField(0).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(1).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(2).extent(0), static_cast<size_t>(size));
}

TEST_F(WavefieldViewBackwardAcousticTest, GetNumFieldsReturns3) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getNumFields(), 3);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldNameByIndex) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getFieldName(0), "qn");
  EXPECT_EQ(view.getFieldName(1), "qnPrev");
  EXPECT_EQ(view.getFieldName(2), "qnPrevPrev");
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldZeroReturnsQn) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto field = toHost(view.getField(0));

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 2.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldOneReturnsQnPrev) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto field = toHost(view.getField(1));

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 3.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldTwoReturnsQnPrevPrev) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto field = toHost(view.getField(2));

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i) EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 4.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldOutOfBoundsFallsBackToQn) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto fallback = toHost(view.getField(99));

  EXPECT_EQ(fallback.extent(0), static_cast<size_t>(size));
  EXPECT_FLOAT_EQ(fallback(0), getAt(qn, 0));
}

TEST_F(WavefieldViewBackwardAcousticTest, QnViewIsShallowCopy) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  setAt(qn, 555.0f, 1);
  EXPECT_FLOAT_EQ(getAt(view.getField(0), 1), 555.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, QnPrevViewIsShallowCopy) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  setAt(qnPrev, 666.0f, 3);
  EXPECT_FLOAT_EQ(getAt(view.getField(1), 3), 666.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, QnPrevPrevViewIsShallowCopy) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  setAt(qnPrevPrev, 777.0f, 5);
  EXPECT_FLOAT_EQ(getAt(view.getField(2), 5), 777.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, PrintDoesNotThrow) {
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_NO_THROW(view.print());
}

TEST_F(WavefieldViewBackwardAcousticTest, PolymorphicInterface) {
  std::unique_ptr<WavefieldView> view = std::make_unique<WavefieldViewBackwardAcoustic>(qn, qnPrev, qnPrevPrev);

  EXPECT_EQ(view->getNumFields(), 3);
  EXPECT_EQ(view->getFieldName(0), "qn");
  EXPECT_EQ(view->getField(2).extent(0), static_cast<size_t>(size));
}

}  // namespace test
}  // namespace gradient

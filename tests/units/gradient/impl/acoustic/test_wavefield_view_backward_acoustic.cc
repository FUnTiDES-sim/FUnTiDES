#include <gtest/gtest.h>

#include "data_type.h"
#include "wavefield_view_backward_acoustic.h"

namespace gradient
{
namespace test
{

class WavefieldViewBackwardAcousticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    size = 50;
    qn   = allocateVector<VECTOR_REAL_VIEW>(size, "qn");
    qdt2 = allocateVector<VECTOR_REAL_VIEW>(size, "qdt2");

    for (int i = 0; i < size; ++i)
    {
      qn(i)   = static_cast<float>(i) * 2.0f;
      qdt2(i) = static_cast<float>(i) * 3.0f;
    }
  }

  int size;
  VECTOR_REAL_VIEW qn;
  VECTOR_REAL_VIEW qdt2;
};

TEST_F(WavefieldViewBackwardAcousticTest, NumFieldsConstant)
{
  EXPECT_EQ(WavefieldViewBackwardAcoustic::kNumFields, 2);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldNamesCorrect)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  EXPECT_EQ(view.getFieldName(0), "qn");
  EXPECT_EQ(view.getFieldName(1), "qdt2");
}

TEST_F(WavefieldViewBackwardAcousticTest, ConstructorStoresFields)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  EXPECT_EQ(view.getField(0).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(1).extent(0), static_cast<size_t>(size));
}

TEST_F(WavefieldViewBackwardAcousticTest, GetNumFieldsReturns2)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  EXPECT_EQ(view.getNumFields(), 2);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldNameByIndex)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  EXPECT_EQ(view.getFieldName(0), "qn");
  EXPECT_EQ(view.getFieldName(1), "qdt2");
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldZeroReturnsQn)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  auto field = view.getField(0);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 2.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldOneReturnsQdt2)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  auto field = view.getField(1);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 3.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldOutOfBoundsFallsBackToQn)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  auto fallback = view.getField(99);

  EXPECT_EQ(fallback.extent(0), static_cast<size_t>(size));
  EXPECT_FLOAT_EQ(fallback(0), qn(0));
}

TEST_F(WavefieldViewBackwardAcousticTest, QnViewIsShallowCopy)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  qn(1) = 555.0f;
  EXPECT_FLOAT_EQ(view.getField(0)(1), 555.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, Qdt2ViewIsShallowCopy)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  qdt2(3) = 666.0f;
  EXPECT_FLOAT_EQ(view.getField(1)(3), 666.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, PrintDoesNotThrow)
{
  WavefieldViewBackwardAcoustic view(qn, qdt2);
  EXPECT_NO_THROW(view.print());
}

TEST_F(WavefieldViewBackwardAcousticTest, PolymorphicInterface)
{
  std::unique_ptr<WavefieldView> view =
      std::make_unique<WavefieldViewBackwardAcoustic>(qn, qdt2);

  EXPECT_EQ(view->getNumFields(), 2);
  EXPECT_EQ(view->getFieldName(0), "qn");
  EXPECT_EQ(view->getField(1).extent(0), static_cast<size_t>(size));
}

}  // namespace test
}  // namespace gradient

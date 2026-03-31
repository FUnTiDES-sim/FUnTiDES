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
    qn = allocateVector<VECTOR_REAL_VIEW>(size, "qn");
    qnPrev = allocateVector<VECTOR_REAL_VIEW>(size, "qnPrev");
    qnPrevPrev = allocateVector<VECTOR_REAL_VIEW>(size, "qnPrevPrev");
    dt = 0.001f;

    for (int i = 0; i < size; ++i)
    {
      qn(i) = static_cast<float>(i) * 2.0f;
      qnPrev(i) = static_cast<float>(i) * 3.0f;
      qnPrevPrev(i) = static_cast<float>(i) * 4.0f;
    }
  }

  int size;
  float dt;
  VECTOR_REAL_VIEW qn;
  VECTOR_REAL_VIEW qnPrev;
  VECTOR_REAL_VIEW qnPrevPrev;
};

TEST_F(WavefieldViewBackwardAcousticTest, NumFieldsConstant)
{
  EXPECT_EQ(WavefieldViewBackwardAcoustic::kNumFields, 3);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldNamesCorrect)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getFieldName(0), "qn");
  EXPECT_EQ(view.getFieldName(1), "qnPrev");
  EXPECT_EQ(view.getFieldName(2), "qnPrevPrev");
}

TEST_F(WavefieldViewBackwardAcousticTest, ConstructorStoresFields)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getField(0).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(1).extent(0), static_cast<size_t>(size));
  EXPECT_EQ(view.getField(2).extent(0), static_cast<size_t>(size));
}

TEST_F(WavefieldViewBackwardAcousticTest, GetNumFieldsReturns3)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getNumFields(), 3);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldNameByIndex)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_EQ(view.getFieldName(0), "qn");
  EXPECT_EQ(view.getFieldName(1), "qnPrev");
  EXPECT_EQ(view.getFieldName(2), "qnPrevPrev");
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldZeroReturnsQn)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto field = view.getField(0);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 2.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldOneReturnsQnPrev)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto field = view.getField(1);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 3.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldTwoReturnsQnPrevPrev)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto field = view.getField(2);

  ASSERT_EQ(field.extent(0), static_cast<size_t>(size));
  for (int i = 0; i < size; ++i)
    EXPECT_FLOAT_EQ(field(i), static_cast<float>(i) * 4.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, GetFieldOutOfBoundsFallsBackToQn)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  auto fallback = view.getField(99);

  EXPECT_EQ(fallback.extent(0), static_cast<size_t>(size));
  EXPECT_FLOAT_EQ(fallback(0), qn(0));
}

TEST_F(WavefieldViewBackwardAcousticTest, QnViewIsShallowCopy)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  qn(1) = 555.0f;
  EXPECT_FLOAT_EQ(view.getField(0)(1), 555.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, QnPrevViewIsShallowCopy)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  qnPrev(3) = 666.0f;
  EXPECT_FLOAT_EQ(view.getField(1)(3), 666.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, QnPrevPrevViewIsShallowCopy)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  qnPrevPrev(5) = 777.0f;
  EXPECT_FLOAT_EQ(view.getField(2)(5), 777.0f);
}

TEST_F(WavefieldViewBackwardAcousticTest, PrintDoesNotThrow)
{
  WavefieldViewBackwardAcoustic view(qn, qnPrev, qnPrevPrev);
  EXPECT_NO_THROW(view.print());
}

TEST_F(WavefieldViewBackwardAcousticTest, PolymorphicInterface)
{
  std::unique_ptr<WavefieldView> view =
      std::make_unique<WavefieldViewBackwardAcoustic>(qn, qnPrev, qnPrevPrev);

  EXPECT_EQ(view->getNumFields(), 3);
  EXPECT_EQ(view->getFieldName(0), "qn");
  EXPECT_EQ(view->getField(2).extent(0), static_cast<size_t>(size));
}

}  // namespace test
}  // namespace gradient

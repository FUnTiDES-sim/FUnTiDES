#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "wavefield_acoustoelastic.h"

namespace solver {
namespace fe {
namespace test {

class WavefieldAcoustoElasticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size1 = 100;
    size2 = 200;

    pPrev = allocateVector<vectorReal>(size1, "pPrev");
    pCurr = allocateVector<vectorReal>(size1, "pCurr");
    uxPrev = allocateVector<vectorReal>(size1, "uxPrev");
    uxCurr = allocateVector<vectorReal>(size1, "uxCurr");
    uyPrev = allocateVector<vectorReal>(size1, "uyPrev");
    uyCurr = allocateVector<vectorReal>(size1, "uyCurr");
    uzPrev = allocateVector<vectorReal>(size1, "uzPrev");
    uzCurr = allocateVector<vectorReal>(size1, "uzCurr");

    pPrev2 = allocateVector<vectorReal>(size2, "pPrev2");
    pCurr2 = allocateVector<vectorReal>(size2, "pCurr2");
    uxPrev2 = allocateVector<vectorReal>(size2, "uxPrev2");
    uxCurr2 = allocateVector<vectorReal>(size2, "uxCurr2");
    uyPrev2 = allocateVector<vectorReal>(size2, "uyPrev2");
    uyCurr2 = allocateVector<vectorReal>(size2, "uyCurr2");
    uzPrev2 = allocateVector<vectorReal>(size2, "uzPrev2");
    uzCurr2 = allocateVector<vectorReal>(size2, "uzCurr2");

    for (size_t i = 0; i < size1; ++i) {
      pPrev(i) = i;
      pCurr(i) = i * 2;
      uxPrev(i) = i * 3;
      uxCurr(i) = i * 4;
      uyPrev(i) = i * 5;
      uyCurr(i) = i * 6;
      uzPrev(i) = i * 7;
      uzCurr(i) = i * 8;
    }

    for (size_t i = 0; i < size2; ++i) {
      pPrev2(i) = i * 9;
      pCurr2(i) = i * 10;
      uxPrev2(i) = i * 11;
      uxCurr2(i) = i * 12;
      uyPrev2(i) = i * 13;
      uyCurr2(i) = i * 14;
      uzPrev2(i) = i * 15;
      uzCurr2(i) = i * 16;
    }
  }

  size_t size1;
  size_t size2;
  vectorReal pPrev, pCurr;
  vectorReal uxPrev, uxCurr;
  vectorReal uyPrev, uyCurr;
  vectorReal uzPrev, uzCurr;
  vectorReal pPrev2, pCurr2;
  vectorReal uxPrev2, uxCurr2;
  vectorReal uyPrev2, uyCurr2;
  vectorReal uzPrev2, uzCurr2;
};

TEST_F(WavefieldAcoustoElasticTest, Constructor) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  EXPECT_EQ(wf.m_acoustic.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf.m_acoustic.m_pnGlobalCurr.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uxnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uxnGlobalCurr.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uynGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uynGlobalCurr.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uznGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uznGlobalCurr.extent(0), size1);

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalPrev(i), i);
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalCurr(i), i * 2);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalPrev(i), i * 3);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalCurr(i), i * 4);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uynGlobalPrev(i), i * 5);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uynGlobalCurr(i), i * 6);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalPrev(i), i * 7);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalCurr(i), i * 8);
  }
}

TEST_F(WavefieldAcoustoElasticTest, GetNumFields) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  EXPECT_EQ(wf.getNumFields(), 4);
}

TEST_F(WavefieldAcoustoElasticTest, GetFieldNames) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  const char* const* names = wf.getFieldNames();
  EXPECT_STREQ(names[0], "pressure");
  EXPECT_STREQ(names[1], "ux");
  EXPECT_STREQ(names[2], "uy");
  EXPECT_STREQ(names[3], "uz");
}

TEST_F(WavefieldAcoustoElasticTest, GetCurrentField) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  auto p = wf.getCurrentField(0);
  auto ux = wf.getCurrentField(1);
  auto uy = wf.getCurrentField(2);
  auto uz = wf.getCurrentField(3);

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p(i), i * 2);
    EXPECT_FLOAT_EQ(ux(i), i * 4);
    EXPECT_FLOAT_EQ(uy(i), i * 6);
    EXPECT_FLOAT_EQ(uz(i), i * 8);
  }
}

TEST_F(WavefieldAcoustoElasticTest, GetPreviousField) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  auto p = wf.getPreviousField(0);
  auto ux = wf.getPreviousField(1);
  auto uy = wf.getPreviousField(2);
  auto uz = wf.getPreviousField(3);

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p(i), i);
    EXPECT_FLOAT_EQ(ux(i), i * 3);
    EXPECT_FLOAT_EQ(uy(i), i * 5);
    EXPECT_FLOAT_EQ(uz(i), i * 7);
  }
}

TEST_F(WavefieldAcoustoElasticTest, Swap) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  wf.swap();

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalCurr(i), i);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalPrev(i), i * 4);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalCurr(i), i * 3);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uynGlobalPrev(i), i * 6);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uynGlobalCurr(i), i * 5);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalPrev(i), i * 8);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalCurr(i), i * 7);
  }
}

TEST_F(WavefieldAcoustoElasticTest, SwapTwiceRestoresState) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  wf.swap();
  wf.swap();

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalPrev(i), i);
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalCurr(i), i * 2);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalPrev(i), i * 3);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalCurr(i), i * 4);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uynGlobalPrev(i), i * 5);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uynGlobalCurr(i), i * 6);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalPrev(i), i * 7);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalCurr(i), i * 8);
  }
}

TEST_F(WavefieldAcoustoElasticTest, SwapWithThreeBuffersAcoustic) {
  auto pPrevPrev = allocateVector<vectorReal>(size1, "pPrevPrev");
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  for (size_t i = 0; i < size1; ++i) {
    pPrevPrev(i) = 10.0f;
    uxPrevPrev(i) = 100.0f;
    uyPrevPrev(i) = 200.0f;
    uzPrevPrev(i) = 300.0f;
  }

  WavefieldAcoustoElastic wf(pPrevPrev, pPrev, pCurr,
                             uxPrevPrev, uxPrev, uxCurr,
                             uyPrevPrev, uyPrev, uyCurr,
                             uzPrevPrev, uzPrev, uzCurr);
  EXPECT_TRUE(wf.hasPrevPrev());

  wf.swap();

  // After 3-way rotation:
  // curr ← old prevPrev, prev ← old curr, prevPrev ← old prev
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalCurr(i), 10.0f);
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(wf.getPrevPrevField(0)(i), static_cast<float>(i));

    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalCurr(i), 100.0f);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalPrev(i), i * 4);
    EXPECT_FLOAT_EQ(wf.getPrevPrevField(1)(i), i * 3);
  }
}

TEST_F(WavefieldAcoustoElasticTest, SwapWithThreeBuffersElastic) {
  auto pPrevPrev = allocateVector<vectorReal>(size1, "pPrevPrev");
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  for (size_t i = 0; i < size1; ++i) {
    pPrevPrev(i) = 10.0f;
    uxPrevPrev(i) = 20.0f;
    uyPrevPrev(i) = 30.0f;
    uzPrevPrev(i) = 40.0f;
  }

  WavefieldAcoustoElastic wf(pPrevPrev, pPrev, pCurr,
                             uxPrevPrev, uxPrev, uxCurr,
                             uyPrevPrev, uyPrev, uyCurr,
                             uzPrevPrev, uzPrev, uzCurr);
  EXPECT_TRUE(wf.hasPrevPrev());

  wf.swap();

  // After 3-way rotation, verify all elastic components rotated correctly
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalCurr(i), 20.0f);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalPrev(i), i * 4);
    EXPECT_FLOAT_EQ(wf.getPrevPrevField(1)(i), i * 3);

    EXPECT_FLOAT_EQ(wf.m_elastic.m_uynGlobalCurr(i), 30.0f);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalCurr(i), 40.0f);
  }
  // Acoustic also rotated
  EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalCurr(1), 10.0f);
  EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalPrev(1), 2.0f);
}

TEST_F(WavefieldAcoustoElasticTest, CopyConstructor) {
  WavefieldAcoustoElastic original(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  WavefieldAcoustoElastic copy(original);

  EXPECT_EQ(copy.m_acoustic.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(copy.m_elastic.m_uxnGlobalPrev.extent(0), size1);

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy.m_acoustic.m_pnGlobalPrev(i), original.m_acoustic.m_pnGlobalPrev(i));
    EXPECT_FLOAT_EQ(copy.m_elastic.m_uxnGlobalCurr(i), original.m_elastic.m_uxnGlobalCurr(i));
  }

  // Kokkos views are shallow: modifying original affects copy
  original.m_acoustic.m_pnGlobalPrev(0) = 999.0f;
  EXPECT_FLOAT_EQ(copy.m_acoustic.m_pnGlobalPrev(0), 999.0f);
  original.m_elastic.m_uxnGlobalCurr(0) = 888.0f;
  EXPECT_FLOAT_EQ(copy.m_elastic.m_uxnGlobalCurr(0), 888.0f);
  original.m_elastic.m_uynGlobalPrev(0) = 777.0f;
  EXPECT_FLOAT_EQ(copy.m_elastic.m_uynGlobalPrev(0), 777.0f);
}

TEST_F(WavefieldAcoustoElasticTest, CopyAssignmentOperator) {
  WavefieldAcoustoElastic wf1(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  WavefieldAcoustoElastic wf2(pPrev2, pCurr2, uxPrev2, uxCurr2, uyPrev2, uyCurr2, uzPrev2, uzCurr2);

  EXPECT_EQ(wf2.m_acoustic.m_pnGlobalPrev.extent(0), size2);

  wf2 = wf1;

  EXPECT_EQ(wf2.m_acoustic.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf2.m_elastic.m_uxnGlobalPrev.extent(0), size1);

  // Shallow copy after assignment
  wf1.m_acoustic.m_pnGlobalPrev(0) = 777.0f;
  EXPECT_FLOAT_EQ(wf2.m_acoustic.m_pnGlobalPrev(0), 777.0f);
  wf1.m_elastic.m_uznGlobalCurr(0) = 666.0f;
  EXPECT_FLOAT_EQ(wf2.m_elastic.m_uznGlobalCurr(0), 666.0f);
}

TEST_F(WavefieldAcoustoElasticTest, CopyAssignmentSelfAssignment) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  wf = wf;

  EXPECT_EQ(wf.m_acoustic.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uxnGlobalPrev.extent(0), size1);

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalPrev(i), i);
    EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalCurr(i), i * 2);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uxnGlobalPrev(i), i * 3);
    EXPECT_FLOAT_EQ(wf.m_elastic.m_uznGlobalCurr(i), i * 8);
  }
}

TEST_F(WavefieldAcoustoElasticTest, EmptyFields) {
  auto empty = allocateVector<vectorReal>(0, "empty");
  WavefieldAcoustoElastic wf(empty, empty, empty, empty, empty, empty, empty, empty);

  EXPECT_EQ(wf.m_acoustic.m_pnGlobalPrev.extent(0), 0);
  EXPECT_EQ(wf.m_acoustic.m_pnGlobalCurr.extent(0), 0);
  EXPECT_EQ(wf.m_elastic.m_uxnGlobalPrev.extent(0), 0);
  EXPECT_EQ(wf.m_elastic.m_uznGlobalCurr.extent(0), 0);

  wf.swap();
  EXPECT_EQ(wf.m_acoustic.m_pnGlobalPrev.extent(0), 0);
  EXPECT_EQ(wf.m_elastic.m_uxnGlobalPrev.extent(0), 0);
}

TEST_F(WavefieldAcoustoElasticTest, SwapThreeTimesRestoresStateWithThreeBuffers) {
  auto pPrevPrev = allocateVector<vectorReal>(size1, "pPrevPrev");
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  for (size_t i = 0; i < size1; ++i) {
    pPrevPrev(i) = 10.0f;
    uxPrevPrev(i) = 20.0f;
    uyPrevPrev(i) = 30.0f;
    uzPrevPrev(i) = 40.0f;
  }

  float initialPPrev0 = pPrev(0);
  float initialPCurr0 = pCurr(0);
  float initialPPrevPrev0 = pPrevPrev(0);

  WavefieldAcoustoElastic wf(pPrevPrev, pPrev, pCurr,
                             uxPrevPrev, uxPrev, uxCurr,
                             uyPrevPrev, uyPrev, uyCurr,
                             uzPrevPrev, uzPrev, uzCurr);

  wf.swap();
  wf.swap();
  wf.swap();

  EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalPrev(0), initialPPrev0);
  EXPECT_FLOAT_EQ(wf.m_acoustic.m_pnGlobalCurr(0), initialPCurr0);
  EXPECT_FLOAT_EQ(wf.getPrevPrevField(0)(0), initialPPrevPrev0);
}

}  // namespace test
}  // namespace fe
}  // namespace solver

#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "device_view_access.h"
#include "wavefield_acoustoelastic.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromHost;
using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

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

    fillFromHost(pPrev, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i;
    });
    fillFromHost(pCurr, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 2;
    });
    fillFromHost(uxPrev, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 3;
    });
    fillFromHost(uxCurr, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 4;
    });
    fillFromHost(uyPrev, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 5;
    });
    fillFromHost(uyCurr, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 6;
    });
    fillFromHost(uzPrev, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 7;
    });
    fillFromHost(uzCurr, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 8;
    });

    fillFromHost(pPrev2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 9;
    });
    fillFromHost(pCurr2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 10;
    });
    fillFromHost(uxPrev2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 11;
    });
    fillFromHost(uxCurr2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 12;
    });
    fillFromHost(uyPrev2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 13;
    });
    fillFromHost(uyCurr2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 14;
    });
    fillFromHost(uzPrev2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 15;
    });
    fillFromHost(uzCurr2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 16;
    });
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

  auto p_prev = toHost(wf.m_acoustic.m_pnGlobalPrev);
  auto p_curr = toHost(wf.m_acoustic.m_pnGlobalCurr);
  auto ux_prev = toHost(wf.m_elastic.m_uxnGlobalPrev);
  auto ux_curr = toHost(wf.m_elastic.m_uxnGlobalCurr);
  auto uy_prev = toHost(wf.m_elastic.m_uynGlobalPrev);
  auto uy_curr = toHost(wf.m_elastic.m_uynGlobalCurr);
  auto uz_prev = toHost(wf.m_elastic.m_uznGlobalPrev);
  auto uz_curr = toHost(wf.m_elastic.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p_prev(i), i);
    EXPECT_FLOAT_EQ(p_curr(i), i * 2);
    EXPECT_FLOAT_EQ(ux_prev(i), i * 3);
    EXPECT_FLOAT_EQ(ux_curr(i), i * 4);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 5);
    EXPECT_FLOAT_EQ(uy_curr(i), i * 6);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 7);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 8);
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

  auto p = toHost(wf.getCurrentField(0));
  auto ux = toHost(wf.getCurrentField(1));
  auto uy = toHost(wf.getCurrentField(2));
  auto uz = toHost(wf.getCurrentField(3));

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p(i), i * 2);
    EXPECT_FLOAT_EQ(ux(i), i * 4);
    EXPECT_FLOAT_EQ(uy(i), i * 6);
    EXPECT_FLOAT_EQ(uz(i), i * 8);
  }
}

TEST_F(WavefieldAcoustoElasticTest, GetPreviousField) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  auto p = toHost(wf.getPreviousField(0));
  auto ux = toHost(wf.getPreviousField(1));
  auto uy = toHost(wf.getPreviousField(2));
  auto uz = toHost(wf.getPreviousField(3));

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

  auto p_prev = toHost(wf.m_acoustic.m_pnGlobalPrev);
  auto p_curr = toHost(wf.m_acoustic.m_pnGlobalCurr);
  auto ux_prev = toHost(wf.m_elastic.m_uxnGlobalPrev);
  auto ux_curr = toHost(wf.m_elastic.m_uxnGlobalCurr);
  auto uy_prev = toHost(wf.m_elastic.m_uynGlobalPrev);
  auto uy_curr = toHost(wf.m_elastic.m_uynGlobalCurr);
  auto uz_prev = toHost(wf.m_elastic.m_uznGlobalPrev);
  auto uz_curr = toHost(wf.m_elastic.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p_prev(i), i * 2);
    EXPECT_FLOAT_EQ(p_curr(i), i);
    EXPECT_FLOAT_EQ(ux_prev(i), i * 4);
    EXPECT_FLOAT_EQ(ux_curr(i), i * 3);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 6);
    EXPECT_FLOAT_EQ(uy_curr(i), i * 5);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 8);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 7);
  }
}

TEST_F(WavefieldAcoustoElasticTest, SwapTwiceRestoresState) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  wf.swap();
  wf.swap();

  auto p_prev = toHost(wf.m_acoustic.m_pnGlobalPrev);
  auto p_curr = toHost(wf.m_acoustic.m_pnGlobalCurr);
  auto ux_prev = toHost(wf.m_elastic.m_uxnGlobalPrev);
  auto ux_curr = toHost(wf.m_elastic.m_uxnGlobalCurr);
  auto uy_prev = toHost(wf.m_elastic.m_uynGlobalPrev);
  auto uy_curr = toHost(wf.m_elastic.m_uynGlobalCurr);
  auto uz_prev = toHost(wf.m_elastic.m_uznGlobalPrev);
  auto uz_curr = toHost(wf.m_elastic.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p_prev(i), i);
    EXPECT_FLOAT_EQ(p_curr(i), i * 2);
    EXPECT_FLOAT_EQ(ux_prev(i), i * 3);
    EXPECT_FLOAT_EQ(ux_curr(i), i * 4);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 5);
    EXPECT_FLOAT_EQ(uy_curr(i), i * 6);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 7);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 8);
  }
}

TEST_F(WavefieldAcoustoElasticTest, SwapWithThreeBuffersAcoustic) {
  auto pPrevPrev = allocateVector<vectorReal>(size1, "pPrevPrev");
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  Kokkos::deep_copy(pPrevPrev, 10.0f);
  Kokkos::deep_copy(uxPrevPrev, 100.0f);
  Kokkos::deep_copy(uyPrevPrev, 200.0f);
  Kokkos::deep_copy(uzPrevPrev, 300.0f);

  WavefieldAcoustoElastic wf(pPrevPrev, pPrev, pCurr, uxPrevPrev, uxPrev, uxCurr, uyPrevPrev, uyPrev, uyCurr,
                             uzPrevPrev, uzPrev, uzCurr);
  EXPECT_TRUE(wf.hasPrevPrev());

  wf.swap();

  // After 3-way rotation:
  // curr ← old prevPrev, prev ← old curr, prevPrev ← old prev
  auto p_curr = toHost(wf.m_acoustic.m_pnGlobalCurr);
  auto p_prev = toHost(wf.m_acoustic.m_pnGlobalPrev);
  auto p_prevprev = toHost(wf.getPrevPrevField(0));
  auto ux_curr = toHost(wf.m_elastic.m_uxnGlobalCurr);
  auto ux_prev = toHost(wf.m_elastic.m_uxnGlobalPrev);
  auto ux_prevprev = toHost(wf.getPrevPrevField(1));
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p_curr(i), 10.0f);
    EXPECT_FLOAT_EQ(p_prev(i), i * 2);
    EXPECT_FLOAT_EQ(p_prevprev(i), static_cast<float>(i));

    EXPECT_FLOAT_EQ(ux_curr(i), 100.0f);
    EXPECT_FLOAT_EQ(ux_prev(i), i * 4);
    EXPECT_FLOAT_EQ(ux_prevprev(i), i * 3);
  }
}

TEST_F(WavefieldAcoustoElasticTest, SwapWithThreeBuffersElastic) {
  auto pPrevPrev = allocateVector<vectorReal>(size1, "pPrevPrev");
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  Kokkos::deep_copy(pPrevPrev, 10.0f);
  Kokkos::deep_copy(uxPrevPrev, 20.0f);
  Kokkos::deep_copy(uyPrevPrev, 30.0f);
  Kokkos::deep_copy(uzPrevPrev, 40.0f);

  WavefieldAcoustoElastic wf(pPrevPrev, pPrev, pCurr, uxPrevPrev, uxPrev, uxCurr, uyPrevPrev, uyPrev, uyCurr,
                             uzPrevPrev, uzPrev, uzCurr);
  EXPECT_TRUE(wf.hasPrevPrev());

  wf.swap();

  // After 3-way rotation, verify all elastic components rotated correctly
  auto ux_curr = toHost(wf.m_elastic.m_uxnGlobalCurr);
  auto ux_prev = toHost(wf.m_elastic.m_uxnGlobalPrev);
  auto ux_prevprev = toHost(wf.getPrevPrevField(1));
  auto uy_curr = toHost(wf.m_elastic.m_uynGlobalCurr);
  auto uz_curr = toHost(wf.m_elastic.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_curr(i), 20.0f);
    EXPECT_FLOAT_EQ(ux_prev(i), i * 4);
    EXPECT_FLOAT_EQ(ux_prevprev(i), i * 3);

    EXPECT_FLOAT_EQ(uy_curr(i), 30.0f);
    EXPECT_FLOAT_EQ(uz_curr(i), 40.0f);
  }
  // Acoustic also rotated
  EXPECT_FLOAT_EQ(getAt(wf.m_acoustic.m_pnGlobalCurr, 1), 10.0f);
  EXPECT_FLOAT_EQ(getAt(wf.m_acoustic.m_pnGlobalPrev, 1), 2.0f);
}

TEST_F(WavefieldAcoustoElasticTest, CopyConstructor) {
  WavefieldAcoustoElastic original(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  WavefieldAcoustoElastic copy(original);

  EXPECT_EQ(copy.m_acoustic.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(copy.m_elastic.m_uxnGlobalPrev.extent(0), size1);

  auto copy_p_prev = toHost(copy.m_acoustic.m_pnGlobalPrev);
  auto original_p_prev = toHost(original.m_acoustic.m_pnGlobalPrev);
  auto copy_ux_curr = toHost(copy.m_elastic.m_uxnGlobalCurr);
  auto original_ux_curr = toHost(original.m_elastic.m_uxnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy_p_prev(i), original_p_prev(i));
    EXPECT_FLOAT_EQ(copy_ux_curr(i), original_ux_curr(i));
  }

  // Kokkos views are shallow: modifying original affects copy
  setAt(original.m_acoustic.m_pnGlobalPrev, 999.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_acoustic.m_pnGlobalPrev, 0), 999.0f);
  setAt(original.m_elastic.m_uxnGlobalCurr, 888.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_elastic.m_uxnGlobalCurr, 0), 888.0f);
  setAt(original.m_elastic.m_uynGlobalPrev, 777.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_elastic.m_uynGlobalPrev, 0), 777.0f);
}

TEST_F(WavefieldAcoustoElasticTest, CopyAssignmentOperator) {
  WavefieldAcoustoElastic wf1(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  WavefieldAcoustoElastic wf2(pPrev2, pCurr2, uxPrev2, uxCurr2, uyPrev2, uyCurr2, uzPrev2, uzCurr2);

  EXPECT_EQ(wf2.m_acoustic.m_pnGlobalPrev.extent(0), size2);

  wf2 = wf1;

  EXPECT_EQ(wf2.m_acoustic.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf2.m_elastic.m_uxnGlobalPrev.extent(0), size1);

  // Shallow copy after assignment
  setAt(wf1.m_acoustic.m_pnGlobalPrev, 777.0f, 0);
  EXPECT_FLOAT_EQ(getAt(wf2.m_acoustic.m_pnGlobalPrev, 0), 777.0f);
  setAt(wf1.m_elastic.m_uznGlobalCurr, 666.0f, 0);
  EXPECT_FLOAT_EQ(getAt(wf2.m_elastic.m_uznGlobalCurr, 0), 666.0f);
}

TEST_F(WavefieldAcoustoElasticTest, CopyAssignmentSelfAssignment) {
  WavefieldAcoustoElastic wf(pPrev, pCurr, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);

  wf = wf;

  EXPECT_EQ(wf.m_acoustic.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wf.m_elastic.m_uxnGlobalPrev.extent(0), size1);

  auto p_prev = toHost(wf.m_acoustic.m_pnGlobalPrev);
  auto p_curr = toHost(wf.m_acoustic.m_pnGlobalCurr);
  auto ux_prev = toHost(wf.m_elastic.m_uxnGlobalPrev);
  auto uz_curr = toHost(wf.m_elastic.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(p_prev(i), i);
    EXPECT_FLOAT_EQ(p_curr(i), i * 2);
    EXPECT_FLOAT_EQ(ux_prev(i), i * 3);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 8);
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
  Kokkos::deep_copy(pPrevPrev, 10.0f);
  Kokkos::deep_copy(uxPrevPrev, 20.0f);
  Kokkos::deep_copy(uyPrevPrev, 30.0f);
  Kokkos::deep_copy(uzPrevPrev, 40.0f);

  float initialPPrev0 = getAt(pPrev, 0);
  float initialPCurr0 = getAt(pCurr, 0);
  float initialPPrevPrev0 = getAt(pPrevPrev, 0);

  WavefieldAcoustoElastic wf(pPrevPrev, pPrev, pCurr, uxPrevPrev, uxPrev, uxCurr, uyPrevPrev, uyPrev, uyCurr,
                             uzPrevPrev, uzPrev, uzCurr);

  wf.swap();
  wf.swap();
  wf.swap();

  EXPECT_FLOAT_EQ(getAt(wf.m_acoustic.m_pnGlobalPrev, 0), initialPPrev0);
  EXPECT_FLOAT_EQ(getAt(wf.m_acoustic.m_pnGlobalCurr, 0), initialPCurr0);
  EXPECT_FLOAT_EQ(getAt(wf.getPrevPrevField(0), 0), initialPPrevPrev0);
}

}  // namespace test
}  // namespace fe
}  // namespace solver

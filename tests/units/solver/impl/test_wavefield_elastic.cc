#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "device_view_access.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromHost;
using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

class WavefieldElasticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create test data with specific sizes
    size1 = 100;
    size2 = 200;

    uxPrevField = allocateVector<vectorReal>(size1, "uxPrevField");
    uxCurrField = allocateVector<vectorReal>(size1, "uxCurrField");
    uyPrevField = allocateVector<vectorReal>(size1, "uyPrevField");
    uyCurrField = allocateVector<vectorReal>(size1, "uyCurrField");
    uzPrevField = allocateVector<vectorReal>(size1, "uzPrevField");
    uzCurrField = allocateVector<vectorReal>(size1, "uzCurrField");

    uxPrevField2 = allocateVector<vectorReal>(size2, "uxPrevField2");
    uxCurrField2 = allocateVector<vectorReal>(size2, "uxCurrField2");
    uyPrevField2 = allocateVector<vectorReal>(size2, "uyPrevField2");
    uyCurrField2 = allocateVector<vectorReal>(size2, "uyCurrField2");
    uzPrevField2 = allocateVector<vectorReal>(size2, "uzPrevField2");
    uzCurrField2 = allocateVector<vectorReal>(size2, "uzCurrField2");

    // Initialize with test values
    fillFromHost(uxPrevField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i;
    });
    fillFromHost(uxCurrField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 2;
    });
    fillFromHost(uyPrevField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 3;
    });
    fillFromHost(uyCurrField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 4;
    });
    fillFromHost(uzPrevField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 5;
    });
    fillFromHost(uzCurrField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 6;
    });

    fillFromHost(uxPrevField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 7;
    });
    fillFromHost(uxCurrField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 8;
    });
    fillFromHost(uyPrevField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 9;
    });
    fillFromHost(uyCurrField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 10;
    });
    fillFromHost(uzPrevField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 11;
    });
    fillFromHost(uzCurrField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 12;
    });
  }

  size_t size1;
  size_t size2;
  vectorReal uxPrevField, uxCurrField;
  vectorReal uyPrevField, uyCurrField;
  vectorReal uzPrevField, uzCurrField;
  vectorReal uxPrevField2, uxCurrField2;
  vectorReal uyPrevField2, uyCurrField2;
  vectorReal uzPrevField2, uzCurrField2;
};

TEST_F(WavefieldElasticTest, Constructor) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  EXPECT_EQ(wavefield.m_uxnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_uxnGlobalCurr.extent(0), size1);
  EXPECT_EQ(wavefield.m_uynGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_uynGlobalCurr.extent(0), size1);
  EXPECT_EQ(wavefield.m_uznGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_uznGlobalCurr.extent(0), size1);

  // Verify data is correctly stored
  auto ux_prev = toHost(wavefield.m_uxnGlobalPrev);
  auto ux_curr = toHost(wavefield.m_uxnGlobalCurr);
  auto uy_prev = toHost(wavefield.m_uynGlobalPrev);
  auto uy_curr = toHost(wavefield.m_uynGlobalCurr);
  auto uz_prev = toHost(wavefield.m_uznGlobalPrev);
  auto uz_curr = toHost(wavefield.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_prev(i), i);
    EXPECT_FLOAT_EQ(ux_curr(i), i * 2);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 3);
    EXPECT_FLOAT_EQ(uy_curr(i), i * 4);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 5);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 6);
  }
}

TEST_F(WavefieldElasticTest, CopyConstructor) {
  WavefieldElastic original(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);
  WavefieldElastic copy(original);

  // Check that copy has the same extent
  EXPECT_EQ(copy.m_uxnGlobalPrev.extent(0), original.m_uxnGlobalPrev.extent(0));
  EXPECT_EQ(copy.m_uxnGlobalCurr.extent(0), original.m_uxnGlobalCurr.extent(0));
  EXPECT_EQ(copy.m_uynGlobalPrev.extent(0), original.m_uynGlobalPrev.extent(0));
  EXPECT_EQ(copy.m_uynGlobalCurr.extent(0), original.m_uynGlobalCurr.extent(0));
  EXPECT_EQ(copy.m_uznGlobalPrev.extent(0), original.m_uznGlobalPrev.extent(0));
  EXPECT_EQ(copy.m_uznGlobalCurr.extent(0), original.m_uznGlobalCurr.extent(0));

  // Check that copy has the same data
  auto copy_ux_prev = toHost(copy.m_uxnGlobalPrev), orig_ux_prev = toHost(original.m_uxnGlobalPrev);
  auto copy_ux_curr = toHost(copy.m_uxnGlobalCurr), orig_ux_curr = toHost(original.m_uxnGlobalCurr);
  auto copy_uy_prev = toHost(copy.m_uynGlobalPrev), orig_uy_prev = toHost(original.m_uynGlobalPrev);
  auto copy_uy_curr = toHost(copy.m_uynGlobalCurr), orig_uy_curr = toHost(original.m_uynGlobalCurr);
  auto copy_uz_prev = toHost(copy.m_uznGlobalPrev), orig_uz_prev = toHost(original.m_uznGlobalPrev);
  auto copy_uz_curr = toHost(copy.m_uznGlobalCurr), orig_uz_curr = toHost(original.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy_ux_prev(i), orig_ux_prev(i));
    EXPECT_FLOAT_EQ(copy_ux_curr(i), orig_ux_curr(i));
    EXPECT_FLOAT_EQ(copy_uy_prev(i), orig_uy_prev(i));
    EXPECT_FLOAT_EQ(copy_uy_curr(i), orig_uy_curr(i));
    EXPECT_FLOAT_EQ(copy_uz_prev(i), orig_uz_prev(i));
    EXPECT_FLOAT_EQ(copy_uz_curr(i), orig_uz_curr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(original.m_uxnGlobalPrev, 999.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_uxnGlobalPrev, 0), 999.0f);
  setAt(original.m_uxnGlobalCurr, 888.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_uxnGlobalCurr, 0), 888.0f);
  setAt(original.m_uynGlobalPrev, 777.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_uynGlobalPrev, 0), 777.0f);
  setAt(original.m_uynGlobalCurr, 666.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_uynGlobalCurr, 0), 666.0f);
  setAt(original.m_uznGlobalPrev, 555.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_uznGlobalPrev, 0), 555.0f);
  setAt(original.m_uznGlobalCurr, 444.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_uznGlobalCurr, 0), 444.0f);
}

TEST_F(WavefieldElasticTest, CopyAssignmentOperator) {
  WavefieldElastic wavefield1(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);
  WavefieldElastic wavefield2(uxPrevField2, uxCurrField2, uyPrevField2, uyCurrField2, uzPrevField2, uzCurrField2);

  // Verify initial state
  EXPECT_EQ(wavefield2.m_uxnGlobalPrev.extent(0), size2);
  EXPECT_EQ(wavefield2.m_uxnGlobalCurr.extent(0), size2);

  // Perform assignment
  wavefield2 = wavefield1;

  // Check that wavefield2 now has the same extent as wavefield1
  EXPECT_EQ(wavefield2.m_uxnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield2.m_uxnGlobalCurr.extent(0), size1);
  EXPECT_EQ(wavefield2.m_uynGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield2.m_uynGlobalCurr.extent(0), size1);
  EXPECT_EQ(wavefield2.m_uznGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield2.m_uznGlobalCurr.extent(0), size1);

  // Check that data matches
  auto w2_ux_prev = toHost(wavefield2.m_uxnGlobalPrev), w1_ux_prev = toHost(wavefield1.m_uxnGlobalPrev);
  auto w2_ux_curr = toHost(wavefield2.m_uxnGlobalCurr), w1_ux_curr = toHost(wavefield1.m_uxnGlobalCurr);
  auto w2_uy_prev = toHost(wavefield2.m_uynGlobalPrev), w1_uy_prev = toHost(wavefield1.m_uynGlobalPrev);
  auto w2_uy_curr = toHost(wavefield2.m_uynGlobalCurr), w1_uy_curr = toHost(wavefield1.m_uynGlobalCurr);
  auto w2_uz_prev = toHost(wavefield2.m_uznGlobalPrev), w1_uz_prev = toHost(wavefield1.m_uznGlobalPrev);
  auto w2_uz_curr = toHost(wavefield2.m_uznGlobalCurr), w1_uz_curr = toHost(wavefield1.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(w2_ux_prev(i), w1_ux_prev(i));
    EXPECT_FLOAT_EQ(w2_ux_curr(i), w1_ux_curr(i));
    EXPECT_FLOAT_EQ(w2_uy_prev(i), w1_uy_prev(i));
    EXPECT_FLOAT_EQ(w2_uy_curr(i), w1_uy_curr(i));
    EXPECT_FLOAT_EQ(w2_uz_prev(i), w1_uz_prev(i));
    EXPECT_FLOAT_EQ(w2_uz_curr(i), w1_uz_curr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(wavefield1.m_uxnGlobalPrev, 777.0f, 0);
  EXPECT_FLOAT_EQ(getAt(wavefield2.m_uxnGlobalPrev, 0), 777.0f);
  setAt(wavefield1.m_uxnGlobalCurr, 666.0f, 0);
  EXPECT_FLOAT_EQ(getAt(wavefield2.m_uxnGlobalCurr, 0), 666.0f);
}

TEST_F(WavefieldElasticTest, CopyAssignmentSelfAssignment) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Self-assignment should not cause issues
  wavefield = wavefield;

  // Verify data is unchanged
  EXPECT_EQ(wavefield.m_uxnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_uxnGlobalCurr.extent(0), size1);
  EXPECT_EQ(wavefield.m_uynGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_uynGlobalCurr.extent(0), size1);
  EXPECT_EQ(wavefield.m_uznGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_uznGlobalCurr.extent(0), size1);

  auto ux_prev = toHost(wavefield.m_uxnGlobalPrev);
  auto ux_curr = toHost(wavefield.m_uxnGlobalCurr);
  auto uy_prev = toHost(wavefield.m_uynGlobalPrev);
  auto uy_curr = toHost(wavefield.m_uynGlobalCurr);
  auto uz_prev = toHost(wavefield.m_uznGlobalPrev);
  auto uz_curr = toHost(wavefield.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_prev(i), i);
    EXPECT_FLOAT_EQ(ux_curr(i), i * 2);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 3);
    EXPECT_FLOAT_EQ(uy_curr(i), i * 4);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 5);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 6);
  }
}

TEST_F(WavefieldElasticTest, GetCurrentField) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  auto currentX = wavefield.getCurrentField(0);
  auto currentY = wavefield.getCurrentField(1);
  auto currentZ = wavefield.getCurrentField(2);

  EXPECT_EQ(currentX.extent(0), size1);
  EXPECT_EQ(currentY.extent(0), size1);
  EXPECT_EQ(currentZ.extent(0), size1);

  auto h_x = toHost(currentX), h_y = toHost(currentY), h_z = toHost(currentZ);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(h_x(i), i * 2);
    EXPECT_FLOAT_EQ(h_y(i), i * 4);
    EXPECT_FLOAT_EQ(h_z(i), i * 6);
  }
}

TEST_F(WavefieldElasticTest, GetPreviousField) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  auto previousX = wavefield.getPreviousField(0);
  auto previousY = wavefield.getPreviousField(1);
  auto previousZ = wavefield.getPreviousField(2);

  EXPECT_EQ(previousX.extent(0), size1);
  EXPECT_EQ(previousY.extent(0), size1);
  EXPECT_EQ(previousZ.extent(0), size1);

  auto h_x = toHost(previousX), h_y = toHost(previousY), h_z = toHost(previousZ);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(h_x(i), i);
    EXPECT_FLOAT_EQ(h_y(i), i * 3);
    EXPECT_FLOAT_EQ(h_z(i), i * 5);
  }
}

TEST_F(WavefieldElasticTest, Swap) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Store original values
  float originalUxPrev0 = getAt(wavefield.m_uxnGlobalPrev, 0);
  float originalUxCurr0 = getAt(wavefield.m_uxnGlobalCurr, 0);
  float originalUyPrev0 = getAt(wavefield.m_uynGlobalPrev, 0);
  float originalUyCurr0 = getAt(wavefield.m_uynGlobalCurr, 0);
  float originalUzPrev0 = getAt(wavefield.m_uznGlobalPrev, 0);
  float originalUzCurr0 = getAt(wavefield.m_uznGlobalCurr, 0);

  // Perform swap
  wavefield.swap();

  // Verify that prev and curr have been swapped
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uxnGlobalPrev, 0), originalUxCurr0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uxnGlobalCurr, 0), originalUxPrev0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uynGlobalPrev, 0), originalUyCurr0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uynGlobalCurr, 0), originalUyPrev0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uznGlobalPrev, 0), originalUzCurr0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uznGlobalCurr, 0), originalUzPrev0);

  // Verify all elements were swapped
  auto ux_prev = toHost(wavefield.m_uxnGlobalPrev);
  auto ux_curr = toHost(wavefield.m_uxnGlobalCurr);
  auto uy_prev = toHost(wavefield.m_uynGlobalPrev);
  auto uy_curr = toHost(wavefield.m_uynGlobalCurr);
  auto uz_prev = toHost(wavefield.m_uznGlobalPrev);
  auto uz_curr = toHost(wavefield.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_prev(i), i * 2);
    EXPECT_FLOAT_EQ(ux_curr(i), i);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 4);
    EXPECT_FLOAT_EQ(uy_curr(i), i * 3);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 6);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 5);
  }
}

TEST_F(WavefieldElasticTest, SwapTwice) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Store original values
  auto original_ux_prev = toHost(wavefield.m_uxnGlobalPrev);
  auto original_ux_curr = toHost(wavefield.m_uxnGlobalCurr);
  auto original_uy_prev = toHost(wavefield.m_uynGlobalPrev);
  auto original_uy_curr = toHost(wavefield.m_uynGlobalCurr);
  auto original_uz_prev = toHost(wavefield.m_uznGlobalPrev);
  auto original_uz_curr = toHost(wavefield.m_uznGlobalCurr);

  // Swap twice should restore original state
  wavefield.swap();
  wavefield.swap();

  // Verify restoration
  auto ux_prev = toHost(wavefield.m_uxnGlobalPrev);
  auto ux_curr = toHost(wavefield.m_uxnGlobalCurr);
  auto uy_prev = toHost(wavefield.m_uynGlobalPrev);
  auto uy_curr = toHost(wavefield.m_uynGlobalCurr);
  auto uz_prev = toHost(wavefield.m_uznGlobalPrev);
  auto uz_curr = toHost(wavefield.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_prev(i), original_ux_prev(i));
    EXPECT_FLOAT_EQ(ux_curr(i), original_ux_curr(i));
    EXPECT_FLOAT_EQ(uy_prev(i), original_uy_prev(i));
    EXPECT_FLOAT_EQ(uy_curr(i), original_uy_curr(i));
    EXPECT_FLOAT_EQ(uz_prev(i), original_uz_prev(i));
    EXPECT_FLOAT_EQ(uz_curr(i), original_uz_curr(i));
  }
}

TEST_F(WavefieldElasticTest, SwapWithModification) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Modify current fields
  setAt(wavefield.m_uxnGlobalCurr, 123.456f, 5);
  setAt(wavefield.m_uynGlobalCurr, 234.567f, 7);
  setAt(wavefield.m_uznGlobalCurr, 345.678f, 9);

  // Swap
  wavefield.swap();

  // The modified values should now be in the previous fields
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uxnGlobalPrev, 5), 123.456f);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uxnGlobalCurr, 5), 5.0f);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uynGlobalPrev, 7), 234.567f);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uynGlobalCurr, 7), 21.0f);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uznGlobalPrev, 9), 345.678f);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uznGlobalCurr, 9), 45.0f);
}

TEST_F(WavefieldElasticTest, CopyConstructorAfterSwap) {
  WavefieldElastic original(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Swap the original
  original.swap();

  // Create a copy after swap
  WavefieldElastic copy(original);

  // Verify copy has the swapped state
  auto ux_prev = toHost(copy.m_uxnGlobalPrev);
  auto ux_curr = toHost(copy.m_uxnGlobalCurr);
  auto uy_prev = toHost(copy.m_uynGlobalPrev);
  auto uy_curr = toHost(copy.m_uynGlobalCurr);
  auto uz_prev = toHost(copy.m_uznGlobalPrev);
  auto uz_curr = toHost(copy.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_prev(i), i * 2);
    EXPECT_FLOAT_EQ(ux_curr(i), i);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 4);
    EXPECT_FLOAT_EQ(uy_curr(i), i * 3);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 6);
    EXPECT_FLOAT_EQ(uz_curr(i), i * 5);
  }
}

TEST_F(WavefieldElasticTest, EmptyFields) {
  auto emptyUxPrev = allocateVector<vectorReal>(0, "emptyUxPrev");
  auto emptyUxCurr = allocateVector<vectorReal>(0, "emptyUxCurr");
  auto emptyUyPrev = allocateVector<vectorReal>(0, "emptyUyPrev");
  auto emptyUyCurr = allocateVector<vectorReal>(0, "emptyUyCurr");
  auto emptyUzPrev = allocateVector<vectorReal>(0, "emptyUzPrev");
  auto emptyUzCurr = allocateVector<vectorReal>(0, "emptyUzCurr");

  WavefieldElastic wavefield(emptyUxPrev, emptyUxCurr, emptyUyPrev, emptyUyCurr, emptyUzPrev, emptyUzCurr);

  EXPECT_EQ(wavefield.m_uxnGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_uxnGlobalCurr.extent(0), 0);
  EXPECT_EQ(wavefield.m_uynGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_uynGlobalCurr.extent(0), 0);
  EXPECT_EQ(wavefield.m_uznGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_uznGlobalCurr.extent(0), 0);

  // Swap should work with empty fields
  wavefield.swap();
  EXPECT_EQ(wavefield.m_uxnGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_uxnGlobalCurr.extent(0), 0);
  EXPECT_EQ(wavefield.m_uynGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_uynGlobalCurr.extent(0), 0);
  EXPECT_EQ(wavefield.m_uznGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_uznGlobalCurr.extent(0), 0);
}

TEST_F(WavefieldElasticTest, CopyInContainerClass) {
  // Create a simple container class that stores wavefield by copy
  struct WavefieldContainer {
    WavefieldElastic wavefield;

    WavefieldContainer(const WavefieldElastic& wf) : wavefield(wf) {}

    void swap() { wavefield.swap(); }
  };

  WavefieldElastic original(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Store original values
  auto original_ux_prev = toHost(original.m_uxnGlobalPrev);
  auto original_ux_curr = toHost(original.m_uxnGlobalCurr);
  auto original_uy_prev = toHost(original.m_uynGlobalPrev);
  auto original_uy_curr = toHost(original.m_uynGlobalCurr);
  auto original_uz_prev = toHost(original.m_uznGlobalPrev);
  auto original_uz_curr = toHost(original.m_uznGlobalCurr);

  // Create container with wavefield copy
  WavefieldContainer container(original);

  // Verify container has correct initial state
  auto ux_prev = toHost(container.wavefield.m_uxnGlobalPrev);
  auto ux_curr = toHost(container.wavefield.m_uxnGlobalCurr);
  auto uy_prev = toHost(container.wavefield.m_uynGlobalPrev);
  auto uy_curr = toHost(container.wavefield.m_uynGlobalCurr);
  auto uz_prev = toHost(container.wavefield.m_uznGlobalPrev);
  auto uz_curr = toHost(container.wavefield.m_uznGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_prev(i), original_ux_prev(i));
    EXPECT_FLOAT_EQ(ux_curr(i), original_ux_curr(i));
    EXPECT_FLOAT_EQ(uy_prev(i), original_uy_prev(i));
    EXPECT_FLOAT_EQ(uy_curr(i), original_uy_curr(i));
    EXPECT_FLOAT_EQ(uz_prev(i), original_uz_prev(i));
    EXPECT_FLOAT_EQ(uz_curr(i), original_uz_curr(i));
  }

  // Modify original wavefield
  setAt(original.m_uxnGlobalCurr, 999.0f, 10);
  setAt(original.m_uxnGlobalPrev, 888.0f, 20);
  setAt(original.m_uynGlobalCurr, 777.0f, 15);

  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalCurr, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalPrev, 20), 888.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uynGlobalCurr, 15), 777.0f);

  // Perform multiple swaps on container, checking after each
  container.swap();
  // After first swap
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalPrev, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalCurr, 20), 888.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uynGlobalPrev, 15), 777.0f);

  container.swap();
  // After second swap (back to original)
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalCurr, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalPrev, 20), 888.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uynGlobalCurr, 15), 777.0f);

  container.swap();
  // After third swap
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalPrev, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalCurr, 20), 888.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uynGlobalPrev, 15), 777.0f);

  // Swap container's wavefield independently
  container.swap();

  // Now container should be back to original order
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalCurr, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uxnGlobalPrev, 20), 888.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_uynGlobalCurr, 15), 777.0f);

  // Verify they still share the same underlying data
  setAt(container.wavefield.m_uznGlobalCurr, 666.0f, 30);
  EXPECT_FLOAT_EQ(getAt(original.m_uznGlobalCurr, 30), 666.0f);
}

TEST_F(WavefieldElasticTest, SwapWithThreeBuffersRotatesCorrectly) {
  // Allocate prevprev buffers initialised to 10*component_index
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  Kokkos::deep_copy(uxPrevPrev, 10.0f);
  Kokkos::deep_copy(uyPrevPrev, 20.0f);
  Kokkos::deep_copy(uzPrevPrev, 30.0f);

  WavefieldElastic wavefield(uxPrevPrev, uxPrevField, uxCurrField, uyPrevPrev, uyPrevField, uyCurrField, uzPrevPrev,
                             uzPrevField, uzCurrField);
  EXPECT_TRUE(wavefield.hasPrevPrev());

  wavefield.swap();

  // After 3-way rotation:
  //   curr      ← old prevPrev  (ux=10, uy=20, uz=30)
  //   prev      ← old curr      (ux=i*2, uy=i*4, uz=i*6)
  //   prevPrev  ← old prev      (ux=i,   uy=i*3, uz=i*5)
  auto ux_curr = toHost(wavefield.m_uxnGlobalCurr);
  auto uy_curr = toHost(wavefield.m_uynGlobalCurr);
  auto uz_curr = toHost(wavefield.m_uznGlobalCurr);
  auto ux_prev = toHost(wavefield.m_uxnGlobalPrev);
  auto uy_prev = toHost(wavefield.m_uynGlobalPrev);
  auto uz_prev = toHost(wavefield.m_uznGlobalPrev);
  auto ux_prevprev = toHost(wavefield.getPrevPrevField(0));
  auto uy_prevprev = toHost(wavefield.getPrevPrevField(1));
  auto uz_prevprev = toHost(wavefield.getPrevPrevField(2));
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(ux_curr(i), 10.0f);
    EXPECT_FLOAT_EQ(uy_curr(i), 20.0f);
    EXPECT_FLOAT_EQ(uz_curr(i), 30.0f);

    EXPECT_FLOAT_EQ(ux_prev(i), i * 2);
    EXPECT_FLOAT_EQ(uy_prev(i), i * 4);
    EXPECT_FLOAT_EQ(uz_prev(i), i * 6);

    EXPECT_FLOAT_EQ(ux_prevprev(i), static_cast<float>(i));
    EXPECT_FLOAT_EQ(uy_prevprev(i), static_cast<float>(i) * 3);
    EXPECT_FLOAT_EQ(uz_prevprev(i), static_cast<float>(i) * 5);
  }
}

TEST_F(WavefieldElasticTest, SwapThreeTimesRestoresStateWithThreeBuffers) {
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  Kokkos::deep_copy(uxPrevPrev, 10.0f);
  Kokkos::deep_copy(uyPrevPrev, 20.0f);
  Kokkos::deep_copy(uzPrevPrev, 30.0f);

  float initialUxPrev0 = getAt(uxPrevField, 0);
  float initialUxCurr0 = getAt(uxCurrField, 0);
  float initialUxPP0 = getAt(uxPrevPrev, 0);

  WavefieldElastic wavefield(uxPrevPrev, uxPrevField, uxCurrField, uyPrevPrev, uyPrevField, uyCurrField, uzPrevPrev,
                             uzPrevField, uzCurrField);
  EXPECT_TRUE(wavefield.hasPrevPrev());

  wavefield.swap();
  wavefield.swap();
  wavefield.swap();

  EXPECT_FLOAT_EQ(getAt(wavefield.m_uxnGlobalPrev, 0), initialUxPrev0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_uxnGlobalCurr, 0), initialUxCurr0);
  EXPECT_FLOAT_EQ(getAt(wavefield.getPrevPrevField(0), 0), initialUxPP0);
}

TEST_F(WavefieldElasticTest, SwapThreeBuffersNoDataCopy) {
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  Kokkos::deep_copy(uxPrevPrev, 10.0f);
  Kokkos::deep_copy(uyPrevPrev, 20.0f);
  Kokkos::deep_copy(uzPrevPrev, 30.0f);

  WavefieldElastic wavefield(uxPrevPrev, uxPrevField, uxCurrField, uyPrevPrev, uyPrevField, uyCurrField, uzPrevPrev,
                             uzPrevField, uzCurrField);
  wavefield.swap();

  // getPrevPrevField(0) now aliases uxPrevField; write through one and read through the other
  setAt(wavefield.getPrevPrevField(0), 999.0f, 0);
  EXPECT_FLOAT_EQ(getAt(uxPrevField, 0), 999.0f);
}

TEST_F(WavefieldElasticTest, GetNumFieldsReturnsThree) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);
  EXPECT_EQ(wavefield.getNumFields(), 3);
}

TEST_F(WavefieldElasticTest, GetFieldNamesMatchComponents) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);
  const char* const* names = wavefield.getFieldNames();
  ASSERT_NE(names, nullptr);
  EXPECT_STREQ(names[0], "ux");
  EXPECT_STREQ(names[1], "uy");
  EXPECT_STREQ(names[2], "uz");
}

TEST_F(WavefieldElasticTest, PrintDoesNotCrash) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);
  testing::internal::CaptureStdout();
  EXPECT_NO_THROW(wavefield.print());
  testing::internal::GetCapturedStdout();
}

}  // namespace test
}  // namespace fe
}  // namespace solver

#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {
namespace test {

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
    for (size_t i = 0; i < size1; ++i) {
      uxPrevField(i) = i;
      uxCurrField(i) = i * 2;
      uyPrevField(i) = i * 3;
      uyCurrField(i) = i * 4;
      uzPrevField(i) = i * 5;
      uzCurrField(i) = i * 6;
    }

    for (size_t i = 0; i < size2; ++i) {
      uxPrevField2(i) = i * 7;
      uxCurrField2(i) = i * 8;
      uyPrevField2(i) = i * 9;
      uyCurrField2(i) = i * 10;
      uzPrevField2(i) = i * 11;
      uzCurrField2(i) = i * 12;
    }
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
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(i), i);
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(i), i * 2);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalPrev(i), i * 3);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalCurr(i), i * 4);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalPrev(i), i * 5);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalCurr(i), i * 6);
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
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy.m_uxnGlobalPrev(i), original.m_uxnGlobalPrev(i));
    EXPECT_FLOAT_EQ(copy.m_uxnGlobalCurr(i), original.m_uxnGlobalCurr(i));
    EXPECT_FLOAT_EQ(copy.m_uynGlobalPrev(i), original.m_uynGlobalPrev(i));
    EXPECT_FLOAT_EQ(copy.m_uynGlobalCurr(i), original.m_uynGlobalCurr(i));
    EXPECT_FLOAT_EQ(copy.m_uznGlobalPrev(i), original.m_uznGlobalPrev(i));
    EXPECT_FLOAT_EQ(copy.m_uznGlobalCurr(i), original.m_uznGlobalCurr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  original.m_uxnGlobalPrev(0) = 999.0f;
  EXPECT_FLOAT_EQ(copy.m_uxnGlobalPrev(0), 999.0f);
  original.m_uxnGlobalCurr(0) = 888.0f;
  EXPECT_FLOAT_EQ(copy.m_uxnGlobalCurr(0), 888.0f);
  original.m_uynGlobalPrev(0) = 777.0f;
  EXPECT_FLOAT_EQ(copy.m_uynGlobalPrev(0), 777.0f);
  original.m_uynGlobalCurr(0) = 666.0f;
  EXPECT_FLOAT_EQ(copy.m_uynGlobalCurr(0), 666.0f);
  original.m_uznGlobalPrev(0) = 555.0f;
  EXPECT_FLOAT_EQ(copy.m_uznGlobalPrev(0), 555.0f);
  original.m_uznGlobalCurr(0) = 444.0f;
  EXPECT_FLOAT_EQ(copy.m_uznGlobalCurr(0), 444.0f);
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
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield2.m_uxnGlobalPrev(i), wavefield1.m_uxnGlobalPrev(i));
    EXPECT_FLOAT_EQ(wavefield2.m_uxnGlobalCurr(i), wavefield1.m_uxnGlobalCurr(i));
    EXPECT_FLOAT_EQ(wavefield2.m_uynGlobalPrev(i), wavefield1.m_uynGlobalPrev(i));
    EXPECT_FLOAT_EQ(wavefield2.m_uynGlobalCurr(i), wavefield1.m_uynGlobalCurr(i));
    EXPECT_FLOAT_EQ(wavefield2.m_uznGlobalPrev(i), wavefield1.m_uznGlobalPrev(i));
    EXPECT_FLOAT_EQ(wavefield2.m_uznGlobalCurr(i), wavefield1.m_uznGlobalCurr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  wavefield1.m_uxnGlobalPrev(0) = 777.0f;
  EXPECT_FLOAT_EQ(wavefield2.m_uxnGlobalPrev(0), 777.0f);
  wavefield1.m_uxnGlobalCurr(0) = 666.0f;
  EXPECT_FLOAT_EQ(wavefield2.m_uxnGlobalCurr(0), 666.0f);
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

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(i), i);
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(i), i * 2);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalPrev(i), i * 3);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalCurr(i), i * 4);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalPrev(i), i * 5);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalCurr(i), i * 6);
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

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(currentX(i), i * 2);
    EXPECT_FLOAT_EQ(currentY(i), i * 4);
    EXPECT_FLOAT_EQ(currentZ(i), i * 6);
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

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(previousX(i), i);
    EXPECT_FLOAT_EQ(previousY(i), i * 3);
    EXPECT_FLOAT_EQ(previousZ(i), i * 5);
  }
}

TEST_F(WavefieldElasticTest, Swap) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Store original values
  float originalUxPrev0 = wavefield.m_uxnGlobalPrev(0);
  float originalUxCurr0 = wavefield.m_uxnGlobalCurr(0);
  float originalUyPrev0 = wavefield.m_uynGlobalPrev(0);
  float originalUyCurr0 = wavefield.m_uynGlobalCurr(0);
  float originalUzPrev0 = wavefield.m_uznGlobalPrev(0);
  float originalUzCurr0 = wavefield.m_uznGlobalCurr(0);

  // Perform swap
  wavefield.swap();

  // Verify that prev and curr have been swapped
  EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(0), originalUxCurr0);
  EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(0), originalUxPrev0);
  EXPECT_FLOAT_EQ(wavefield.m_uynGlobalPrev(0), originalUyCurr0);
  EXPECT_FLOAT_EQ(wavefield.m_uynGlobalCurr(0), originalUyPrev0);
  EXPECT_FLOAT_EQ(wavefield.m_uznGlobalPrev(0), originalUzCurr0);
  EXPECT_FLOAT_EQ(wavefield.m_uznGlobalCurr(0), originalUzPrev0);

  // Verify all elements were swapped
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(i), i);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalPrev(i), i * 4);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalCurr(i), i * 3);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalPrev(i), i * 6);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalCurr(i), i * 5);
  }
}

TEST_F(WavefieldElasticTest, SwapTwice) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Store original values
  std::vector<float> originalUxPrev(size1), originalUxCurr(size1);
  std::vector<float> originalUyPrev(size1), originalUyCurr(size1);
  std::vector<float> originalUzPrev(size1), originalUzCurr(size1);

  for (size_t i = 0; i < size1; ++i) {
    originalUxPrev[i] = wavefield.m_uxnGlobalPrev(i);
    originalUxCurr[i] = wavefield.m_uxnGlobalCurr(i);
    originalUyPrev[i] = wavefield.m_uynGlobalPrev(i);
    originalUyCurr[i] = wavefield.m_uynGlobalCurr(i);
    originalUzPrev[i] = wavefield.m_uznGlobalPrev(i);
    originalUzCurr[i] = wavefield.m_uznGlobalCurr(i);
  }

  // Swap twice should restore original state
  wavefield.swap();
  wavefield.swap();

  // Verify restoration
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(i), originalUxPrev[i]);
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(i), originalUxCurr[i]);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalPrev(i), originalUyPrev[i]);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalCurr(i), originalUyCurr[i]);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalPrev(i), originalUzPrev[i]);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalCurr(i), originalUzCurr[i]);
  }
}

TEST_F(WavefieldElasticTest, SwapWithModification) {
  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Modify current fields
  wavefield.m_uxnGlobalCurr(5) = 123.456f;
  wavefield.m_uynGlobalCurr(7) = 234.567f;
  wavefield.m_uznGlobalCurr(9) = 345.678f;

  // Swap
  wavefield.swap();

  // The modified values should now be in the previous fields
  EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(5), 123.456f);
  EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(5), 5.0f);
  EXPECT_FLOAT_EQ(wavefield.m_uynGlobalPrev(7), 234.567f);
  EXPECT_FLOAT_EQ(wavefield.m_uynGlobalCurr(7), 21.0f);
  EXPECT_FLOAT_EQ(wavefield.m_uznGlobalPrev(9), 345.678f);
  EXPECT_FLOAT_EQ(wavefield.m_uznGlobalCurr(9), 45.0f);
}

TEST_F(WavefieldElasticTest, CopyConstructorAfterSwap) {
  WavefieldElastic original(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  // Swap the original
  original.swap();

  // Create a copy after swap
  WavefieldElastic copy(original);

  // Verify copy has the swapped state
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy.m_uxnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(copy.m_uxnGlobalCurr(i), i);
    EXPECT_FLOAT_EQ(copy.m_uynGlobalPrev(i), i * 4);
    EXPECT_FLOAT_EQ(copy.m_uynGlobalCurr(i), i * 3);
    EXPECT_FLOAT_EQ(copy.m_uznGlobalPrev(i), i * 6);
    EXPECT_FLOAT_EQ(copy.m_uznGlobalCurr(i), i * 5);
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
  std::vector<float> originalUxPrev(size1), originalUxCurr(size1);
  std::vector<float> originalUyPrev(size1), originalUyCurr(size1);
  std::vector<float> originalUzPrev(size1), originalUzCurr(size1);

  for (size_t i = 0; i < size1; ++i) {
    originalUxPrev[i] = original.m_uxnGlobalPrev(i);
    originalUxCurr[i] = original.m_uxnGlobalCurr(i);
    originalUyPrev[i] = original.m_uynGlobalPrev(i);
    originalUyCurr[i] = original.m_uynGlobalCurr(i);
    originalUzPrev[i] = original.m_uznGlobalPrev(i);
    originalUzCurr[i] = original.m_uznGlobalCurr(i);
  }

  // Create container with wavefield copy
  WavefieldContainer container(original);

  // Verify container has correct initial state
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalPrev(i), originalUxPrev[i]);
    EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalCurr(i), originalUxCurr[i]);
    EXPECT_FLOAT_EQ(container.wavefield.m_uynGlobalPrev(i), originalUyPrev[i]);
    EXPECT_FLOAT_EQ(container.wavefield.m_uynGlobalCurr(i), originalUyCurr[i]);
    EXPECT_FLOAT_EQ(container.wavefield.m_uznGlobalPrev(i), originalUzPrev[i]);
    EXPECT_FLOAT_EQ(container.wavefield.m_uznGlobalCurr(i), originalUzCurr[i]);
  }

  // Modify original wavefield
  original.m_uxnGlobalCurr(10) = 999.0f;
  original.m_uxnGlobalPrev(20) = 888.0f;
  original.m_uynGlobalCurr(15) = 777.0f;

  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalCurr(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalPrev(20), 888.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uynGlobalCurr(15), 777.0f);

  // Perform multiple swaps on container, checking after each
  container.swap();
  // After first swap
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalPrev(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalCurr(20), 888.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uynGlobalPrev(15), 777.0f);

  container.swap();
  // After second swap (back to original)
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalCurr(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalPrev(20), 888.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uynGlobalCurr(15), 777.0f);

  container.swap();
  // After third swap
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalPrev(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalCurr(20), 888.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uynGlobalPrev(15), 777.0f);

  // Swap container's wavefield independently
  container.swap();

  // Now container should be back to original order
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalCurr(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uxnGlobalPrev(20), 888.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_uynGlobalCurr(15), 777.0f);

  // Verify they still share the same underlying data
  container.wavefield.m_uznGlobalCurr(30) = 666.0f;
  EXPECT_FLOAT_EQ(original.m_uznGlobalCurr(30), 666.0f);
}

TEST_F(WavefieldElasticTest, RotateRotatesThreeBuffers) {
  // Allocate prevprev buffers initialised to 10*component_index
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  for (size_t i = 0; i < size1; ++i) {
    uxPrevPrev(i) = 10.0f;
    uyPrevPrev(i) = 20.0f;
    uzPrevPrev(i) = 30.0f;
  }

  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);
  wavefield.rotate(uxPrevPrev, 0);
  wavefield.rotate(uyPrevPrev, 1);
  wavefield.rotate(uzPrevPrev, 2);

  // After rotation:
  //   curr      ← old prevPrev  (ux=10, uy=20, uz=30)
  //   prev      ← old curr      (ux=i*2, uy=i*4, uz=i*6)
  //   prevPrev  ← old prev      (ux=i,   uy=i*3, uz=i*5)
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(i), 10.0f);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalCurr(i), 20.0f);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalCurr(i), 30.0f);

    EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(wavefield.m_uynGlobalPrev(i), i * 4);
    EXPECT_FLOAT_EQ(wavefield.m_uznGlobalPrev(i), i * 6);

    EXPECT_FLOAT_EQ(uxPrevPrev(i), static_cast<float>(i));
    EXPECT_FLOAT_EQ(uyPrevPrev(i), static_cast<float>(i) * 3);
    EXPECT_FLOAT_EQ(uzPrevPrev(i), static_cast<float>(i) * 5);
  }
}

TEST_F(WavefieldElasticTest, RotateThreeTimesRestoresState) {
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  for (size_t i = 0; i < size1; ++i) {
    uxPrevPrev(i) = 10.0f;
    uyPrevPrev(i) = 20.0f;
    uzPrevPrev(i) = 30.0f;
  }

  float initialUxPrev0 = uxPrevField(0);
  float initialUxCurr0 = uxCurrField(0);
  float initialUxPP0 = uxPrevPrev(0);

  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);

  wavefield.rotate(uxPrevPrev, 0);
  wavefield.rotate(uyPrevPrev, 1);
  wavefield.rotate(uzPrevPrev, 2);
  wavefield.rotate(uxPrevPrev, 0);
  wavefield.rotate(uyPrevPrev, 1);
  wavefield.rotate(uzPrevPrev, 2);
  wavefield.rotate(uxPrevPrev, 0);
  wavefield.rotate(uyPrevPrev, 1);
  wavefield.rotate(uzPrevPrev, 2);

  EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalPrev(0), initialUxPrev0);
  EXPECT_FLOAT_EQ(wavefield.m_uxnGlobalCurr(0), initialUxCurr0);
  EXPECT_FLOAT_EQ(uxPrevPrev(0), initialUxPP0);
}

TEST_F(WavefieldElasticTest, RotateNoDataCopy) {
  auto uxPrevPrev = allocateVector<vectorReal>(size1, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size1, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size1, "uzPrevPrev");
  for (size_t i = 0; i < size1; ++i) {
    uxPrevPrev(i) = 10.0f;
    uyPrevPrev(i) = 20.0f;
    uzPrevPrev(i) = 30.0f;
  }

  WavefieldElastic wavefield(uxPrevField, uxCurrField, uyPrevField, uyCurrField, uzPrevField, uzCurrField);
  wavefield.rotate(uxPrevPrev, 0);
  wavefield.rotate(uyPrevPrev, 1);
  wavefield.rotate(uzPrevPrev, 2);

  // uxPrevPrev now aliases uxPrevField; write through one and read through the
  // other
  uxPrevPrev(0) = 999.0f;
  EXPECT_FLOAT_EQ(uxPrevField(0), 999.0f);
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

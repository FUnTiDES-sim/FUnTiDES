#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "device_view_access.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromHost;
using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

class WavefieldAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create test data with specific sizes
    size1 = 100;
    size2 = 200;

    prevField = allocateVector<vectorReal>(size1, "prevField");
    currField = allocateVector<vectorReal>(size1, "currField");
    prevField2 = allocateVector<vectorReal>(size2, "prevField2");
    currField2 = allocateVector<vectorReal>(size2, "currField2");

    // Initialize with test values
    fillFromHost(prevField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i;
    });
    fillFromHost(currField, [&](const auto& h) {
      for (size_t i = 0; i < size1; ++i) h(i) = i * 2;
    });
    fillFromHost(prevField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 3;
    });
    fillFromHost(currField2, [&](const auto& h) {
      for (size_t i = 0; i < size2; ++i) h(i) = i * 4;
    });
  }

  size_t size1;
  size_t size2;
  vectorReal prevField;
  vectorReal currField;
  vectorReal prevField2;
  vectorReal currField2;
};

TEST_F(WavefieldAcousticTest, Constructor) {
  WavefieldAcoustic wavefield(prevField, currField);

  EXPECT_EQ(wavefield.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_pnGlobalCurr.extent(0), size1);

  // Verify data is correctly stored
  auto prev = toHost(wavefield.m_pnGlobalPrev);
  auto curr = toHost(wavefield.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(prev(i), i);
    EXPECT_FLOAT_EQ(curr(i), i * 2);
  }
}

TEST_F(WavefieldAcousticTest, CopyConstructor) {
  WavefieldAcoustic original(prevField, currField);
  WavefieldAcoustic copy(original);

  // Check that copy has the same extent
  EXPECT_EQ(copy.m_pnGlobalPrev.extent(0), original.m_pnGlobalPrev.extent(0));
  EXPECT_EQ(copy.m_pnGlobalCurr.extent(0), original.m_pnGlobalCurr.extent(0));

  // Check that copy has the same data
  auto copy_prev = toHost(copy.m_pnGlobalPrev), original_prev = toHost(original.m_pnGlobalPrev);
  auto copy_curr = toHost(copy.m_pnGlobalCurr), original_curr = toHost(original.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy_prev(i), original_prev(i));
    EXPECT_FLOAT_EQ(copy_curr(i), original_curr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(original.m_pnGlobalPrev, 999.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_pnGlobalPrev, 0), 999.0f);
  setAt(original.m_pnGlobalCurr, 888.0f, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_pnGlobalCurr, 0), 888.0f);
}

TEST_F(WavefieldAcousticTest, CopyAssignmentOperator) {
  WavefieldAcoustic wavefield1(prevField, currField);
  WavefieldAcoustic wavefield2(prevField2, currField2);

  // Verify initial state
  EXPECT_EQ(wavefield2.m_pnGlobalPrev.extent(0), size2);
  EXPECT_EQ(wavefield2.m_pnGlobalCurr.extent(0), size2);

  // Perform assignment
  wavefield2 = wavefield1;

  // Check that wavefield2 now has the same extent as wavefield1
  EXPECT_EQ(wavefield2.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield2.m_pnGlobalCurr.extent(0), size1);

  // Check that data matches
  auto w2_prev = toHost(wavefield2.m_pnGlobalPrev), w1_prev = toHost(wavefield1.m_pnGlobalPrev);
  auto w2_curr = toHost(wavefield2.m_pnGlobalCurr), w1_curr = toHost(wavefield1.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(w2_prev(i), w1_prev(i));
    EXPECT_FLOAT_EQ(w2_curr(i), w1_curr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(wavefield1.m_pnGlobalPrev, 777.0f, 0);
  EXPECT_FLOAT_EQ(getAt(wavefield2.m_pnGlobalPrev, 0), 777.0f);
  setAt(wavefield1.m_pnGlobalCurr, 666.0f, 0);
  EXPECT_FLOAT_EQ(getAt(wavefield2.m_pnGlobalCurr, 0), 666.0f);
}

TEST_F(WavefieldAcousticTest, CopyAssignmentSelfAssignment) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Self-assignment should not cause issues
  wavefield = wavefield;

  // Verify data is unchanged
  EXPECT_EQ(wavefield.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_pnGlobalCurr.extent(0), size1);

  auto prev = toHost(wavefield.m_pnGlobalPrev);
  auto curr = toHost(wavefield.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(prev(i), i);
    EXPECT_FLOAT_EQ(curr(i), i * 2);
  }
}

TEST_F(WavefieldAcousticTest, GetCurrentField) {
  WavefieldAcoustic wavefield(prevField, currField);

  auto current = wavefield.getCurrentField(0);

  EXPECT_EQ(current.extent(0), size1);
  auto h_current = toHost(current);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(h_current(i), i * 2);
  }
}

TEST_F(WavefieldAcousticTest, GetPreviousField) {
  WavefieldAcoustic wavefield(prevField, currField);

  auto previous = wavefield.getPreviousField(0);

  EXPECT_EQ(previous.extent(0), size1);
  auto h_previous = toHost(previous);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(h_previous(i), i);
  }
}

TEST_F(WavefieldAcousticTest, Swap) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Store original values
  float originalPrev0 = getAt(wavefield.m_pnGlobalPrev, 0);
  float originalCurr0 = getAt(wavefield.m_pnGlobalCurr, 0);

  // Perform swap
  wavefield.swap();

  // Verify that prev and curr have been swapped
  EXPECT_FLOAT_EQ(getAt(wavefield.m_pnGlobalPrev, 0), originalCurr0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_pnGlobalCurr, 0), originalPrev0);

  // Verify all elements were swapped
  auto prev = toHost(wavefield.m_pnGlobalPrev);
  auto curr = toHost(wavefield.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(prev(i), i * 2);
    EXPECT_FLOAT_EQ(curr(i), i);
  }
}

TEST_F(WavefieldAcousticTest, SwapTwice) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Store original values
  auto original_prev = toHost(wavefield.m_pnGlobalPrev);
  auto original_curr = toHost(wavefield.m_pnGlobalCurr);

  // Swap twice should restore original state
  wavefield.swap();
  wavefield.swap();

  // Verify restoration
  auto prev = toHost(wavefield.m_pnGlobalPrev);
  auto curr = toHost(wavefield.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(prev(i), original_prev(i));
    EXPECT_FLOAT_EQ(curr(i), original_curr(i));
  }
}

TEST_F(WavefieldAcousticTest, SwapWithModification) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Modify current field
  setAt(wavefield.m_pnGlobalCurr, 123.456f, 5);

  // Swap
  wavefield.swap();

  // The modified value should now be in the previous field
  EXPECT_FLOAT_EQ(getAt(wavefield.m_pnGlobalPrev, 5), 123.456f);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_pnGlobalCurr, 5), 5.0f);
}

TEST_F(WavefieldAcousticTest, CopyConstructorAfterSwap) {
  WavefieldAcoustic original(prevField, currField);

  // Swap the original
  original.swap();

  // Create a copy after swap
  WavefieldAcoustic copy(original);

  // Verify copy has the swapped state
  auto prev = toHost(copy.m_pnGlobalPrev);
  auto curr = toHost(copy.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(prev(i), i * 2);
    EXPECT_FLOAT_EQ(curr(i), i);
  }
}

TEST_F(WavefieldAcousticTest, EmptyFields) {
  auto emptyPrev = allocateVector<vectorReal>(0, "emptyPrev");
  auto emptyCurr = allocateVector<vectorReal>(0, "emptyCurr");

  WavefieldAcoustic wavefield(emptyPrev, emptyCurr);

  EXPECT_EQ(wavefield.m_pnGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_pnGlobalCurr.extent(0), 0);

  // Swap should work with empty fields
  wavefield.swap();
  EXPECT_EQ(wavefield.m_pnGlobalPrev.extent(0), 0);
  EXPECT_EQ(wavefield.m_pnGlobalCurr.extent(0), 0);
}

TEST_F(WavefieldAcousticTest, CopyInContainerClass) {
  // Create a simple container class that stores wavefield by copy
  struct WavefieldContainer {
    WavefieldAcoustic wavefield;

    WavefieldContainer(const WavefieldAcoustic& wf) : wavefield(wf) {}

    void swap() { wavefield.swap(); }
  };

  WavefieldAcoustic original(prevField, currField);

  // Store original values
  auto original_prev = toHost(original.m_pnGlobalPrev);
  auto original_curr = toHost(original.m_pnGlobalCurr);

  // Create container with wavefield copy
  WavefieldContainer container(original);

  // Verify container has correct initial state
  auto prev = toHost(container.wavefield.m_pnGlobalPrev);
  auto curr = toHost(container.wavefield.m_pnGlobalCurr);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(prev(i), original_prev(i));
    EXPECT_FLOAT_EQ(curr(i), original_curr(i));
  }

  // Modify original wavefield
  setAt(original.m_pnGlobalCurr, 999.0f, 10);
  setAt(original.m_pnGlobalPrev, 888.0f, 20);

  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalCurr, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalPrev, 20), 888.0f);

  // Perform multiple swaps on original, checking after each
  container.swap();
  // After first swap
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalPrev, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalCurr, 20), 888.0f);

  container.swap();
  // After second swap (back to original)
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalCurr, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalPrev, 20), 888.0f);

  container.swap();
  // After third swap
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalPrev, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalCurr, 20), 888.0f);

  // Swap container's wavefield independently
  container.swap();

  // Now container should be back to original order
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalCurr, 10), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.wavefield.m_pnGlobalPrev, 20), 888.0f);

  // Verify they still share the same underlying data
  setAt(container.wavefield.m_pnGlobalCurr, 777.0f, 30);
  EXPECT_FLOAT_EQ(getAt(original.m_pnGlobalCurr, 30), 777.0f);
}

TEST_F(WavefieldAcousticTest, SwapWithThreeBuffersRotatesCorrectly) {
  // Test backward mode (3-buffer): prevprev = {10}, prev = {i}, curr = {i*2}
  auto prevPrev = allocateVector<vectorReal>(size1, "prevPrev");
  Kokkos::deep_copy(prevPrev, 10.0f);

  WavefieldAcoustic wavefield(prevPrev, prevField, currField);
  EXPECT_TRUE(wavefield.hasPrevPrev());

  wavefield.swap();

  // After 3-way rotation:
  //   curr      ← old prevPrev  (value = 10.0)
  //   prev      ← old curr      (value = i*2)
  //   prevPrev  ← old prev      (value = i)
  auto curr = toHost(wavefield.m_pnGlobalCurr);
  auto prev = toHost(wavefield.m_pnGlobalPrev);
  auto prevprev = toHost(wavefield.getPrevPrevField(0));
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(curr(i), 10.0f);
    EXPECT_FLOAT_EQ(prev(i), i * 2);
    EXPECT_FLOAT_EQ(prevprev(i), static_cast<float>(i));
  }
}

TEST_F(WavefieldAcousticTest, SwapThreeTimesRestoresStateWithThreeBuffers) {
  auto prevPrev = allocateVector<vectorReal>(size1, "prevPrev");
  Kokkos::deep_copy(prevPrev, 10.0f);

  // Record initial values
  float initialPrev0 = getAt(prevField, 0);     // i=0 → 0
  float initialCurr0 = getAt(currField, 0);     // i=0 → 0
  float initialPrevPrev0 = getAt(prevPrev, 0);  // 10.0

  WavefieldAcoustic wavefield(prevPrev, prevField, currField);
  EXPECT_TRUE(wavefield.hasPrevPrev());

  // Three swaps on a 3-element ring must restore the original assignment
  wavefield.swap();
  wavefield.swap();
  wavefield.swap();

  EXPECT_FLOAT_EQ(getAt(wavefield.m_pnGlobalPrev, 0), initialPrev0);
  EXPECT_FLOAT_EQ(getAt(wavefield.m_pnGlobalCurr, 0), initialCurr0);
  EXPECT_FLOAT_EQ(getAt(wavefield.getPrevPrevField(0), 0), initialPrevPrev0);
}

TEST_F(WavefieldAcousticTest, SwapThreeBuffersNoDataCopy) {
  // Verifies that swap is view-handle only: no data copying
  auto prevPrev = allocateVector<vectorReal>(size1, "prevPrev");
  Kokkos::deep_copy(prevPrev, 10.0f);

  WavefieldAcoustic wavefield(prevPrev, prevField, currField);
  wavefield.swap();

  // After swap, getPrevPrevField(0) now points to what was prevField
  setAt(wavefield.getPrevPrevField(0), 999.0f, 0);
  // prevField is the same underlying allocation → must see the change
  EXPECT_FLOAT_EQ(getAt(prevField, 0), 999.0f);
}

TEST_F(WavefieldAcousticTest, GetNumFieldsReturnsOne) {
  WavefieldAcoustic wavefield(prevField, currField);
  EXPECT_EQ(wavefield.getNumFields(), 1);
}

TEST_F(WavefieldAcousticTest, GetFieldNamesMatchComponents) {
  WavefieldAcoustic wavefield(prevField, currField);
  const char* const* names = wavefield.getFieldNames();
  ASSERT_NE(names, nullptr);
  EXPECT_STREQ(names[0], "pressure");
}

TEST_F(WavefieldAcousticTest, PrintDoesNotCrash) {
  WavefieldAcoustic wavefield(prevField, currField);
  testing::internal::CaptureStdout();
  EXPECT_NO_THROW(wavefield.print());
  testing::internal::GetCapturedStdout();
}

}  // namespace test
}  // namespace fe
}  // namespace solver
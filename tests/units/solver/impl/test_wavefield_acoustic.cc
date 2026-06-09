#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {
namespace test {

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
    for (size_t i = 0; i < size1; ++i) {
      prevField(i) = i;
      currField(i) = i * 2;
    }

    for (size_t i = 0; i < size2; ++i) {
      prevField2(i) = i * 3;
      currField2(i) = i * 4;
    }
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
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(i), i);
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(i), i * 2);
  }
}

TEST_F(WavefieldAcousticTest, CopyConstructor) {
  WavefieldAcoustic original(prevField, currField);
  WavefieldAcoustic copy(original);

  // Check that copy has the same extent
  EXPECT_EQ(copy.m_pnGlobalPrev.extent(0), original.m_pnGlobalPrev.extent(0));
  EXPECT_EQ(copy.m_pnGlobalCurr.extent(0), original.m_pnGlobalCurr.extent(0));

  // Check that copy has the same data
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy.m_pnGlobalPrev(i), original.m_pnGlobalPrev(i));
    EXPECT_FLOAT_EQ(copy.m_pnGlobalCurr(i), original.m_pnGlobalCurr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  original.m_pnGlobalPrev(0) = 999.0f;
  EXPECT_FLOAT_EQ(copy.m_pnGlobalPrev(0), 999.0f);
  original.m_pnGlobalCurr(0) = 888.0f;
  EXPECT_FLOAT_EQ(copy.m_pnGlobalCurr(0), 888.0f);
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
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield2.m_pnGlobalPrev(i), wavefield1.m_pnGlobalPrev(i));
    EXPECT_FLOAT_EQ(wavefield2.m_pnGlobalCurr(i), wavefield1.m_pnGlobalCurr(i));
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  wavefield1.m_pnGlobalPrev(0) = 777.0f;
  EXPECT_FLOAT_EQ(wavefield2.m_pnGlobalPrev(0), 777.0f);
  wavefield1.m_pnGlobalCurr(0) = 666.0f;
  EXPECT_FLOAT_EQ(wavefield2.m_pnGlobalCurr(0), 666.0f);
}

TEST_F(WavefieldAcousticTest, CopyAssignmentSelfAssignment) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Self-assignment should not cause issues
  wavefield = wavefield;

  // Verify data is unchanged
  EXPECT_EQ(wavefield.m_pnGlobalPrev.extent(0), size1);
  EXPECT_EQ(wavefield.m_pnGlobalCurr.extent(0), size1);

  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(i), i);
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(i), i * 2);
  }
}

TEST_F(WavefieldAcousticTest, GetCurrentField) {
  WavefieldAcoustic wavefield(prevField, currField);

  auto current = wavefield.getCurrentField(0);

  EXPECT_EQ(current.extent(0), size1);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(current(i), i * 2);
  }
}

TEST_F(WavefieldAcousticTest, GetPreviousField) {
  WavefieldAcoustic wavefield(prevField, currField);

  auto previous = wavefield.getPreviousField(0);

  EXPECT_EQ(previous.extent(0), size1);
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(previous(i), i);
  }
}

TEST_F(WavefieldAcousticTest, Swap) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Store original values
  float originalPrev0 = wavefield.m_pnGlobalPrev(0);
  float originalCurr0 = wavefield.m_pnGlobalCurr(0);

  // Perform swap
  wavefield.swap();

  // Verify that prev and curr have been swapped
  EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(0), originalCurr0);
  EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(0), originalPrev0);

  // Verify all elements were swapped
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(i), i);
  }
}

TEST_F(WavefieldAcousticTest, SwapTwice) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Store original values
  std::vector<float> originalPrev(size1);
  std::vector<float> originalCurr(size1);
  for (size_t i = 0; i < size1; ++i) {
    originalPrev[i] = wavefield.m_pnGlobalPrev(i);
    originalCurr[i] = wavefield.m_pnGlobalCurr(i);
  }

  // Swap twice should restore original state
  wavefield.swap();
  wavefield.swap();

  // Verify restoration
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(i), originalPrev[i]);
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(i), originalCurr[i]);
  }
}

TEST_F(WavefieldAcousticTest, SwapWithModification) {
  WavefieldAcoustic wavefield(prevField, currField);

  // Modify current field
  wavefield.m_pnGlobalCurr(5) = 123.456f;

  // Swap
  wavefield.swap();

  // The modified value should now be in the previous field
  EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(5), 123.456f);
  EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(5), 5.0f);
}

TEST_F(WavefieldAcousticTest, CopyConstructorAfterSwap) {
  WavefieldAcoustic original(prevField, currField);

  // Swap the original
  original.swap();

  // Create a copy after swap
  WavefieldAcoustic copy(original);

  // Verify copy has the swapped state
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(copy.m_pnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(copy.m_pnGlobalCurr(i), i);
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
  std::vector<float> originalPrev(size1);
  std::vector<float> originalCurr(size1);
  for (size_t i = 0; i < size1; ++i) {
    originalPrev[i] = original.m_pnGlobalPrev(i);
    originalCurr[i] = original.m_pnGlobalCurr(i);
  }

  // Create container with wavefield copy
  WavefieldContainer container(original);

  // Verify container has correct initial state
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalPrev(i), originalPrev[i]);
    EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalCurr(i), originalCurr[i]);
  }

  // Modify original wavefield
  original.m_pnGlobalCurr(10) = 999.0f;
  original.m_pnGlobalPrev(20) = 888.0f;

  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalCurr(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalPrev(20), 888.0f);

  // Perform multiple swaps on original, checking after each
  container.swap();
  // After first swap
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalPrev(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalCurr(20), 888.0f);

  container.swap();
  // After second swap (back to original)
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalCurr(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalPrev(20), 888.0f);

  container.swap();
  // After third swap
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalPrev(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalCurr(20), 888.0f);

  // Swap container's wavefield independently
  container.swap();

  // Now container should be back to original order
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalCurr(10), 999.0f);
  EXPECT_FLOAT_EQ(container.wavefield.m_pnGlobalPrev(20), 888.0f);

  // Verify they still share the same underlying data
  container.wavefield.m_pnGlobalCurr(30) = 777.0f;
  EXPECT_FLOAT_EQ(original.m_pnGlobalCurr(30), 777.0f);
}

TEST_F(WavefieldAcousticTest, RotateRotatesThreeBuffers) {
  // prevprev = {10, 10, ...}, prev = {i}, curr = {i*2}
  auto prevPrev = allocateVector<vectorReal>(size1, "prevPrev");
  for (size_t i = 0; i < size1; ++i) prevPrev(i) = 10.0f;

  WavefieldAcoustic wavefield(prevField, currField);
  wavefield.rotate(prevPrev, 0);

  // After rotation:
  //   curr      ← old prevPrev  (value = 10.0)
  //   prev      ← old curr      (value = i*2)
  //   prevPrev  ← old prev      (value = i)
  for (size_t i = 0; i < size1; ++i) {
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(i), 10.0f);
    EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(i), i * 2);
    EXPECT_FLOAT_EQ(prevPrev(i), static_cast<float>(i));
  }
}

TEST_F(WavefieldAcousticTest, RotateThreeTimesRestoresState) {
  auto prevPrev = allocateVector<vectorReal>(size1, "prevPrev");
  for (size_t i = 0; i < size1; ++i) prevPrev(i) = 10.0f;

  // Record identity of underlying data before rotation
  float initialPrev0 = prevField(0);     // i=0 → 0
  float initialCurr0 = currField(0);     // i=0 → 0
  float initialPrevPrev0 = prevPrev(0);  // 10.0

  WavefieldAcoustic wavefield(prevField, currField);

  // Three rotations on a 3-element ring must restore the original assignment
  wavefield.rotate(prevPrev, 0);
  wavefield.rotate(prevPrev, 0);
  wavefield.rotate(prevPrev, 0);

  EXPECT_FLOAT_EQ(wavefield.m_pnGlobalPrev(0), initialPrev0);
  EXPECT_FLOAT_EQ(wavefield.m_pnGlobalCurr(0), initialCurr0);
  EXPECT_FLOAT_EQ(prevPrev(0), initialPrevPrev0);
}

TEST_F(WavefieldAcousticTest, RotateNoDataCopy) {
  // Verifies that the rotation is view-handle only: mutating the underlying
  // buffer is visible through the rotated handle.
  auto prevPrev = allocateVector<vectorReal>(size1, "prevPrev");
  for (size_t i = 0; i < size1; ++i) prevPrev(i) = 10.0f;

  WavefieldAcoustic wavefield(prevField, currField);
  wavefield.rotate(prevPrev, 0);

  // prevPrev now points to what was prevField; write through prevPrev
  prevPrev(0) = 999.0f;
  // prevField is the same underlying allocation → must see the change
  EXPECT_FLOAT_EQ(prevField(0), 999.0f);
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
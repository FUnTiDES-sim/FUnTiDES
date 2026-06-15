#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

#include "data_type.h"
#include "rhs_acoustic.h"
#include "sem_solver_data.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {
namespace testing {

// ============================================================================
// Print Method Tests (Using standard cout redirection)
// ============================================================================

class SEMSolverDataAcousticPrintTest : public ::testing::Test {
 protected:
  std::stringstream buffer;
  std::streambuf* originalCoutBuffer;

  void SetUp() override {
    // Redirect std::cout to the local buffer before each test
    originalCoutBuffer = std::cout.rdbuf(buffer.rdbuf());
  }

  void TearDown() override {
    // Restore the original std::cout buffer after each test
    std::cout.rdbuf(originalCoutBuffer);
  }
};

TEST_F(SEMSolverDataAcousticPrintTest, PrintMethodOutputsAcousticInfo) {
  WavefieldAcoustic mock_wavefield;
  RhsAcoustic mock_rhs;

  SEMsolverDataAcoustic solver_data(mock_wavefield, mock_rhs);

  solver_data.print();

  std::string output = buffer.str();

  // Verify key components of the output string
  EXPECT_TRUE(output.find("SEMsolverData<") != std::string::npos);
  EXPECT_TRUE(output.find("Field[0]") != std::string::npos);
  EXPECT_TRUE(output.find("size:") != std::string::npos);
  EXPECT_TRUE(output.find("RHS[0] size:") != std::string::npos);
  EXPECT_TRUE(output.find("RHS Element size:") != std::string::npos);
  EXPECT_TRUE(output.find("RHS Weights size:") != std::string::npos);
}

TEST_F(SEMSolverDataAcousticPrintTest, PrintMethodWithThreeBuffers) {
  size_t size = 50;
  auto prevPrev = allocateVector<vectorReal>(size, "prevPrev");
  auto prevField = allocateVector<vectorReal>(size, "prevField");
  auto currField = allocateVector<vectorReal>(size, "currField");

  WavefieldAcoustic wavefield(prevPrev, prevField, currField);
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  EXPECT_NO_THROW(solver_data.print());

  std::string output = buffer.str();

  // Verify key components are present
  EXPECT_TRUE(output.find("SEMsolverData<") != std::string::npos);
  EXPECT_TRUE(output.find("Field[0]") != std::string::npos);
  EXPECT_TRUE(output.find("size:") != std::string::npos);
}

// ============================================================================
// General Acoustic Tests
// ============================================================================

TEST(SEMsolverDataAcousticTest, InitializesAndDelegatesCorrectly) {
  WavefieldAcoustic wavefield;
  RhsAcoustic rhs;

  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  // Verify delegation by comparing the extent of the returned views
  EXPECT_EQ(solver_data.getCurrentField(0).extent(0), wavefield.getCurrentField(0).extent(0));

  EXPECT_EQ(solver_data.getPreviousField(0).extent(0), wavefield.getPreviousField(0).extent(0));

  EXPECT_EQ(solver_data.getRhsTerm(0).extent(0), rhs.getTerm(0).extent(0));

  EXPECT_EQ(solver_data.getRhsElement().extent(0), rhs.getElement().extent(0));

  EXPECT_EQ(solver_data.getRhsWeights().extent(0), rhs.getWeights().extent(0));
}

// ============================================================================
// Acoustic Swap and Rotate Tests
// ============================================================================

class SEMSolverDataAcousticSwapRotateTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size = 50;
    prevField = allocateVector<vectorReal>(size, "prevField");
    currField = allocateVector<vectorReal>(size, "currField");
    for (size_t i = 0; i < size; ++i) {
      prevField(i) = static_cast<float>(i);
      currField(i) = static_cast<float>(i) * 2.0f;
    }
  }

  size_t size;
  vectorReal prevField;
  vectorReal currField;
};

TEST_F(SEMSolverDataAcousticSwapRotateTest, SwapExchangesPrevAndCurrFields) {
  WavefieldAcoustic wavefield(prevField, currField);
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  solver_data.swapWavefields();

  auto curr = solver_data.getCurrentField(0);
  auto prev = solver_data.getPreviousField(0);
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(curr(i), static_cast<float>(i));         // old prev
    EXPECT_FLOAT_EQ(prev(i), static_cast<float>(i) * 2.0f);  // old curr
  }
}

TEST_F(SEMSolverDataAcousticSwapRotateTest, SwapTwiceRestoresOriginalValues) {
  WavefieldAcoustic wavefield(prevField, currField);
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  solver_data.swapWavefields();
  solver_data.swapWavefields();

  auto curr = solver_data.getCurrentField(0);
  auto prev = solver_data.getPreviousField(0);
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(curr(i), static_cast<float>(i) * 2.0f);  // original curr
    EXPECT_FLOAT_EQ(prev(i), static_cast<float>(i));         // original prev
  }
}

TEST_F(SEMSolverDataAcousticSwapRotateTest, SwapWithThreeBuffersRotatesCorrectly) {
  auto prevPrev = allocateVector<vectorReal>(size, "prevPrev");
  for (size_t i = 0; i < size; ++i) prevPrev(i) = 10.0f;

  WavefieldAcoustic wavefield(prevPrev, prevField, currField);
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  solver_data.swapWavefields();

  // After 3-way rotation:
  //   curr     <- old prevPrev (10.0)
  //   prev     <- old curr     (i*2)
  //   prevPrev <- old prev     (i)
  auto curr = solver_data.getCurrentField(0);
  auto prev = solver_data.getPreviousField(0);
  auto pp = solver_data.getPrevPrevField(0);
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(curr(i), 10.0f);
    EXPECT_FLOAT_EQ(prev(i), static_cast<float>(i) * 2.0f);
    EXPECT_FLOAT_EQ(pp(i), static_cast<float>(i));
  }
}

TEST_F(SEMSolverDataAcousticSwapRotateTest, SwapThreeTimesRestoresStateWithThreeBuffers) {
  auto prevPrev = allocateVector<vectorReal>(size, "prevPrev");
  for (size_t i = 0; i < size; ++i) prevPrev(i) = 10.0f;

  WavefieldAcoustic wavefield(prevPrev, prevField, currField);
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  solver_data.swapWavefields();
  solver_data.swapWavefields();
  solver_data.swapWavefields();

  auto curr = solver_data.getCurrentField(0);
  auto prev = solver_data.getPreviousField(0);
  auto pp = solver_data.getPrevPrevField(0);
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(curr(i), static_cast<float>(i) * 2.0f);  // original curr
    EXPECT_FLOAT_EQ(prev(i), static_cast<float>(i));         // original prev
    EXPECT_FLOAT_EQ(pp(i), 10.0f);                           // original prevPrev
  }
}

}  // namespace testing
}  // namespace fe
}  // namespace solver

#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

#include "data_type.h"
#include "rhs_acoustic.h"
#include "rhs_elastic.h"
#include "sem_solver_data.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {
namespace testing {

// ============================================================================
// Print Method Tests (Using standard cout redirection)
// ============================================================================

class SEMSolverDataPrintTest : public ::testing::Test {
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

TEST_F(SEMSolverDataPrintTest, PrintMethodOutputsAcousticInfo) {
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

TEST_F(SEMSolverDataPrintTest, PrintMethodOutputsElasticInfo) {
  WavefieldElastic mock_wavefield;
  RhsElastic mock_rhs;
  bool is_distributed = true;

  SEMsolverDataElastic solver_data(mock_wavefield, mock_rhs, is_distributed);

  solver_data.print();

  std::string output = buffer.str();

  EXPECT_TRUE(output.find("SEMsolverData<") != std::string::npos);
  EXPECT_TRUE(output.find("Field[0]") != std::string::npos);
  EXPECT_TRUE(output.find("Field[1]") != std::string::npos);
  EXPECT_TRUE(output.find("Field[2]") != std::string::npos);
  EXPECT_TRUE(output.find("RHS Element size:") != std::string::npos);
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
// General Elastic Tests
// ============================================================================

TEST(SEMsolverDataElasticTest, InitializesAndDelegatesCorrectly) {
  WavefieldElastic wavefield;
  RhsElastic rhs;

  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  for (int i = 0; i < SEMsolverDataElastic::kNumFields; ++i) {
    EXPECT_EQ(solver_data.getCurrentField(i).extent(0), wavefield.getCurrentField(i).extent(0));
  }

  EXPECT_FALSE(solver_data.isDistributed);
}

TEST(SEMsolverDataElasticTest, RespectsDistributedFlag) {
  WavefieldElastic wavefield;
  RhsElastic rhs;

  SEMsolverDataElastic solver_data(wavefield, rhs, true);

  EXPECT_TRUE(solver_data.isDistributed);
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

TEST_F(SEMSolverDataAcousticSwapRotateTest, RotateRotatesThreeBuffers) {
  auto prevPrev = allocateVector<vectorReal>(size, "prevPrev");
  for (size_t i = 0; i < size; ++i) prevPrev(i) = 10.0f;

  WavefieldAcoustic wavefield(prevField, currField);
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  solver_data.rotateWavefields(prevPrev, 0);

  // After rotation:
  //   curr     <- old prevPrev (10.0)
  //   prev     <- old curr     (i*2)
  //   prevPrev <- old prev     (i)
  auto curr = solver_data.getCurrentField(0);
  auto prev = solver_data.getPreviousField(0);
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(curr(i), 10.0f);
    EXPECT_FLOAT_EQ(prev(i), static_cast<float>(i) * 2.0f);
    EXPECT_FLOAT_EQ(prevPrev(i), static_cast<float>(i));
  }
}

TEST_F(SEMSolverDataAcousticSwapRotateTest, RotateThreeTimesRestoresState) {
  auto prevPrev = allocateVector<vectorReal>(size, "prevPrev");
  for (size_t i = 0; i < size; ++i) prevPrev(i) = 10.0f;

  WavefieldAcoustic wavefield(prevField, currField);
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  solver_data.rotateWavefields(prevPrev, 0);
  solver_data.rotateWavefields(prevPrev, 0);
  solver_data.rotateWavefields(prevPrev, 0);

  auto curr = solver_data.getCurrentField(0);
  auto prev = solver_data.getPreviousField(0);
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(curr(i), static_cast<float>(i) * 2.0f);  // original curr
    EXPECT_FLOAT_EQ(prev(i), static_cast<float>(i));         // original prev
    EXPECT_FLOAT_EQ(prevPrev(i), 10.0f);                     // original prevPrev
  }
}

// ============================================================================
// Elastic Swap and Rotate Tests
// ============================================================================

class SEMSolverDataElasticSwapRotateTest : public ::testing::Test {
 protected:
  void SetUp() override {
    size = 50;
    uxPrev = allocateVector<vectorReal>(size, "uxPrev");
    uxCurr = allocateVector<vectorReal>(size, "uxCurr");
    uyPrev = allocateVector<vectorReal>(size, "uyPrev");
    uyCurr = allocateVector<vectorReal>(size, "uyCurr");
    uzPrev = allocateVector<vectorReal>(size, "uzPrev");
    uzCurr = allocateVector<vectorReal>(size, "uzCurr");
    for (size_t i = 0; i < size; ++i) {
      uxPrev(i) = static_cast<float>(i);
      uxCurr(i) = static_cast<float>(i) * 2.0f;
      uyPrev(i) = static_cast<float>(i) * 3.0f;
      uyCurr(i) = static_cast<float>(i) * 4.0f;
      uzPrev(i) = static_cast<float>(i) * 5.0f;
      uzCurr(i) = static_cast<float>(i) * 6.0f;
    }
  }

  size_t size;
  vectorReal uxPrev, uxCurr;
  vectorReal uyPrev, uyCurr;
  vectorReal uzPrev, uzCurr;
};

TEST_F(SEMSolverDataElasticSwapRotateTest, SwapExchangesPrevAndCurrForAllComponents) {
  WavefieldElastic wavefield(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  solver_data.swapWavefields();

  // Each component: curr <- old prev, prev <- old curr
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(0)(i), static_cast<float>(i));          // old uxPrev
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(0)(i), static_cast<float>(i) * 2.0f);  // old uxCurr
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(1)(i), static_cast<float>(i) * 3.0f);   // old uyPrev
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(1)(i), static_cast<float>(i) * 4.0f);  // old uyCurr
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(2)(i), static_cast<float>(i) * 5.0f);   // old uzPrev
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(2)(i), static_cast<float>(i) * 6.0f);  // old uzCurr
  }
}

TEST_F(SEMSolverDataElasticSwapRotateTest, SwapTwiceRestoresOriginalValues) {
  WavefieldElastic wavefield(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  solver_data.swapWavefields();
  solver_data.swapWavefields();

  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(0)(i), static_cast<float>(i) * 2.0f);   // original uxCurr
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(0)(i), static_cast<float>(i));         // original uxPrev
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(1)(i), static_cast<float>(i) * 4.0f);   // original uyCurr
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(1)(i), static_cast<float>(i) * 3.0f);  // original uyPrev
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(2)(i), static_cast<float>(i) * 6.0f);   // original uzCurr
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(2)(i), static_cast<float>(i) * 5.0f);  // original uzPrev
  }
}

TEST_F(SEMSolverDataElasticSwapRotateTest, RotateRotatesThreeBuffersPerComponent) {
  auto uxPrevPrev = allocateVector<vectorReal>(size, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size, "uzPrevPrev");
  for (size_t i = 0; i < size; ++i) {
    uxPrevPrev(i) = 10.0f;
    uyPrevPrev(i) = 20.0f;
    uzPrevPrev(i) = 30.0f;
  }

  WavefieldElastic wavefield(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  solver_data.rotateWavefields(uxPrevPrev, 0);
  solver_data.rotateWavefields(uyPrevPrev, 1);
  solver_data.rotateWavefields(uzPrevPrev, 2);

  // After rotation for each component:
  //   curr     <- old prevPrev
  //   prev     <- old curr
  //   prevPrev <- old prev
  for (size_t i = 0; i < size; ++i) {
    // ux
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(0)(i), 10.0f);
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(0)(i), static_cast<float>(i) * 2.0f);
    EXPECT_FLOAT_EQ(uxPrevPrev(i), static_cast<float>(i));
    // uy
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(1)(i), 20.0f);
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(1)(i), static_cast<float>(i) * 4.0f);
    EXPECT_FLOAT_EQ(uyPrevPrev(i), static_cast<float>(i) * 3.0f);
    // uz
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(2)(i), 30.0f);
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(2)(i), static_cast<float>(i) * 6.0f);
    EXPECT_FLOAT_EQ(uzPrevPrev(i), static_cast<float>(i) * 5.0f);
  }
}

TEST_F(SEMSolverDataElasticSwapRotateTest, RotateThreeTimesRestoresStatePerComponent) {
  auto uxPrevPrev = allocateVector<vectorReal>(size, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size, "uzPrevPrev");
  for (size_t i = 0; i < size; ++i) {
    uxPrevPrev(i) = 10.0f;
    uyPrevPrev(i) = 20.0f;
    uzPrevPrev(i) = 30.0f;
  }

  WavefieldElastic wavefield(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  for (int r = 0; r < 3; ++r) {
    solver_data.rotateWavefields(uxPrevPrev, 0);
    solver_data.rotateWavefields(uyPrevPrev, 1);
    solver_data.rotateWavefields(uzPrevPrev, 2);
  }

  for (size_t i = 0; i < size; ++i) {
    // ux restored
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(0)(i), static_cast<float>(i) * 2.0f);
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(0)(i), static_cast<float>(i));
    EXPECT_FLOAT_EQ(uxPrevPrev(i), 10.0f);
    // uy restored
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(1)(i), static_cast<float>(i) * 4.0f);
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(1)(i), static_cast<float>(i) * 3.0f);
    EXPECT_FLOAT_EQ(uyPrevPrev(i), 20.0f);
    // uz restored
    EXPECT_FLOAT_EQ(solver_data.getCurrentField(2)(i), static_cast<float>(i) * 6.0f);
    EXPECT_FLOAT_EQ(solver_data.getPreviousField(2)(i), static_cast<float>(i) * 5.0f);
    EXPECT_FLOAT_EQ(uzPrevPrev(i), 30.0f);
  }
}

}  // namespace testing
}  // namespace fe
}  // namespace solver

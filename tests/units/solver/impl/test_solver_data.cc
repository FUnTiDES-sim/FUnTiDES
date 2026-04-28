#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

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

TEST(SEMsolverDataAcousticTest, SwapsWavefieldsWithoutThrowing) {
  WavefieldAcoustic wavefield;
  RhsAcoustic rhs;
  SEMsolverDataAcoustic solver_data(wavefield, rhs);

  EXPECT_NO_THROW(solver_data.swapWavefields());
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

}  // namespace testing
}  // namespace fe
}  // namespace solver

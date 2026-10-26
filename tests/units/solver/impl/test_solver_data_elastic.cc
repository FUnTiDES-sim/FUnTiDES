#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

#include "data_type.h"
#include "device_view_access.h"
#include "rhs_elastic.h"
#include "sem_solver_data.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {
namespace testing {

using test_utils::fillFromHost;
using test_utils::toHost;

// ============================================================================
// Print Method Tests (Using standard cout redirection)
// ============================================================================

class SEMSolverDataElasticPrintTest : public ::testing::Test {
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

TEST_F(SEMSolverDataElasticPrintTest, PrintMethodOutputsElasticInfo) {
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

TEST_F(SEMSolverDataElasticPrintTest, PrintMethodWithThreeBuffers) {
  size_t size = 50;
  auto uxPrevPrev = allocateVector<vectorReal>(size, "uxPrevPrev");
  auto uxPrev = allocateVector<vectorReal>(size, "uxPrev");
  auto uxCurr = allocateVector<vectorReal>(size, "uxCurr");
  auto uyPrevPrev = allocateVector<vectorReal>(size, "uyPrevPrev");
  auto uyPrev = allocateVector<vectorReal>(size, "uyPrev");
  auto uyCurr = allocateVector<vectorReal>(size, "uyCurr");
  auto uzPrevPrev = allocateVector<vectorReal>(size, "uzPrevPrev");
  auto uzPrev = allocateVector<vectorReal>(size, "uzPrev");
  auto uzCurr = allocateVector<vectorReal>(size, "uzCurr");

  WavefieldElastic wavefield(uxPrevPrev, uxPrev, uxCurr, uyPrevPrev, uyPrev, uyCurr, uzPrevPrev, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  EXPECT_NO_THROW(solver_data.print());

  std::string output = buffer.str();

  // Verify all three components are present
  EXPECT_TRUE(output.find("SEMsolverData<") != std::string::npos);
  EXPECT_TRUE(output.find("Field[0]") != std::string::npos);
  EXPECT_TRUE(output.find("Field[1]") != std::string::npos);
  EXPECT_TRUE(output.find("Field[2]") != std::string::npos);
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
    fillFromHost(uxPrev, [&](const auto& h) {
      for (size_t i = 0; i < size; ++i) h(i) = static_cast<float>(i);
    });
    fillFromHost(uxCurr, [&](const auto& h) {
      for (size_t i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 2.0f;
    });
    fillFromHost(uyPrev, [&](const auto& h) {
      for (size_t i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 3.0f;
    });
    fillFromHost(uyCurr, [&](const auto& h) {
      for (size_t i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 4.0f;
    });
    fillFromHost(uzPrev, [&](const auto& h) {
      for (size_t i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 5.0f;
    });
    fillFromHost(uzCurr, [&](const auto& h) {
      for (size_t i = 0; i < size; ++i) h(i) = static_cast<float>(i) * 6.0f;
    });
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
  auto h_curr0 = toHost(solver_data.getCurrentField(0));
  auto h_prev0 = toHost(solver_data.getPreviousField(0));
  auto h_curr1 = toHost(solver_data.getCurrentField(1));
  auto h_prev1 = toHost(solver_data.getPreviousField(1));
  auto h_curr2 = toHost(solver_data.getCurrentField(2));
  auto h_prev2 = toHost(solver_data.getPreviousField(2));
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(h_curr0(i), static_cast<float>(i));         // old uxPrev
    EXPECT_FLOAT_EQ(h_prev0(i), static_cast<float>(i) * 2.0f);  // old uxCurr
    EXPECT_FLOAT_EQ(h_curr1(i), static_cast<float>(i) * 3.0f);  // old uyPrev
    EXPECT_FLOAT_EQ(h_prev1(i), static_cast<float>(i) * 4.0f);  // old uyCurr
    EXPECT_FLOAT_EQ(h_curr2(i), static_cast<float>(i) * 5.0f);  // old uzPrev
    EXPECT_FLOAT_EQ(h_prev2(i), static_cast<float>(i) * 6.0f);  // old uzCurr
  }
}

TEST_F(SEMSolverDataElasticSwapRotateTest, SwapTwiceRestoresOriginalValues) {
  WavefieldElastic wavefield(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  solver_data.swapWavefields();
  solver_data.swapWavefields();

  auto h_curr0 = toHost(solver_data.getCurrentField(0));
  auto h_prev0 = toHost(solver_data.getPreviousField(0));
  auto h_curr1 = toHost(solver_data.getCurrentField(1));
  auto h_prev1 = toHost(solver_data.getPreviousField(1));
  auto h_curr2 = toHost(solver_data.getCurrentField(2));
  auto h_prev2 = toHost(solver_data.getPreviousField(2));
  for (size_t i = 0; i < size; ++i) {
    EXPECT_FLOAT_EQ(h_curr0(i), static_cast<float>(i) * 2.0f);  // original uxCurr
    EXPECT_FLOAT_EQ(h_prev0(i), static_cast<float>(i));         // original uxPrev
    EXPECT_FLOAT_EQ(h_curr1(i), static_cast<float>(i) * 4.0f);  // original uyCurr
    EXPECT_FLOAT_EQ(h_prev1(i), static_cast<float>(i) * 3.0f);  // original uyPrev
    EXPECT_FLOAT_EQ(h_curr2(i), static_cast<float>(i) * 6.0f);  // original uzCurr
    EXPECT_FLOAT_EQ(h_prev2(i), static_cast<float>(i) * 5.0f);  // original uzPrev
  }
}

TEST_F(SEMSolverDataElasticSwapRotateTest, SwapWithThreeBuffersRotatesPerComponent) {
  auto uxPrevPrev = allocateVector<vectorReal>(size, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size, "uzPrevPrev");
  Kokkos::deep_copy(uxPrevPrev, 10.0f);
  Kokkos::deep_copy(uyPrevPrev, 20.0f);
  Kokkos::deep_copy(uzPrevPrev, 30.0f);

  WavefieldElastic wavefield(uxPrevPrev, uxPrev, uxCurr, uyPrevPrev, uyPrev, uyCurr, uzPrevPrev, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  solver_data.swapWavefields();

  // After 3-way rotation for each component:
  //   curr     <- old prevPrev
  //   prev     <- old curr
  //   prevPrev <- old prev
  auto h_curr0 = toHost(solver_data.getCurrentField(0));
  auto h_prev0 = toHost(solver_data.getPreviousField(0));
  auto h_prevprev0 = toHost(solver_data.getPrevPrevField(0));
  auto h_curr1 = toHost(solver_data.getCurrentField(1));
  auto h_prev1 = toHost(solver_data.getPreviousField(1));
  auto h_prevprev1 = toHost(solver_data.getPrevPrevField(1));
  auto h_curr2 = toHost(solver_data.getCurrentField(2));
  auto h_prev2 = toHost(solver_data.getPreviousField(2));
  auto h_prevprev2 = toHost(solver_data.getPrevPrevField(2));
  for (size_t i = 0; i < size; ++i) {
    // ux
    EXPECT_FLOAT_EQ(h_curr0(i), 10.0f);
    EXPECT_FLOAT_EQ(h_prev0(i), static_cast<float>(i) * 2.0f);
    EXPECT_FLOAT_EQ(h_prevprev0(i), static_cast<float>(i));
    // uy
    EXPECT_FLOAT_EQ(h_curr1(i), 20.0f);
    EXPECT_FLOAT_EQ(h_prev1(i), static_cast<float>(i) * 4.0f);
    EXPECT_FLOAT_EQ(h_prevprev1(i), static_cast<float>(i) * 3.0f);
    // uz
    EXPECT_FLOAT_EQ(h_curr2(i), 30.0f);
    EXPECT_FLOAT_EQ(h_prev2(i), static_cast<float>(i) * 6.0f);
    EXPECT_FLOAT_EQ(h_prevprev2(i), static_cast<float>(i) * 5.0f);
  }
}

TEST_F(SEMSolverDataElasticSwapRotateTest, SwapThreeTimesRestoresStatePerComponentWithThreeBuffers) {
  auto uxPrevPrev = allocateVector<vectorReal>(size, "uxPrevPrev");
  auto uyPrevPrev = allocateVector<vectorReal>(size, "uyPrevPrev");
  auto uzPrevPrev = allocateVector<vectorReal>(size, "uzPrevPrev");
  Kokkos::deep_copy(uxPrevPrev, 10.0f);
  Kokkos::deep_copy(uyPrevPrev, 20.0f);
  Kokkos::deep_copy(uzPrevPrev, 30.0f);

  WavefieldElastic wavefield(uxPrevPrev, uxPrev, uxCurr, uyPrevPrev, uyPrev, uyCurr, uzPrevPrev, uzPrev, uzCurr);
  RhsElastic rhs;
  SEMsolverDataElastic solver_data(wavefield, rhs, false);

  for (int r = 0; r < 3; ++r) {
    solver_data.swapWavefields();
  }

  auto h_curr0 = toHost(solver_data.getCurrentField(0));
  auto h_prev0 = toHost(solver_data.getPreviousField(0));
  auto h_prevprev0 = toHost(solver_data.getPrevPrevField(0));
  auto h_curr1 = toHost(solver_data.getCurrentField(1));
  auto h_prev1 = toHost(solver_data.getPreviousField(1));
  auto h_prevprev1 = toHost(solver_data.getPrevPrevField(1));
  auto h_curr2 = toHost(solver_data.getCurrentField(2));
  auto h_prev2 = toHost(solver_data.getPreviousField(2));
  auto h_prevprev2 = toHost(solver_data.getPrevPrevField(2));
  for (size_t i = 0; i < size; ++i) {
    // ux restored
    EXPECT_FLOAT_EQ(h_curr0(i), static_cast<float>(i) * 2.0f);
    EXPECT_FLOAT_EQ(h_prev0(i), static_cast<float>(i));
    EXPECT_FLOAT_EQ(h_prevprev0(i), 10.0f);
    // uy restored
    EXPECT_FLOAT_EQ(h_curr1(i), static_cast<float>(i) * 4.0f);
    EXPECT_FLOAT_EQ(h_prev1(i), static_cast<float>(i) * 3.0f);
    EXPECT_FLOAT_EQ(h_prevprev1(i), 20.0f);
    // uz restored
    EXPECT_FLOAT_EQ(h_curr2(i), static_cast<float>(i) * 6.0f);
    EXPECT_FLOAT_EQ(h_prev2(i), static_cast<float>(i) * 5.0f);
    EXPECT_FLOAT_EQ(h_prevprev2(i), 30.0f);
  }
}

}  // namespace testing
}  // namespace fe
}  // namespace solver

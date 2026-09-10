/**
 * @file test_dg_padaptive_solver_data.cc
 * @brief Unit tests for DGPAdaptiveSolverData.
 *
 * Covers constructors, accessors, swapWavefields, and the print helper.
 * No mesh or model is needed — only Kokkos views allocated with the
 * project's allocateArray2D / allocateVector utilities.
 */
#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "dg_padaptive_solver_data.h"

namespace solver {
namespace fe {
namespace test {

class DGPAdaptiveSolverDataTest : public ::testing::Test {
 protected:
  void SetUp() override {
    nElem_ = 4;
    nDof_min_ = 8;
    nDof_max_ = 27;
    nNode_ = 18;
    nSrc_ = 2;
    nSample_ = 10;

    p_min_Curr = allocateArray2D<arrayReal>(nElem_, nDof_min_, "pnPMinDGCurr");
    p_min_Prev = allocateArray2D<arrayReal>(nElem_, nDof_min_, "pnPMinDGPrev");
    p_max_Curr = allocateArray2D<arrayReal>(nElem_, nDof_max_, "pnPMaxDGCurr");
    p_max_Prev = allocateArray2D<arrayReal>(nElem_, nDof_max_, "pnPMaxDGPrev");
    rhsTerm_pmin = allocateArray2D<arrayReal>(nSrc_, nSample_, "rhsTermPMin");
    rhsTerm_pmax = allocateArray2D<arrayReal>(nSrc_, nSample_, "rhsTermPMax");
    rhsElem_ = allocateVector<vectorInt>(nSrc_, "rhsElem");
    rhsWeights_pmin = allocateArray2D<arrayReal>(nSrc_, nDof_min_, "rhsWeightsPMin");
    rhsWeights_pmax = allocateArray2D<arrayReal>(nSrc_, nDof_max_, "rhsWeightsPMax");

    for (int e = 0; e < nElem_; ++e) {
      for (int d = 0; d < nDof_min_; ++d) {
        p_min_Prev(e, d) = static_cast<float>(e * nDof_min_ + d);
        p_min_Curr(e, d) = static_cast<float>(e * nDof_min_ + d) * 2.0f;
      }
      for (int d = 0; d < nDof_max_; ++d) {
        p_max_Prev(e, d) = static_cast<float>(e * nDof_max_ + d);
        p_max_Curr(e, d) = static_cast<float>(e * nDof_max_ + d) * 2.0f;
      }
    }

    for (int s = 0; s < nSrc_; ++s) {
      rhsElem_(s) = s * 3;
      for (int t = 0; t < nSample_; ++t) {
        rhsTerm_pmin(s, t) = static_cast<float>(s * 10 + t);
        rhsTerm_pmax(s, t) = static_cast<float>(s * 10 + t);
      }
      for (int d = 0; d < nDof_min_; ++d) rhsWeights_pmin(s, d) = static_cast<float>(d) * 0.1f;
      for (int d = 0; d < nDof_max_; ++d) rhsWeights_pmax(s, d) = static_cast<float>(d) * 0.1f;
    }
  }

  int nElem_, nDof_min_, nDof_max_, nNode_, nSrc_, nSample_;
  arrayReal p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr, rhsTerm_pmin, rhsTerm_pmax, rhsWeights_pmin,
      rhsWeights_pmax;
  vectorInt rhsElem_;
};

// ============================================================
// Constructor
// ============================================================

TEST_F(DGPAdaptiveSolverDataTest, Constructor_StoresAllExtents) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  EXPECT_EQ(data.m_wavefield.getPMinPreviousField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.m_wavefield.getPMinPreviousField(0).extent(1), static_cast<size_t>(nDof_min_));
  EXPECT_EQ(data.m_wavefield.getPMinCurrentField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.m_wavefield.getPMinCurrentField(0).extent(1), static_cast<size_t>(nDof_min_));
  EXPECT_EQ(data.m_wavefield.getPMaxPreviousField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.m_wavefield.getPMaxPreviousField(0).extent(1), static_cast<size_t>(nDof_max_));
  EXPECT_EQ(data.m_wavefield.getPMaxCurrentField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.m_wavefield.getPMaxCurrentField(0).extent(1), static_cast<size_t>(nDof_max_));
  // PMin term
  EXPECT_EQ(data.m_rhs.getTerm(0).extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getTerm(0).extent(1), static_cast<size_t>(nSample_));
  // PMax term
  EXPECT_EQ(data.m_rhs.getTerm(1).extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getTerm(1).extent(1), static_cast<size_t>(nSample_));

  EXPECT_EQ(data.m_rhs.getElement().extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getWeights(0).extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getWeights(0).extent(1), static_cast<size_t>(nDof_min_));
  EXPECT_EQ(data.m_rhs.getWeights(1).extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getWeights(1).extent(1), static_cast<size_t>(nDof_max_));
}

TEST_F(DGPAdaptiveSolverDataTest, Constructor_DataValuesPreserved) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  for (int e = 0; e < nElem_; ++e) {
    for (int d = 0; d < nDof_min_; ++d) {
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMinPreviousField(0)(e, d), static_cast<float>(e * nDof_min_ + d));
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMinCurrentField(0)(e, d), static_cast<float>(e * nDof_min_ + d) * 2.0f);
    }
    for (int d = 0; d < nDof_max_; ++d) {
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxPreviousField(0)(e, d), static_cast<float>(e * nDof_max_ + d));
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxCurrentField(0)(e, d), static_cast<float>(e * nDof_max_ + d) * 2.0f);
    }
  }

  for (int s = 0; s < nSrc_; ++s) {
    EXPECT_EQ(data.m_rhs.getElement()(s), s * 3);
    for (int t = 0; t < nSample_; ++t) {
      EXPECT_FLOAT_EQ(data.m_rhs.getTerm(0)(s, t), static_cast<float>(s * 10 + t));  // PMin term
      EXPECT_FLOAT_EQ(data.m_rhs.getTerm(1)(s, t), static_cast<float>(s * 10 + t));  // PMax term
    }
    for (int d = 0; d < nDof_min_; ++d) EXPECT_FLOAT_EQ(data.m_rhs.getWeights(0)(s, d), static_cast<float>(d) * 0.1f);
    for (int d = 0; d < nDof_max_; ++d) EXPECT_FLOAT_EQ(data.m_rhs.getWeights(1)(s, d), static_cast<float>(d) * 0.1f);
  }
}

TEST_F(DGPAdaptiveSolverDataTest, IsDistributed_DefaultsFalse) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  EXPECT_FALSE(data.isDistributed);
}

// ============================================================
// Accessors
// ============================================================

TEST_F(DGPAdaptiveSolverDataTest, GetCurrentField_ReturnsPnCurr) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  auto pmin_curr = data.m_wavefield.getPMinCurrentField(0);
  auto pmax_curr = data.m_wavefield.getPMaxCurrentField(0);

  ASSERT_EQ(pmin_curr.extent(0), static_cast<size_t>(nElem_));
  ASSERT_EQ(pmin_curr.extent(1), static_cast<size_t>(nDof_min_));
  ASSERT_EQ(pmax_curr.extent(0), static_cast<size_t>(nElem_));
  ASSERT_EQ(pmax_curr.extent(1), static_cast<size_t>(nDof_max_));
  for (int e = 0; e < nElem_; ++e) {
    for (int d = 0; d < nDof_min_; ++d) EXPECT_FLOAT_EQ(pmin_curr(e, d), static_cast<float>(e * nDof_min_ + d) * 2.0f);
    for (int d = 0; d < nDof_max_; ++d) EXPECT_FLOAT_EQ(pmax_curr(e, d), static_cast<float>(e * nDof_max_ + d) * 2.0f);
  }
}

TEST_F(DGPAdaptiveSolverDataTest, GetPreviousField_ReturnsPnPrev) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  auto pmin_prev = data.m_wavefield.getPMinPreviousField(0);
  auto pmax_prev = data.m_wavefield.getPMaxPreviousField(0);

  ASSERT_EQ(pmin_prev.extent(0), static_cast<size_t>(nElem_));
  ASSERT_EQ(pmin_prev.extent(1), static_cast<size_t>(nDof_min_));
  ASSERT_EQ(pmax_prev.extent(0), static_cast<size_t>(nElem_));
  ASSERT_EQ(pmax_prev.extent(1), static_cast<size_t>(nDof_max_));
  for (int e = 0; e < nElem_; ++e) {
    for (int d = 0; d < nDof_min_; ++d) EXPECT_FLOAT_EQ(pmin_prev(e, d), static_cast<float>(e * nDof_min_ + d));
    for (int d = 0; d < nDof_max_; ++d) EXPECT_FLOAT_EQ(pmax_prev(e, d), static_cast<float>(e * nDof_max_ + d));
  }
}

TEST_F(DGPAdaptiveSolverDataTest, GetRhsTerm_ReturnsTerm) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  auto pmin_term = data.m_rhs.getTerm(0);  // PMin term
  auto pmax_term = data.m_rhs.getTerm(1);  // PMax term

  ASSERT_EQ(pmin_term.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(pmin_term.extent(1), static_cast<size_t>(nSample_));
  ASSERT_EQ(pmax_term.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(pmax_term.extent(1), static_cast<size_t>(nSample_));
  for (int s = 0; s < nSrc_; ++s)
    for (int t = 0; t < nSample_; ++t) {
      EXPECT_FLOAT_EQ(pmin_term(s, t), static_cast<float>(s * 10 + t));
      EXPECT_FLOAT_EQ(pmax_term(s, t), static_cast<float>(s * 10 + t));
    }
}

TEST_F(DGPAdaptiveSolverDataTest, GetRhsElement_ReturnsElem) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  auto elem = data.m_rhs.getElement();

  ASSERT_EQ(elem.extent(0), static_cast<size_t>(nSrc_));
  for (int s = 0; s < nSrc_; ++s) EXPECT_EQ(elem(s), s * 3);
}

TEST_F(DGPAdaptiveSolverDataTest, GetRhsWeights_ReturnsWeights) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  auto w_pmin = data.m_rhs.getWeights(0);
  auto w_pmax = data.m_rhs.getWeights(1);

  ASSERT_EQ(w_pmin.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(w_pmin.extent(1), static_cast<size_t>(nDof_min_));
  ASSERT_EQ(w_pmax.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(w_pmax.extent(1), static_cast<size_t>(nDof_max_));
  for (int s = 0; s < nSrc_; ++s) {
    for (int d = 0; d < nDof_min_; ++d) EXPECT_FLOAT_EQ(w_pmin(s, d), static_cast<float>(d) * 0.1f);
    for (int d = 0; d < nDof_max_; ++d) EXPECT_FLOAT_EQ(w_pmax(s, d), static_cast<float>(d) * 0.1f);
  }
}

// ============================================================
// swapWavefields
// ============================================================

TEST_F(DGPAdaptiveSolverDataTest, SwapWavefields_ExchangesPrevAndCurr) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  float const pmin_prev00 = data.m_wavefield.getPMinPreviousField(0)(0, 0);
  float const pmin_curr00 = data.m_wavefield.getPMinCurrentField(0)(0, 0);
  float const pmax_prev00 = data.m_wavefield.getPMaxPreviousField(0)(0, 0);
  float const pmax_curr00 = data.m_wavefield.getPMaxCurrentField(0)(0, 0);

  data.swapWavefields();

  EXPECT_FLOAT_EQ(data.m_wavefield.getPMinPreviousField(0)(0, 0), pmin_curr00);
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMinCurrentField(0)(0, 0), pmin_prev00);
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxPreviousField(0)(0, 0), pmax_curr00);
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxCurrentField(0)(0, 0), pmax_prev00);
}

TEST_F(DGPAdaptiveSolverDataTest, SwapWavefields_TwiceRestoresOriginal) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  float const pmin_prev00 = data.m_wavefield.getPMinPreviousField(0)(0, 0);
  float const pmin_curr00 = data.m_wavefield.getPMinCurrentField(0)(0, 0);
  float const pmax_prev00 = data.m_wavefield.getPMaxPreviousField(0)(0, 0);
  float const pmax_curr00 = data.m_wavefield.getPMaxCurrentField(0)(0, 0);

  data.swapWavefields();
  data.swapWavefields();

  EXPECT_FLOAT_EQ(data.m_wavefield.getPMinPreviousField(0)(0, 0), pmin_prev00);
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMinCurrentField(0)(0, 0), pmin_curr00);
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxPreviousField(0)(0, 0), pmax_prev00);
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxCurrentField(0)(0, 0), pmax_curr00);
}

TEST_F(DGPAdaptiveSolverDataTest, SwapWavefields_IsShallowHandleSwap) {
  // After swap: pnPrev points to pCurr_ allocation and vice-versa.
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  data.swapWavefields();

  // pnPrev now references the old pCurr_ allocation
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMinPreviousField(0)(0, 0), p_min_Curr(0, 0));
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxPreviousField(0)(0, 0), p_max_Curr(0, 0));
  // pnCurr now references the old pPrev_ allocation
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMinCurrentField(0)(0, 0), p_min_Prev(0, 0));
  EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxCurrentField(0)(0, 0), p_max_Prev(0, 0));
}

TEST_F(DGPAdaptiveSolverDataTest, SwapWavefields_AllElementsSwapped) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  data.swapWavefields();

  for (int e = 0; e < nElem_; ++e) {
    for (int d = 0; d < nDof_min_; ++d) {
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMinPreviousField(0)(e, d), p_min_Curr(e, d));
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMinCurrentField(0)(e, d), p_min_Prev(e, d));
    }
    for (int d = 0; d < nDof_max_; ++d) {
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxPreviousField(0)(e, d), p_max_Curr(e, d));
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxCurrentField(0)(e, d), p_max_Prev(e, d));
    }
  }
}

// ============================================================
// Misc
// ============================================================

TEST_F(DGPAdaptiveSolverDataTest, Print_DoesNotCrash) {
  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem_, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  EXPECT_NO_THROW(data.print());
}

TEST_F(DGPAdaptiveSolverDataTest, EmptyViews_ValidConstruction) {
  auto emptyPrev_pmin = allocateArray2D<arrayReal>(0, 0, "epmin");
  auto emptyCurr_pmin = allocateArray2D<arrayReal>(0, 0, "ecmin");
  auto emptyPrev_pmax = allocateArray2D<arrayReal>(0, 0, "epmax");
  auto emptyCurr_pmax = allocateArray2D<arrayReal>(0, 0, "ecmax");
  auto emptyTerm_pmin = allocateArray2D<arrayReal>(0, 0, "etmin");
  auto emptyTerm_pmax = allocateArray2D<arrayReal>(0, 0, "etmax");
  auto emptyElem = allocateVector<vectorInt>(0, "ee");
  auto emptyW_pmin = allocateArray2D<arrayReal>(0, 0, "ewmin");
  auto emptyW_pmax = allocateArray2D<arrayReal>(0, 0, "ewmax");

  DGPAdaptiveWavefieldAcoustic empty_wavefield(emptyPrev_pmin, emptyCurr_pmin, emptyPrev_pmax, emptyCurr_pmax);
  DGPAdaptiveRhsAcoustic empty_rhs(emptyTerm_pmin, emptyTerm_pmax, emptyElem, emptyW_pmin, emptyW_pmax);

  DGPAdaptiveSolverData data(empty_wavefield, empty_rhs);

  EXPECT_EQ(data.m_wavefield.getPMinPreviousField(0).extent(0), 0u);
  EXPECT_EQ(data.m_wavefield.getPMinCurrentField(0).extent(0), 0u);
  EXPECT_EQ(data.m_wavefield.getPMaxPreviousField(0).extent(0), 0u);
  EXPECT_EQ(data.m_wavefield.getPMaxCurrentField(0).extent(0), 0u);
  EXPECT_EQ(data.m_rhs.getTerm(0).extent(0), 0u);  // PMin term
  EXPECT_EQ(data.m_rhs.getTerm(1).extent(0), 0u);  // PMax term
  EXPECT_EQ(data.m_rhs.getElement().extent(0), 0u);
  EXPECT_FALSE(data.isDistributed);
  EXPECT_NO_THROW(data.swapWavefields());
}

}  // namespace test
}  // namespace fe
}  // namespace solver

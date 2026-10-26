/**
 * @file test_dg-sem_solver_data.cc
 * @brief Unit tests for DGSEMsolverData.
 *
 * Covers constructors, accessors, swapWavefields, and the print helper.
 * No mesh or model is needed — only Kokkos views allocated with the
 * project's allocateArray2D / allocateVector utilities.
 */
#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "device_view_access.h"
#include "dg-sem_solver_data.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromHost;
using test_utils::getAt;
using test_utils::toHost;

class DGSEMsolverDataTest : public ::testing::Test {
 protected:
  void SetUp() override {
    nElem_ = 4;
    nDof_ = 8;
    nNode_ = 18;
    nSrc_ = 2;
    nSample_ = 10;

    p_dg_prev = allocateArray2D<arrayReal>(nElem_, nDof_, "pPrevDG");
    p_dg_curr = allocateArray2D<arrayReal>(nElem_, nDof_, "pCurrDG");
    p_sem_prev = allocateVector<vectorReal>(nNode_, "pPrevSEM");
    p_sem_curr = allocateVector<vectorReal>(nNode_, "pCurrSEM");
    rhsTerm_dg = allocateArray2D<arrayReal>(nSrc_, nSample_, "rhsTermDG");
    rhsTerm_sem = allocateArray2D<arrayReal>(nSrc_, nSample_, "rhsTermSEM");
    rhsElem_ = allocateVector<vectorInt>(nSrc_, "rhsElem");
    rhsWeights_ = allocateArray2D<arrayReal>(nSrc_, nDof_, "rhsWeights");

    fillFromHost(p_dg_prev, [&](const auto& h) {
      for (int e = 0; e < nElem_; ++e)
        for (int d = 0; d < nDof_; ++d) h(e, d) = static_cast<float>(e * nDof_ + d);
    });
    fillFromHost(p_dg_curr, [&](const auto& h) {
      for (int e = 0; e < nElem_; ++e)
        for (int d = 0; d < nDof_; ++d) h(e, d) = static_cast<float>(e * nDof_ + d) * 2.0f;
    });
    fillFromHost(p_sem_prev, [&](const auto& h) {
      for (int n = 0; n < nNode_; ++n) h(n) = static_cast<float>(n);
    });
    fillFromHost(p_sem_curr, [&](const auto& h) {
      for (int n = 0; n < nNode_; ++n) h(n) = static_cast<float>(n) * 2.0f;
    });

    fillFromHost(rhsElem_, [&](const auto& h) {
      for (int s = 0; s < nSrc_; ++s) h(s) = s * 3;
    });
    fillFromHost(rhsTerm_dg, [&](const auto& h) {
      for (int s = 0; s < nSrc_; ++s)
        for (int t = 0; t < nSample_; ++t) h(s, t) = static_cast<float>(s * 10 + t);
    });
    Kokkos::deep_copy(rhsTerm_sem, rhsTerm_dg);
    fillFromHost(rhsWeights_, [&](const auto& h) {
      for (int s = 0; s < nSrc_; ++s)
        for (int d = 0; d < nDof_; ++d) h(s, d) = static_cast<float>(d) * 0.1f;
    });
  }

  int nElem_, nDof_, nNode_, nSrc_, nSample_;
  arrayReal p_dg_prev, p_dg_curr, rhsTerm_dg, rhsTerm_sem, rhsWeights_;
  vectorReal p_sem_prev, p_sem_curr;
  vectorInt rhsElem_;
};

// ============================================================
// Constructor
// ============================================================

TEST_F(DGSEMsolverDataTest, Constructor_StoresAllExtents) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  EXPECT_EQ(data.m_wavefield.getDGPreviousField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.m_wavefield.getDGPreviousField(0).extent(1), static_cast<size_t>(nDof_));
  EXPECT_EQ(data.m_wavefield.getDGCurrentField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.m_wavefield.getDGCurrentField(0).extent(1), static_cast<size_t>(nDof_));
  EXPECT_EQ(data.m_wavefield.getSEMPreviousField(0).extent(0), static_cast<size_t>(nNode_));
  EXPECT_EQ(data.m_wavefield.getSEMCurrentField(0).extent(0), static_cast<size_t>(nNode_));
  // DG term
  EXPECT_EQ(data.m_rhs.getTerm(0).extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getTerm(0).extent(1), static_cast<size_t>(nSample_));
  // SEM term
  EXPECT_EQ(data.m_rhs.getTerm(1).extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getTerm(1).extent(1), static_cast<size_t>(nSample_));

  EXPECT_EQ(data.m_rhs.getElement().extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getWeights().extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.m_rhs.getWeights().extent(1), static_cast<size_t>(nDof_));
}

TEST_F(DGSEMsolverDataTest, Constructor_DataValuesPreserved) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  auto h_dg_prev0 = toHost(data.m_wavefield.getDGPreviousField(0));
  auto h_dg_curr0 = toHost(data.m_wavefield.getDGCurrentField(0));
  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) {
      EXPECT_FLOAT_EQ(h_dg_prev0(e, d), static_cast<float>(e * nDof_ + d));
      EXPECT_FLOAT_EQ(h_dg_curr0(e, d), static_cast<float>(e * nDof_ + d) * 2.0f);
    }
  auto h_sem_prev0 = toHost(data.m_wavefield.getSEMPreviousField(0));
  auto h_sem_curr0 = toHost(data.m_wavefield.getSEMCurrentField(0));
  for (int n = 0; n < nNode_; ++n) {
    EXPECT_FLOAT_EQ(h_sem_prev0(n), static_cast<float>(n));
    EXPECT_FLOAT_EQ(h_sem_curr0(n), static_cast<float>(n) * 2.0f);
  }

  auto h_element = toHost(data.m_rhs.getElement());
  auto h_term0 = toHost(data.m_rhs.getTerm(0));
  auto h_term1 = toHost(data.m_rhs.getTerm(1));
  auto h_weights = toHost(data.m_rhs.getWeights());
  for (int s = 0; s < nSrc_; ++s) {
    EXPECT_EQ(h_element(s), s * 3);
    for (int t = 0; t < nSample_; ++t) {
      EXPECT_FLOAT_EQ(h_term0(s, t), static_cast<float>(s * 10 + t));  // DG term
      EXPECT_FLOAT_EQ(h_term1(s, t), static_cast<float>(s * 10 + t));  // SEM term
    }
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(h_weights(s, d), static_cast<float>(d) * 0.1f);
  }
}

TEST_F(DGSEMsolverDataTest, IsDistributed_DefaultsFalse) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  EXPECT_FALSE(data.isDistributed);
}

// ============================================================
// Accessors
// ============================================================

TEST_F(DGSEMsolverDataTest, GetCurrentField_ReturnsPnCurr) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  auto dg_curr = toHost(data.m_wavefield.getDGCurrentField(0));
  auto sem_curr = toHost(data.m_wavefield.getSEMCurrentField(0));

  ASSERT_EQ(dg_curr.extent(0), static_cast<size_t>(nElem_));
  ASSERT_EQ(dg_curr.extent(1), static_cast<size_t>(nDof_));
  ASSERT_EQ(sem_curr.extent(0), static_cast<size_t>(nNode_));
  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(dg_curr(e, d), static_cast<float>(e * nDof_ + d) * 2.0f);
  for (int n = 0; n < nNode_; ++n) EXPECT_FLOAT_EQ(sem_curr(n), static_cast<float>(n) * 2.0f);
}

TEST_F(DGSEMsolverDataTest, GetPreviousField_ReturnsPnPrev) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  auto dg_prev = toHost(data.m_wavefield.getDGPreviousField(0));
  auto sem_prev = toHost(data.m_wavefield.getSEMPreviousField(0));

  ASSERT_EQ(dg_prev.extent(0), static_cast<size_t>(nElem_));
  ASSERT_EQ(dg_prev.extent(1), static_cast<size_t>(nDof_));
  ASSERT_EQ(sem_prev.extent(0), static_cast<size_t>(nNode_));
  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(dg_prev(e, d), static_cast<float>(e * nDof_ + d));
  for (int n = 0; n < nNode_; ++n) EXPECT_FLOAT_EQ(sem_prev(n), static_cast<float>(n));
}

TEST_F(DGSEMsolverDataTest, GetRhsTerm_ReturnsTerm) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  auto dg_term = toHost(data.m_rhs.getTerm(0));   // DG term
  auto sem_term = toHost(data.m_rhs.getTerm(1));  // SEM term

  ASSERT_EQ(dg_term.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(dg_term.extent(1), static_cast<size_t>(nSample_));
  ASSERT_EQ(sem_term.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(sem_term.extent(1), static_cast<size_t>(nSample_));
  for (int s = 0; s < nSrc_; ++s)
    for (int t = 0; t < nSample_; ++t) {
      EXPECT_FLOAT_EQ(dg_term(s, t), static_cast<float>(s * 10 + t));
      EXPECT_FLOAT_EQ(sem_term(s, t), static_cast<float>(s * 10 + t));
    }
}

TEST_F(DGSEMsolverDataTest, GetRhsElement_ReturnsElem) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  auto elem = toHost(data.m_rhs.getElement());

  ASSERT_EQ(elem.extent(0), static_cast<size_t>(nSrc_));
  for (int s = 0; s < nSrc_; ++s) EXPECT_EQ(elem(s), s * 3);
}

TEST_F(DGSEMsolverDataTest, GetRhsWeights_ReturnsWeights) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  auto w = toHost(data.m_rhs.getWeights());

  ASSERT_EQ(w.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(w.extent(1), static_cast<size_t>(nDof_));
  for (int s = 0; s < nSrc_; ++s)
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(w(s, d), static_cast<float>(d) * 0.1f);
}

// ============================================================
// swapWavefields
// ============================================================

TEST_F(DGSEMsolverDataTest, SwapWavefields_ExchangesPrevAndCurr) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  float const dg_prev00 = getAt(data.m_wavefield.getDGPreviousField(0), 0, 0);
  float const dg_curr00 = getAt(data.m_wavefield.getDGCurrentField(0), 0, 0);
  float const sem_prev00 = getAt(data.m_wavefield.getSEMPreviousField(0), 0);
  float const sem_curr00 = getAt(data.m_wavefield.getSEMCurrentField(0), 0);

  data.swapWavefields();

  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getDGPreviousField(0), 0, 0), dg_curr00);
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getDGCurrentField(0), 0, 0), dg_prev00);
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getSEMPreviousField(0), 0), sem_curr00);
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getSEMCurrentField(0), 0), sem_prev00);
}

TEST_F(DGSEMsolverDataTest, SwapWavefields_TwiceRestoresOriginal) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  float const dg_prev00 = getAt(data.m_wavefield.getDGPreviousField(0), 0, 0);
  float const dg_curr00 = getAt(data.m_wavefield.getDGCurrentField(0), 0, 0);
  float const sem_prev00 = getAt(data.m_wavefield.getSEMPreviousField(0), 0);
  float const sem_curr00 = getAt(data.m_wavefield.getSEMCurrentField(0), 0);

  data.swapWavefields();
  data.swapWavefields();

  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getDGPreviousField(0), 0, 0), dg_prev00);
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getDGCurrentField(0), 0, 0), dg_curr00);
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getSEMPreviousField(0), 0), sem_prev00);
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getSEMCurrentField(0), 0), sem_curr00);
}

TEST_F(DGSEMsolverDataTest, SwapWavefields_IsShallowHandleSwap) {
  // After swap: pnPrev points to pCurr_ allocation and vice-versa.
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  data.swapWavefields();

  // pnPrev now references the old pCurr_ allocation
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getDGPreviousField(0), 0, 0), getAt(p_dg_curr, 0, 0));
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getSEMPreviousField(0), 0), getAt(p_sem_curr, 0));
  // pnCurr now references the old pPrev_ allocation
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getDGCurrentField(0), 0, 0), getAt(p_dg_prev, 0, 0));
  EXPECT_FLOAT_EQ(getAt(data.m_wavefield.getSEMCurrentField(0), 0), getAt(p_sem_prev, 0));
}

TEST_F(DGSEMsolverDataTest, SwapWavefields_AllElementsSwapped) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  data.swapWavefields();

  auto h_dg_prev0 = toHost(data.m_wavefield.getDGPreviousField(0));
  auto h_p_dg_curr = toHost(p_dg_curr);
  auto h_dg_curr0 = toHost(data.m_wavefield.getDGCurrentField(0));
  auto h_p_dg_prev = toHost(p_dg_prev);
  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) {
      EXPECT_FLOAT_EQ(h_dg_prev0(e, d), h_p_dg_curr(e, d));
      EXPECT_FLOAT_EQ(h_dg_curr0(e, d), h_p_dg_prev(e, d));
    }
  auto h_sem_prev0 = toHost(data.m_wavefield.getSEMPreviousField(0));
  auto h_p_sem_curr = toHost(p_sem_curr);
  auto h_sem_curr0 = toHost(data.m_wavefield.getSEMCurrentField(0));
  auto h_p_sem_prev = toHost(p_sem_prev);
  for (int n = 0; n < nNode_; ++n) {
    EXPECT_FLOAT_EQ(h_sem_prev0(n), h_p_sem_curr(n));
    EXPECT_FLOAT_EQ(h_sem_curr0(n), h_p_sem_prev(n));
  }
}

// ============================================================
// Misc
// ============================================================

TEST_F(DGSEMsolverDataTest, Print_DoesNotCrash) {
  DGSEMWavefieldAcoustic wavefield(p_dg_prev, p_dg_curr, p_sem_prev, p_sem_curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem_, rhsWeights_);

  DGSEMsolverData data(wavefield, rhs);

  EXPECT_NO_THROW(data.print());
}

TEST_F(DGSEMsolverDataTest, EmptyViews_ValidConstruction) {
  auto emptyPrev_dg = allocateArray2D<arrayReal>(0, 0, "epdg");
  auto emptyCurr_dg = allocateArray2D<arrayReal>(0, 0, "ec");
  auto emptyPrev_sem = allocateVector<vectorReal>(0, "epsem");
  auto emptyCurr_sem = allocateVector<vectorReal>(0, "ecsel");
  auto emptyTerm_dg = allocateArray2D<arrayReal>(0, 0, "etdg");
  auto emptyTerm_sem = allocateArray2D<arrayReal>(0, 0, "etsem");
  auto emptyElem = allocateVector<vectorInt>(0, "ee");
  auto emptyW = allocateArray2D<arrayReal>(0, 0, "ew");

  DGSEMWavefieldAcoustic empty_wavefield(emptyPrev_dg, emptyCurr_dg, emptyPrev_sem, emptyCurr_sem);
  DGSEMRhsAcoustic empty_rhs(emptyTerm_dg, emptyTerm_sem, emptyElem, emptyW);

  DGSEMsolverData data(empty_wavefield, empty_rhs);

  EXPECT_EQ(data.m_wavefield.getDGPreviousField(0).extent(0), 0u);
  EXPECT_EQ(data.m_wavefield.getDGCurrentField(0).extent(0), 0u);
  EXPECT_EQ(data.m_wavefield.getSEMPreviousField(0).extent(0), 0u);
  EXPECT_EQ(data.m_wavefield.getSEMCurrentField(0).extent(0), 0u);
  EXPECT_EQ(data.m_rhs.getTerm(0).extent(0), 0u);  // DG term
  EXPECT_EQ(data.m_rhs.getTerm(1).extent(0), 0u);  // SEM term
  EXPECT_EQ(data.m_rhs.getElement().extent(0), 0u);
  EXPECT_FALSE(data.isDistributed);
  EXPECT_NO_THROW(data.swapWavefields());
}

}  // namespace test
}  // namespace fe
}  // namespace solver

/**
 * @file test_dg_solver_data.cc
 * @brief Unit tests for DGsolverDataAcoustic.
 *
 * Covers constructors, accessors, swapWavefields, and the print helper.
 * No mesh or model is needed — only Kokkos views allocated with the
 * project's allocateArray2D / allocateVector utilities.
 */
#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "dg_solver_data.h"

namespace solver {
namespace fe {
namespace test {

class DGsolverDataAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    nElem_ = 4;
    nDof_ = 8;
    nSrc_ = 2;
    nSample_ = 10;

    pPrev_ = allocateArray2D<ARRAY_REAL_VIEW>(nElem_, nDof_, "pPrev");
    pCurr_ = allocateArray2D<ARRAY_REAL_VIEW>(nElem_, nDof_, "pCurr");
    rhsTerm_ = allocateArray2D<ARRAY_REAL_VIEW>(nSrc_, nSample_, "rhsTerm");
    rhsElem_ = allocateVector<VECTOR_INT_VIEW>(nSrc_, "rhsElem");
    rhsWeights_ = allocateArray2D<ARRAY_REAL_VIEW>(nSrc_, nDof_, "rhsWeights");

    for (int e = 0; e < nElem_; ++e)
      for (int d = 0; d < nDof_; ++d) {
        pPrev_(e, d) = static_cast<float>(e * nDof_ + d);
        pCurr_(e, d) = static_cast<float>(e * nDof_ + d) * 2.0f;
      }

    for (int s = 0; s < nSrc_; ++s) {
      rhsElem_(s) = s * 3;
      for (int t = 0; t < nSample_; ++t) rhsTerm_(s, t) = static_cast<float>(s * 10 + t);
      for (int d = 0; d < nDof_; ++d) rhsWeights_(s, d) = static_cast<float>(d) * 0.1f;
    }
  }

  int nElem_, nDof_, nSrc_, nSample_;
  ARRAY_REAL_VIEW pPrev_, pCurr_, rhsTerm_, rhsWeights_;
  VECTOR_INT_VIEW rhsElem_;
};

// ============================================================
// Constructor
// ============================================================

TEST_F(DGsolverDataAcousticTest, Constructor_StoresAllExtents) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  EXPECT_EQ(data.getPreviousField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.getPreviousField(0).extent(1), static_cast<size_t>(nDof_));
  EXPECT_EQ(data.getCurrentField(0).extent(0), static_cast<size_t>(nElem_));
  EXPECT_EQ(data.getCurrentField(0).extent(1), static_cast<size_t>(nDof_));
  EXPECT_EQ(data.getRhsTerm(0).extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.getRhsTerm(0).extent(1), static_cast<size_t>(nSample_));
  EXPECT_EQ(data.getRhsElement().extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.getRhsWeights().extent(0), static_cast<size_t>(nSrc_));
  EXPECT_EQ(data.getRhsWeights().extent(1), static_cast<size_t>(nDof_));
}

TEST_F(DGsolverDataAcousticTest, Constructor_DataValuesPreserved) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) {
      EXPECT_FLOAT_EQ(data.getPreviousField(0)(e, d), static_cast<float>(e * nDof_ + d));
      EXPECT_FLOAT_EQ(data.getCurrentField(0)(e, d), static_cast<float>(e * nDof_ + d) * 2.0f);
    }

  for (int s = 0; s < nSrc_; ++s) {
    EXPECT_EQ(data.getRhsElement()(s), s * 3);
    for (int t = 0; t < nSample_; ++t) EXPECT_FLOAT_EQ(data.getRhsTerm(0)(s, t), static_cast<float>(s * 10 + t));
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(data.getRhsWeights()(s, d), static_cast<float>(d) * 0.1f);
  }
}

TEST_F(DGsolverDataAcousticTest, IsDistributed_DefaultsFalse) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  EXPECT_FALSE(data.isDistributed);
}

// ============================================================
// Accessors
// ============================================================

TEST_F(DGsolverDataAcousticTest, GetCurrentField_ReturnsPnCurr) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  auto curr = data.getCurrentField(0);

  ASSERT_EQ(curr.extent(0), static_cast<size_t>(nElem_));
  ASSERT_EQ(curr.extent(1), static_cast<size_t>(nDof_));
  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(curr(e, d), static_cast<float>(e * nDof_ + d) * 2.0f);
}

TEST_F(DGsolverDataAcousticTest, GetPreviousField_ReturnsPnPrev) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  auto prev = data.getPreviousField(0);

  ASSERT_EQ(prev.extent(0), static_cast<size_t>(nElem_));
  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(prev(e, d), static_cast<float>(e * nDof_ + d));
}

TEST_F(DGsolverDataAcousticTest, GetRhsTerm_ReturnsTerm) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  auto term = data.getRhsTerm(0);

  ASSERT_EQ(term.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(term.extent(1), static_cast<size_t>(nSample_));
  for (int s = 0; s < nSrc_; ++s)
    for (int t = 0; t < nSample_; ++t) EXPECT_FLOAT_EQ(term(s, t), static_cast<float>(s * 10 + t));
}

TEST_F(DGsolverDataAcousticTest, GetRhsElement_ReturnsElem) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  auto elem = data.getRhsElement();

  ASSERT_EQ(elem.extent(0), static_cast<size_t>(nSrc_));
  for (int s = 0; s < nSrc_; ++s) EXPECT_EQ(elem(s), s * 3);
}

TEST_F(DGsolverDataAcousticTest, GetRhsWeights_ReturnsWeights) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  auto w = data.getRhsWeights();

  ASSERT_EQ(w.extent(0), static_cast<size_t>(nSrc_));
  ASSERT_EQ(w.extent(1), static_cast<size_t>(nDof_));
  for (int s = 0; s < nSrc_; ++s)
    for (int d = 0; d < nDof_; ++d) EXPECT_FLOAT_EQ(w(s, d), static_cast<float>(d) * 0.1f);
}

// ============================================================
// swapWavefields
// ============================================================

TEST_F(DGsolverDataAcousticTest, SwapWavefields_ExchangesPrevAndCurr) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  float const prev00 = data.getPreviousField(0)(0, 0);
  float const curr00 = data.getCurrentField(0)(0, 0);

  data.swapWavefields();

  EXPECT_FLOAT_EQ(data.getPreviousField(0)(0, 0), curr00);
  EXPECT_FLOAT_EQ(data.getCurrentField(0)(0, 0), prev00);
}

TEST_F(DGsolverDataAcousticTest, SwapWavefields_TwiceRestoresOriginal) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  float const prev00 = data.getPreviousField(0)(0, 0);
  float const curr00 = data.getCurrentField(0)(0, 0);

  data.swapWavefields();
  data.swapWavefields();

  EXPECT_FLOAT_EQ(data.getPreviousField(0)(0, 0), prev00);
  EXPECT_FLOAT_EQ(data.getCurrentField(0)(0, 0), curr00);
}

TEST_F(DGsolverDataAcousticTest, SwapWavefields_IsShallowHandleSwap) {
  // After swap: pnPrev points to pCurr_ allocation and vice-versa.
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  data.swapWavefields();

  // pnPrev now references the old pCurr_ allocation
  EXPECT_FLOAT_EQ(data.getPreviousField(0)(0, 0), pCurr_(0, 0));
  // pnCurr now references the old pPrev_ allocation
  EXPECT_FLOAT_EQ(data.getCurrentField(0)(0, 0), pPrev_(0, 0));
}

TEST_F(DGsolverDataAcousticTest, SwapWavefields_AllElementsSwapped) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  data.swapWavefields();

  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < nDof_; ++d) {
      EXPECT_FLOAT_EQ(data.getPreviousField(0)(e, d), pCurr_(e, d));
      EXPECT_FLOAT_EQ(data.getCurrentField(0)(e, d), pPrev_(e, d));
    }
}

// ============================================================
// Misc
// ============================================================

TEST_F(DGsolverDataAcousticTest, Print_DoesNotCrash) {
  DGWavefieldAcoustic wavefield(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);

  DGsolverDataAcoustic data(wavefield, rhs);

  EXPECT_NO_THROW(data.print());
}

TEST_F(DGsolverDataAcousticTest, EmptyViews_ValidConstruction) {
  auto emptyPrev = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "ep");
  auto emptyCurr = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "ec");
  auto emptyTerm = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "et");
  auto emptyElem = allocateVector<VECTOR_INT_VIEW>(0, "ee");
  auto emptyW = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "ew");

  DGWavefieldAcoustic empty_wavefield(emptyPrev, emptyCurr);
  RhsAcoustic empty_rhs(emptyTerm, emptyElem, emptyW);

  DGsolverDataAcoustic data(empty_wavefield, empty_rhs);

  EXPECT_EQ(data.getPreviousField(0).extent(0), 0u);
  EXPECT_EQ(data.getCurrentField(0).extent(0), 0u);
  EXPECT_EQ(data.getRhsTerm(0).extent(0), 0u);
  EXPECT_EQ(data.getRhsElement().extent(0), 0u);
  EXPECT_FALSE(data.isDistributed);
  EXPECT_NO_THROW(data.swapWavefields());
}

}  // namespace test
}  // namespace fe
}  // namespace solver

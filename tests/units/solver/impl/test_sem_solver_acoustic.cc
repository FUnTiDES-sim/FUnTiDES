/**
 * @file test_sem_solver_acoustic.cc
 * @brief Unit tests for single-physics acoustic SEMsolver.
 *
 * Covers computeOneStep, computeForces, updateSolution, computeFEInit,
 * mass matrix, damping matrix, and resetGlobalVectors — paths not exercised
 * by the attenuation tests (no SLS) and not reachable via the acoustoelastic
 * composite solver.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>

#include "cartesian_struct_builder.h"
#include "common_macros.h"
#include "data_type.h"
#include "rhs_acoustic.h"
#include "sem_solver.h"
#include "sem_solver_data.h"
#include "solver_factory.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {
namespace test {

namespace feenum = utils::enums;

// ======================================================================
// Fixture
// ======================================================================
struct AcousticSolverOrderParam {
  int order;
};

class SemSolverAcousticTest : public ::testing::TestWithParam<AcousticSolverOrderParam> {
 protected:
  void SetUp() override {
    int order = GetParam().order;
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;

    switch (order) {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, false);
        mesh_ = b.getModel(false);
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> b(EX, LX, EY, LY, EZ, LZ, false, false);
        mesh_ = b.getModel(false);
        break;
      }
      default: {
        model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, false);
        mesh_ = b.getModel(false);
        break;
      }
    }

    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, order);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    int npp = (order + 1) * (order + 1) * (order + 1);

    pPrev_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "pPrev");
    pCurr_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "pCurr");
    for (int i = 0; i < numNodes_; ++i) {
      pPrev_(i) = 0.0f;
      pCurr_(i) = 0.0f;
    }
    pCurr_(numNodes_ / 2) = 1.0f;

    rhsTerm_ = allocateArray2D<ARRAY_REAL_VIEW>(1, kNumSteps, "rhsTerm");
    rhsElem_ = allocateVector<VECTOR_INT_VIEW>(1, "rhsElem");
    rhsWeights_ = allocateArray2D<ARRAY_REAL_VIEW>(1, npp, "rhsWeights");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) {
      rhsTerm_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  static constexpr int kNumSteps = 50;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  VECTOR_REAL_VIEW pPrev_;
  VECTOR_REAL_VIEW pCurr_;
  ARRAY_REAL_VIEW rhsTerm_;
  VECTOR_INT_VIEW rhsElem_;
  ARRAY_REAL_VIEW rhsWeights_;
};

INSTANTIATE_TEST_SUITE_P(AcousticOrders, SemSolverAcousticTest,
                         ::testing::Values(AcousticSolverOrderParam{1}, AcousticSolverOrderParam{2}));

// ======================================================================
// computeFEInit
// ======================================================================
TEST_P(SemSolverAcousticTest, ComputeFEInitDoesNotCrash) { EXPECT_NE(solver_, nullptr); }

// ======================================================================
// Mass matrix
// ======================================================================
TEST_P(SemSolverAcousticTest, MassMatrixNonZero) {
  auto& mm = solver_->getMassMatrixAcoustic();
  ASSERT_GT(mm.extent(0), 0u);
  float sum = 0.0f;
  for (size_t i = 0; i < mm.extent(0); ++i) sum += mm(i);
  EXPECT_GT(sum, 0.0f);
}

TEST_P(SemSolverAcousticTest, MassMatrixAllPositive) {
  auto& mm = solver_->getMassMatrixAcoustic();
  for (size_t i = 0; i < mm.extent(0); ++i) EXPECT_GT(mm(i), 0.0f) << "mass matrix zero at node " << i;
}

// ======================================================================
// Damping matrix
// ======================================================================
TEST_P(SemSolverAcousticTest, DampingMatrixNonNegative) {
  auto& dm = solver_->getDampingMatrix(0);
  for (size_t i = 0; i < dm.extent(0); ++i) EXPECT_GE(dm(i), 0.0f) << "negative damping at node " << i;
}

TEST_P(SemSolverAcousticTest, DampingMatrixHasBoundaryContribution) {
  auto& dm = solver_->getDampingMatrix(0);
  float sum = 0.0f;
  for (size_t i = 0; i < dm.extent(0); ++i) sum += dm(i);
  EXPECT_GT(sum, 0.0f);
}

// ======================================================================
// resetGlobalVectors
// ======================================================================
TEST_P(SemSolverAcousticTest, ResetGlobalVectorsZerosForceVector) {
  auto& fv = solver_->getForceVector(0);
  fv(0) = 99.0f;
  solver_->resetGlobalVectors(numNodes_);
  Kokkos::fence();
  EXPECT_FLOAT_EQ(fv(0), 0.0f);
}

// ======================================================================
// computeForces / updateSolution
// ======================================================================
TEST_P(SemSolverAcousticTest, ComputeForcesDoesNotCrash) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeForces(kDt, 0, data));
}

TEST_P(SemSolverAcousticTest, UpdateSolutionDoesNotCrash) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_NO_THROW(solver_->updateSolution(kDt, data));
}

// ======================================================================
// computeOneStep
// ======================================================================
TEST_P(SemSolverAcousticTest, ComputeOneStepDoesNotCrash) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_P(SemSolverAcousticTest, ComputeOneStepProducesFiniteValues) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(data.getCurrentField(0)(i))) << "NaN/Inf at node " << i;
}

TEST_P(SemSolverAcousticTest, ComputeOneStepZeroSourceStaysBounded) {
  for (int i = 0; i < numNodes_; ++i) {
    pPrev_(i) = 0.0f;
    pCurr_(i) = 0.0f;
  }
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  float sum = 0.0f;
  for (int i = 0; i < numNodes_; ++i) sum += std::fabs(data.getCurrentField(0)(i));
  EXPECT_FLOAT_EQ(sum, 0.0f);
}

TEST_P(SemSolverAcousticTest, ComputeOneStepPropagatesEnergy) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  float norm = 0.0f;
  for (int i = 0; i < numNodes_; ++i) norm += data.getCurrentField(0)(i) * data.getCurrentField(0)(i);
  EXPECT_GT(norm, 0.0f);
}

// ======================================================================
// Accessors
// ======================================================================
TEST_P(SemSolverAcousticTest, GetNumComponentsReturnsOne) { EXPECT_EQ(solver_->getNumComponents(), 1); }

TEST_P(SemSolverAcousticTest, MassMatrixSizeMatchesNodeCount) {
  EXPECT_EQ(static_cast<int>(solver_->getMassMatrixAcoustic().extent(0)), numNodes_);
}

// ======================================================================
// initSpongeValues / initFEarrays
// ======================================================================
TEST_P(SemSolverAcousticTest, InitSpongeValuesDoesNotCrash) { EXPECT_NO_THROW(solver_->initSpongeValues()); }

TEST_P(SemSolverAcousticTest, InitFEarraysDoesNotCrash) { EXPECT_NO_THROW(solver_->initFEarrays()); }

// ======================================================================
// outputSolutionValues
// ======================================================================
TEST_P(SemSolverAcousticTest, OutputSolutionValuesDoesNotCrash) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  int e = 0;
  EXPECT_NO_THROW(solver_->outputSolutionValues(0, e, pCurr_, "pressure"));
}

}  // namespace test
}  // namespace fe
}  // namespace solver

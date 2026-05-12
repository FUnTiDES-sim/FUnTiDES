/**
 * @file test_sem_solver_elastic.cc
 * @brief Unit tests for single-physics elastic SEMsolver.
 *
 * Covers computeOneStep, computeForces, updateSolution, computeFEInit,
 * mass matrix, damping matrix, and resetGlobalVectors on the elastic
 * physics path of SEMsolver without SLS attenuation.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>

#include "cartesian_struct_builder.h"
#include "common_macros.h"
#include "data_type.h"
#include "rhs_elastic.h"
#include "sem_solver.h"
#include "sem_solver_data.h"
#include "solver_factory.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {
namespace test {

namespace feenum = utils::enums;

// ======================================================================
// Fixture
// ======================================================================
struct ElasticSolverOrderParam {
  int order;
};

class SemSolverElasticTest : public ::testing::TestWithParam<ElasticSolverOrderParam> {
 protected:
  void SetUp() override {
    int order = GetParam().order;
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;

    switch (order) {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, true);
        mesh_ = b.getModel(false);
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> b(EX, LX, EY, LY, EZ, LZ, false, true);
        mesh_ = b.getModel(false);
        break;
      }
      default: {
        model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, true);
        mesh_ = b.getModel(false);
        break;
      }
    }

    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, order);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    int npp = (order + 1) * (order + 1) * (order + 1);

    uxPrev_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "uxPrev");
    uxCurr_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "uxCurr");
    uyPrev_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "uyPrev");
    uyCurr_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "uyCurr");
    uzPrev_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "uzPrev");
    uzCurr_ = allocateVector<VECTOR_REAL_VIEW>(numNodes_, "uzCurr");
    for (int i = 0; i < numNodes_; ++i) {
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = 0.0f;
      uzPrev_(i) = uzCurr_(i) = 0.0f;
    }
    uzCurr_(numNodes_ / 2) = 1.0f;

    rhsTermx_ = allocateArray2D<ARRAY_REAL_VIEW>(1, kNumSteps, "rhsTermx");
    rhsTermy_ = allocateArray2D<ARRAY_REAL_VIEW>(1, kNumSteps, "rhsTermy");
    rhsTermz_ = allocateArray2D<ARRAY_REAL_VIEW>(1, kNumSteps, "rhsTermz");
    rhsElem_ = allocateVector<VECTOR_INT_VIEW>(1, "rhsElem");
    rhsWeights_ = allocateArray2D<ARRAY_REAL_VIEW>(1, npp, "rhsWeights");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) {
      rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  static constexpr int kNumSteps = 50;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  VECTOR_REAL_VIEW uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  ARRAY_REAL_VIEW rhsTermx_, rhsTermy_, rhsTermz_;
  VECTOR_INT_VIEW rhsElem_;
  ARRAY_REAL_VIEW rhsWeights_;
};

INSTANTIATE_TEST_SUITE_P(ElasticOrders, SemSolverElasticTest,
                         ::testing::Values(ElasticSolverOrderParam{1}, ElasticSolverOrderParam{2}));

// ======================================================================
// computeFEInit
// ======================================================================
TEST_P(SemSolverElasticTest, ComputeFEInitDoesNotCrash) { EXPECT_NE(solver_, nullptr); }

// ======================================================================
// Mass matrix
// ======================================================================
TEST_P(SemSolverElasticTest, MassMatrixNonZero) {
  auto& mm = solver_->getMassMatrixElastic();
  ASSERT_GT(mm.extent(0), 0u);
  float sum = 0.0f;
  for (size_t i = 0; i < mm.extent(0); ++i) sum += mm(i);
  EXPECT_GT(sum, 0.0f);
}

TEST_P(SemSolverElasticTest, MassMatrixAllPositive) {
  auto& mm = solver_->getMassMatrixElastic();
  for (size_t i = 0; i < mm.extent(0); ++i) EXPECT_GT(mm(i), 0.0f) << "mass matrix zero at node " << i;
}

// ======================================================================
// Damping matrix (3 components for elastic)
// ======================================================================
TEST_P(SemSolverElasticTest, DampingMatrixNonNegativeAllComponents) {
  for (int c = 0; c < 3; ++c) {
    auto& dm = solver_->getDampingMatrix(c);
    for (size_t i = 0; i < dm.extent(0); ++i)
      EXPECT_GE(dm(i), 0.0f) << "negative damping component " << c << " at node " << i;
  }
}

TEST_P(SemSolverElasticTest, DampingMatrixHasBoundaryContribution) {
  float sum = 0.0f;
  for (int c = 0; c < 3; ++c) {
    auto& dm = solver_->getDampingMatrix(c);
    for (size_t i = 0; i < dm.extent(0); ++i) sum += dm(i);
  }
  EXPECT_GT(sum, 0.0f);
}

// ======================================================================
// resetGlobalVectors
// ======================================================================
TEST_P(SemSolverElasticTest, ResetGlobalVectorsZerosForceVector) {
  for (int c = 0; c < 3; ++c) {
    auto& fv = solver_->getForceVector(c);
    fv(0) = 99.0f;
  }
  solver_->resetGlobalVectors(numNodes_);
  Kokkos::fence();
  for (int c = 0; c < 3; ++c) EXPECT_FLOAT_EQ(solver_->getForceVector(c)(0), 0.0f) << "component " << c;
}

// ======================================================================
// computeForces / updateSolution
// ======================================================================
TEST_P(SemSolverElasticTest, ComputeForcesDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeForces(kDt, 0, data));
}

TEST_P(SemSolverElasticTest, UpdateSolutionDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_NO_THROW(solver_->updateSolution(kDt, data));
}

// ======================================================================
// computeOneStep
// ======================================================================
TEST_P(SemSolverElasticTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_P(SemSolverElasticTest, ComputeOneStepProducesFiniteValues) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes_; ++i)
      EXPECT_TRUE(std::isfinite(data.getCurrentField(c)(i))) << "NaN/Inf field " << c << " node " << i;
}

TEST_P(SemSolverElasticTest, ComputeOneStepZeroSourceStaysBounded) {
  for (int i = 0; i < numNodes_; ++i)
    uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;

  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int c = 0; c < 3; ++c) {
    float sum = 0.0f;
    for (int i = 0; i < numNodes_; ++i) sum += std::fabs(data.getCurrentField(c)(i));
    EXPECT_FLOAT_EQ(sum, 0.0f) << "non-zero field " << c << " with zero source";
  }
}

TEST_P(SemSolverElasticTest, ComputeOneStepPropagatesEnergy) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  float norm = 0.0f;
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes_; ++i) norm += data.getCurrentField(c)(i) * data.getCurrentField(c)(i);
  EXPECT_GT(norm, 0.0f);
}

// ======================================================================
// Accessors
// ======================================================================
TEST_P(SemSolverElasticTest, GetNumComponentsReturnsThree) { EXPECT_EQ(solver_->getNumComponents(), 3); }

TEST_P(SemSolverElasticTest, MassMatrixSizeMatchesNodeCount) {
  EXPECT_EQ(static_cast<int>(solver_->getMassMatrixElastic().extent(0)), numNodes_);
}

// ======================================================================
// initSpongeValues / initFEarrays
// ======================================================================
TEST_P(SemSolverElasticTest, InitSpongeValuesDoesNotCrash) { EXPECT_NO_THROW(solver_->initSpongeValues()); }

TEST_P(SemSolverElasticTest, InitFEarraysDoesNotCrash) { EXPECT_NO_THROW(solver_->initFEarrays()); }

// ======================================================================
// outputSolutionValues
// ======================================================================
TEST_P(SemSolverElasticTest, OutputSolutionValuesDoesNotCrash) {
  int e = 0;
  EXPECT_NO_THROW(solver_->outputSolutionValues(0, e, uzCurr_, "uz"));
}

}  // namespace test
}  // namespace fe
}  // namespace solver

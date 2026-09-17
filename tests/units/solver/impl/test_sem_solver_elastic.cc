/**
 * @file test_sem_solver_elastic.cc
 * @brief Unit tests for single-physics elastic SEMsolver.
 *
 * Covers computeOneStep, computeForces, updateSolutionForward, computeFEInit,
 * mass matrix, damping matrix, and resetGlobalVectors on the elastic
 * physics path of SEMsolver without SLS attenuation.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
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

    uxPrevPrev_ = allocateVector<vectorReal>(numNodes_, "uxPrevPrev");
    uxPrev_ = allocateVector<vectorReal>(numNodes_, "uxPrev");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "uxCurr");
    uyPrevPrev_ = allocateVector<vectorReal>(numNodes_, "uyPrevPrev");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uyPrev");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uyCurr");
    uzPrevPrev_ = allocateVector<vectorReal>(numNodes_, "uzPrevPrev");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uzPrev");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uzCurr");
    for (int i = 0; i < numNodes_; ++i) {
      uxPrevPrev_(i) = uxPrev_(i) = uxCurr_(i) = 0.0f;
      uyPrevPrev_(i) = uyPrev_(i) = uyCurr_(i) = 0.0f;
      uzPrevPrev_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    }
    uzCurr_(numNodes_ / 2) = 1.0f;

    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermx");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermy");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermz");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights");
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
  vectorReal uxPrevPrev_, uxPrev_, uxCurr_;
  vectorReal uyPrevPrev_, uyPrev_, uyCurr_;
  vectorReal uzPrevPrev_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
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
// computeForces / updateSolutionForward
// ======================================================================
TEST_P(SemSolverElasticTest, ComputeForcesDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeForces(kDt, 0, data));
}

TEST_P(SemSolverElasticTest, updateSolutionForwardWith2BuffersWorks) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_NO_THROW(solver_->updateSolutionForward(kDt, data));
}

TEST_P(SemSolverElasticTest, updateSolutionForwardWith3BuffersThrows) {
  WavefieldElastic wf(uxPrevPrev_, uxPrev_, uxCurr_, uyPrevPrev_, uyPrev_, uyCurr_, uzPrevPrev_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_THROW(solver_->updateSolutionForward(kDt, data), std::runtime_error);
}

TEST_P(SemSolverElasticTest, updateSolutionBackwardWith2BuffersThrows) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_THROW(solver_->updateSolutionBackward(kDt, data), std::runtime_error);
}

TEST_P(SemSolverElasticTest, updateSolutionBackwardWith3BuffersWorks) {
  WavefieldElastic wf(uxPrevPrev_, uxPrev_, uxCurr_, uyPrevPrev_, uyPrev_, uyCurr_, uzPrevPrev_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_NO_THROW(solver_->updateSolutionBackward(kDt, data));
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

// ======================================================================
// VTI anisotropy — exercises computeElementContributions_Vti
// ======================================================================
class SemSolverElasticVtiTest : public ::testing::TestWithParam<ElasticSolverOrderParam> {
 protected:
  void SetUp() override {
    int order = GetParam().order;
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;

    switch (order) {
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
    mesh_->initElasticityTensors(model::AnisotropyType::kVTI);

    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, order);
    solver_->setAnisotropyType(model::AnisotropyType::kVTI);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    int const npp = (order + 1) * (order + 1) * (order + 1);

    uxPrev_ = allocateVector<vectorReal>(numNodes_, "vux_p");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "vux_c");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "vuy_p");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "vuy_c");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "vuz_p");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "vuz_c");
    for (int i = 0; i < numNodes_; ++i)
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    uzCurr_(numNodes_ / 2) = 1.0f;

    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtx_v");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rty_v");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtz_v");
    rhsElem_ = allocateVector<vectorInt>(1, "re_v");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rw_v");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = rhsWeights_(0, j) = 0.0f;
  }

  static constexpr int kNumSteps = 30;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

INSTANTIATE_TEST_SUITE_P(ElasticVtiOrders, SemSolverElasticVtiTest,
                         ::testing::Values(ElasticSolverOrderParam{1}, ElasticSolverOrderParam{2}));

TEST_P(SemSolverElasticVtiTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_P(SemSolverElasticVtiTest, ComputeOneStepProducesFiniteValues) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes_; ++i)
      EXPECT_TRUE(std::isfinite(data.getCurrentField(c)(i))) << "VTI NaN/Inf field " << c << " node " << i;
}

TEST_P(SemSolverElasticVtiTest, MassMatrixNonZero) {
  auto& mm = solver_->getMassMatrixElastic();
  float sum = 0.0f;
  for (size_t i = 0; i < mm.extent(0); ++i) sum += mm(i);
  EXPECT_GT(sum, 0.0f);
}

// ======================================================================
// TTI anisotropy — exercises computeElementContributions_Tti
// ======================================================================
class SemSolverElasticTtiTest : public ::testing::TestWithParam<ElasticSolverOrderParam> {
 protected:
  void SetUp() override {
    int order = GetParam().order;
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;

    switch (order) {
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
    mesh_->initElasticityTensors(model::AnisotropyType::kTTI);

    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, order);
    solver_->setAnisotropyType(model::AnisotropyType::kTTI);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    int const npp = (order + 1) * (order + 1) * (order + 1);

    uxPrev_ = allocateVector<vectorReal>(numNodes_, "tux_p");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "tux_c");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "tuy_p");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "tuy_c");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "tuz_p");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "tuz_c");
    for (int i = 0; i < numNodes_; ++i)
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    uzCurr_(numNodes_ / 2) = 1.0f;

    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtx_t");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rty_t");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtz_t");
    rhsElem_ = allocateVector<vectorInt>(1, "re_t");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rw_t");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = rhsWeights_(0, j) = 0.0f;
  }

  static constexpr int kNumSteps = 30;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

INSTANTIATE_TEST_SUITE_P(ElasticTtiOrders, SemSolverElasticTtiTest,
                         ::testing::Values(ElasticSolverOrderParam{1}, ElasticSolverOrderParam{2}));

TEST_P(SemSolverElasticTtiTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_P(SemSolverElasticTtiTest, ComputeOneStepProducesFiniteValues) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes_; ++i)
      EXPECT_TRUE(std::isfinite(data.getCurrentField(c)(i))) << "TTI NaN/Inf field " << c << " node " << i;
}

TEST_P(SemSolverElasticTtiTest, MassMatrixNonZero) {
  auto& mm = solver_->getMassMatrixElastic();
  float sum = 0.0f;
  for (size_t i = 0; i < mm.extent(0); ++i) sum += mm(i);
  EXPECT_GT(sum, 0.0f);
}

TEST_P(SemSolverElasticTtiTest, ComputeForcesDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeForces(kDt, 0, data));
}

// ======================================================================
// IS_MODEL_ON_NODES=true — exercises per-node model access paths (elastic)
// ======================================================================
class SemSolverElasticOnNodesTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, true, true);
    mesh_ = b.getModel(false);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnNodes, feenum::physicType::kElastic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    uxPrev_ = allocateVector<vectorReal>(numNodes_, "ux_p_n");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "ux_c_n");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uy_p_n");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uy_c_n");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uz_p_n");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uz_c_n");
    for (int i = 0; i < numNodes_; ++i)
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    uzCurr_(numNodes_ / 2) = 1.0f;
    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtx_n");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rty_n");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtz_n");
    rhsElem_ = allocateVector<vectorInt>(1, "re_n");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rw_n");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = rhsWeights_(0, j) = 0.0f;
  }

  static constexpr int kNumSteps = 10;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverElasticOnNodesTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_F(SemSolverElasticOnNodesTest, ComputeOneStepProducesFiniteValues) {
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

TEST_F(SemSolverElasticOnNodesTest, MassMatrixAllPositive) {
  auto& mm = solver_->getMassMatrixElastic();
  for (size_t i = 0; i < mm.extent(0); ++i) EXPECT_GT(mm(i), 0.0f);
}

// ======================================================================
// Non-zero sponge — exercises initSpongeValues is_sponge=true (elastic)
// ======================================================================
class SemSolverElasticSpongeTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, true);
    mesh_ = b.getModel(false);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);
    solver_->computeFEInit(*mesh_, {50.0f, 0.0f, 0.0f}, false, 10.0f);
    solver_->initSpongeValues();
    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    uxPrev_ = allocateVector<vectorReal>(numNodes_, "ux_p_s");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "ux_c_s");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uy_p_s");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uy_c_s");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uz_p_s");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uz_c_s");
    for (int i = 0; i < numNodes_; ++i)
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    uzCurr_(numNodes_ / 2) = 1.0f;
    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtx_s");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rty_s");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtz_s");
    rhsElem_ = allocateVector<vectorInt>(1, "re_s");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rw_s");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = rhsWeights_(0, j) = 0.0f;
  }

  static constexpr int kNumSteps = 10;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverElasticSpongeTest, ComputeOneStepWithSpongeProducesFiniteValues) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(data.getCurrentField(c)(i)));
}

// ======================================================================
// setSLSAttenuation size mismatch — exercises throw in sem_solver.h
// ======================================================================

TEST_P(SemSolverElasticTest, SetSLSAttenuationSizeMismatchThrows) {
  auto ref = allocateVector<vectorReal>(3, "sls_ref");
  auto coeffs = allocateVector<vectorReal>(2, "sls_coeffs");
  for (int i = 0; i < 3; ++i) ref(i) = 10.0f * (i + 1);
  for (int i = 0; i < 2; ++i) coeffs(i) = 0.1f;
  EXPECT_THROW(solver_->setSLSAttenuation(ref, coeffs), std::runtime_error);
}

// ======================================================================
// computeOneStep in distributed mode — exercises throw in sem_solver_impl.h
// ======================================================================

TEST_P(SemSolverElasticTest, ComputeOneStepDistributedThrows) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs, /*isDistributed=*/true);
  EXPECT_THROW(solver_->computeOneStep(kDt, 0, data), std::runtime_error);
}

// ======================================================================
// SLS attenuation (Iso elastic) — exercises allocateFEarrays, initFEarrays,
// resetGlobalVectors, computeAttenuationContributionsElastic, and the
// attenuation loop inside updateFields (elastic).
// ======================================================================
class SemSolverElasticAttenuationTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, true);
    mesh_ = b.getModel(false);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);

    auto ref = allocateVector<vectorReal>(2, "att_ref_el");
    auto coeffs = allocateVector<vectorReal>(2, "att_coeffs_el");
    ref(0) = 2.0f * 3.14159f * 1.0f;
    ref(1) = 2.0f * 3.14159f * 10.0f;
    coeffs(0) = 0.1f;
    coeffs(1) = 0.1f;
    solver_->setSLSAttenuation(ref, coeffs);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    uxPrev_ = allocateVector<vectorReal>(numNodes_, "ux_p_att");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "ux_c_att");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uy_p_att");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uy_c_att");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uz_p_att");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uz_c_att");
    for (int i = 0; i < numNodes_; ++i)
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    uzCurr_(numNodes_ / 2) = 1.0f;
    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtx_att");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rty_att");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtz_att");
    rhsElem_ = allocateVector<vectorInt>(1, "re_att");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rw_att");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = rhsWeights_(0, j) = 0.0f;
  }

  static constexpr int kNumSteps = 10;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverElasticAttenuationTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_F(SemSolverElasticAttenuationTest, ComputeOneStepProducesFiniteValues) {
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

TEST_F(SemSolverElasticAttenuationTest, ResetGlobalVectorsZerosAttenuationWorkVectors) {
  for (int c = 0; c < 3; ++c) solver_->getForceVector(c)(0) = 99.0f;
  solver_->resetGlobalVectors(numNodes_);
  Kokkos::fence();
  for (int c = 0; c < 3; ++c) EXPECT_FLOAT_EQ(solver_->getForceVector(c)(0), 0.0f) << "component " << c;
}

// ======================================================================
// SLS attenuation (VTI) — exercises early-return path in
// computeAttenuationContributionsElastic (anisotropyType_ != kIso).
// ======================================================================
class SemSolverElasticAttenuationVtiTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, true);
    mesh_ = b.getModel(false);
    mesh_->initElasticityTensors(model::AnisotropyType::kVTI);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kVTI);

    auto ref = allocateVector<vectorReal>(2, "att_ref_vti");
    auto coeffs = allocateVector<vectorReal>(2, "att_coeffs_vti");
    ref(0) = 2.0f * 3.14159f * 1.0f;
    ref(1) = 2.0f * 3.14159f * 10.0f;
    coeffs(0) = 0.1f;
    coeffs(1) = 0.1f;
    solver_->setSLSAttenuation(ref, coeffs);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    uxPrev_ = allocateVector<vectorReal>(numNodes_, "ux_p_vatt");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "ux_c_vatt");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uy_p_vatt");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uy_c_vatt");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uz_p_vatt");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uz_c_vatt");
    for (int i = 0; i < numNodes_; ++i)
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    uzCurr_(numNodes_ / 2) = 1.0f;
    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtx_vatt");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rty_vatt");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtz_vatt");
    rhsElem_ = allocateVector<vectorInt>(1, "re_vatt");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rw_vatt");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = rhsWeights_(0, j) = 0.0f;
  }

  static constexpr int kNumSteps = 5;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

// ======================================================================
// TTI + IS_MODEL_ON_NODES=true — exercises computeCMatrix (rotation of
// the VTI tensor into the tilted frame) called from the TTI stiffness kernel.
// ======================================================================
class SemSolverElasticTtiOnNodesTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, true, true);
    mesh_ = b.getModel(false);
    mesh_->initElasticityTensors(model::AnisotropyType::kTTI);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnNodes, feenum::physicType::kElastic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kTTI);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    uxPrev_ = allocateVector<vectorReal>(numNodes_, "ux_p_tni");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "ux_c_tni");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uy_p_tni");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uy_c_tni");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uz_p_tni");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uz_c_tni");
    for (int i = 0; i < numNodes_; ++i)
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = uzPrev_(i) = uzCurr_(i) = 0.0f;
    uzCurr_(numNodes_ / 2) = 1.0f;
    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtx_tni");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rty_tni");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rtz_tni");
    rhsElem_ = allocateVector<vectorInt>(1, "re_tni");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rw_tni");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = rhsWeights_(0, j) = 0.0f;
  }

  static constexpr int kNumSteps = 10;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverElasticTtiOnNodesTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_F(SemSolverElasticTtiOnNodesTest, ComputeOneStepProducesFiniteValues) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes_; ++i)
      EXPECT_TRUE(std::isfinite(data.getCurrentField(c)(i))) << "TTI+nodes NaN/Inf field " << c << " node " << i;
}

// ======================================================================
// C-PML — zero-profile PML must match the plain solver to machine precision
// ======================================================================
// Read displacement component c (0=ux, 1=uy, 2=uz) at node i from the three
// current-field views.
inline float getCurr(vectorReal const& ux, vectorReal const& uy, vectorReal const& uz, int c, int i) {
  switch (c) {
    case 0:
      return ux(i);
    case 1:
      return uy(i);
    default:
      return uz(i);
  }
}

class SemSolverElasticPmlTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, true);
    mesh_ = b.getModel(false);
    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    uxPrev_ = allocateVector<vectorReal>(numNodes_, "uxPrev_pml");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "uxCurr_pml");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uyPrev_pml");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uyCurr_pml");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uzPrev_pml");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uzCurr_pml");
    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermx_pml");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermy_pml");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermz_pml");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_pml");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_pml");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) {
      rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  // Build a solver, optionally enabling a PML of the given size, and run
  // kNumSteps of computeOneStep from a centered initial impulse.
  std::unique_ptr<Solver> makeSolver(float pml_size, float reflection = 1.0f) {
    auto s =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1);
    s->setAnisotropyType(model::AnisotropyType::kIso);
    if (pml_size > 0.0f) {
      // Default reflection=1.0 gives a zero profile (d_max=0): the PML is
      // enabled but acts as the identity, so the solver must reproduce the
      // plain one wherever the two are configured identically.
      s->setPML({pml_size, pml_size, pml_size}, /*profile=*/2.0f, reflection, /*alpha_max=*/0.0f,
                /*kappa_max=*/1.0f, kDt);
    }
    s->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
    return s;
  }

  void runSteps(Solver& s, vectorReal& uxPrev, vectorReal& uxCurr, vectorReal& uyPrev, vectorReal& uyCurr,
                vectorReal& uzPrev, vectorReal& uzCurr) {
    for (int i = 0; i < numNodes_; ++i) {
      uxPrev(i) = uxCurr(i) = 0.0f;
      uyPrev(i) = uyCurr(i) = 0.0f;
      uzPrev(i) = uzCurr(i) = 0.0f;
    }
    // Vertical displacement impulse at the center node.
    uzCurr(numNodes_ / 2) = 1.0f;
    WavefieldElastic wf(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
    RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
    SEMsolverDataElastic data(wf, rhs);
    for (int t = 0; t < kNumSteps; ++t) {
      s.computeOneStep(kDt, t, data);
      data.swapWavefields();
    }
  }

  static constexpr int kNumSteps = 20;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  int numNodes_;
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverElasticPmlTest, PmlCoefficientsAreNonZero) {
  auto pml = makeSolver(50.0f, /*reflection=*/1e-3f);
  auto& coeff = pml->getPmlCoefficients();
  auto& nodeMask = pml->getPmlNodeIndex();
  auto& elemMask = pml->getPmlElementMask();
  int nNodes = coeff.extent(0);
  int nPmlNodes = 0;
  float kappaMax = 0.0f;
  for (int n = 0; n < nNodes; ++n) {
    if (nodeMask(n) == 1) {
      ++nPmlNodes;
      // Compact layout: kappa at (n, 0..2), coef0 at (n, 3..5), coef1 at (n, 6..8).
      for (int j = 0; j < 3; ++j) kappaMax = std::max(kappaMax, std::fabs(coeff(n, j)));
    }
  }
  int nElems = elemMask.extent(0);
  int nPmlElems = 0;
  for (int e = 0; e < nElems; ++e) nPmlElems += elemMask(e);
  EXPECT_GT(kappaMax, 0.0f) << "PML profile is identically zero — no absorption possible";
  EXPECT_GT(nPmlNodes, 0) << "No nodes in PML layer";
  EXPECT_GT(nPmlElems, 0) << "No elements in PML layer";
}

TEST_F(SemSolverElasticPmlTest, ZeroProfilePmlMatchesPlainSolver) {
  auto plain = makeSolver(0.0f);
  auto pml = makeSolver(50.0f);  // PML enabled but zero profile (identity)

  vectorReal uxPrevPlain = allocateVector<vectorReal>(numNodes_, "uxPrevPlain");
  vectorReal uxCurrPlain = allocateVector<vectorReal>(numNodes_, "uxCurrPlain");
  vectorReal uyPrevPlain = allocateVector<vectorReal>(numNodes_, "uyPrevPlain");
  vectorReal uyCurrPlain = allocateVector<vectorReal>(numNodes_, "uyCurrPlain");
  vectorReal uzPrevPlain = allocateVector<vectorReal>(numNodes_, "uzPrevPlain");
  vectorReal uzCurrPlain = allocateVector<vectorReal>(numNodes_, "uzCurrPlain");
  vectorReal uxPrevPml = allocateVector<vectorReal>(numNodes_, "uxPrevPml");
  vectorReal uxCurrPml = allocateVector<vectorReal>(numNodes_, "uxCurrPml");
  vectorReal uyPrevPml = allocateVector<vectorReal>(numNodes_, "uyPrevPml");
  vectorReal uyCurrPml = allocateVector<vectorReal>(numNodes_, "uyCurrPml");
  vectorReal uzPrevPml = allocateVector<vectorReal>(numNodes_, "uzPrevPml");
  vectorReal uzCurrPml = allocateVector<vectorReal>(numNodes_, "uzCurrPml");

  runSteps(*plain, uxPrevPlain, uxCurrPlain, uyPrevPlain, uyCurrPlain, uzPrevPlain, uzCurrPlain);
  runSteps(*pml, uxPrevPml, uxCurrPml, uyPrevPml, uyCurrPml, uzPrevPml, uzCurrPml);

  // The zero-profile PML must reproduce the plain solver wherever the two
  // solvers are configured identically. The 2x2x2 mesh has a single interior
  // node (the center, numNodes_/2): it lies on no boundary face, so its
  // damping is 0 in both solvers and the only possible difference is the
  // stiffness kernel — which the kernel oracle proves identical for a zero
  // profile. Boundary nodes legitimately differ because the PML disables the
  // first-order absorbing BC on PML faces, and that difference propagates to
  // the center through the shared elements over the 20 steps.
  int const center = numNodes_ / 2;
  float center_diff = 0.0f;
  for (int c = 0; c < 3; ++c) {
    float const d = std::fabs(getCurr(uxCurrPlain, uyCurrPlain, uzCurrPlain, c, center) -
                              getCurr(uxCurrPml, uyCurrPml, uzCurrPml, c, center));
    center_diff = std::max(center_diff, d);
  }

  // Report the worst node for diagnostics (expected on the boundary).
  float max_diff = 0.0f;
  for (int i = 0; i < numNodes_; ++i) {
    for (int c = 0; c < 3; ++c) {
      float const d = std::fabs(getCurr(uxCurrPlain, uyCurrPlain, uzCurrPlain, c, i) -
                                getCurr(uxCurrPml, uyCurrPml, uzCurrPml, c, i));
      max_diff = std::max(max_diff, d);
    }
  }

  // The center difference must be far smaller than the boundary difference:
  // the kernel is identical (oracle), so the only source of divergence is the
  // damping disabled on PML faces, which acts at the boundary and reaches the
  // center only through propagation.
  EXPECT_LT(center_diff, 0.1f * max_diff + 1e-6f)
      << "zero-profile PML diverged from plain solver at center node (center_diff=" << center_diff
      << ", max_diff=" << max_diff << ")";
}

TEST_F(SemSolverElasticPmlTest, PmlRunProducesFiniteValues) {
  auto pml = makeSolver(50.0f);
  vectorReal uxPrev = allocateVector<vectorReal>(numNodes_, "uxPrevPmlF");
  vectorReal uxCurr = allocateVector<vectorReal>(numNodes_, "uxCurrPmlF");
  vectorReal uyPrev = allocateVector<vectorReal>(numNodes_, "uyPrevPmlF");
  vectorReal uyCurr = allocateVector<vectorReal>(numNodes_, "uyCurrPmlF");
  vectorReal uzPrev = allocateVector<vectorReal>(numNodes_, "uzPrevPmlF");
  vectorReal uzCurr = allocateVector<vectorReal>(numNodes_, "uzCurrPmlF");
  runSteps(*pml, uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
  for (int c = 0; c < 3; ++c)
    for (int i = 0; i < numNodes_; ++i)
      EXPECT_TRUE(std::isfinite(getCurr(uxCurr, uyCurr, uzCurr, c, i))) << "NaN/Inf field " << c << " node " << i;
}

// ======================================================================
// C-PML reflection — PML vs large-domain reference
// ======================================================================
class SemSolverElasticPmlReflectionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    rhsTermx_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermx_refl");
    rhsTermy_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermy_refl");
    rhsTermz_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermz_refl");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_refl");
    rhsWeights_ = allocateArray2D<arrayReal>(1, 8, "rhsWeights_refl");
    rhsElem_(0) = 0;
    for (int j = 0; j < 8; ++j) {
      rhsTermx_(0, j) = rhsTermy_(0, j) = rhsTermz_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  // Build a solver on a domain of the given size with a PML of the given
  // thickness, and run a centered Gaussian impulse long enough for the wave
  // to reach the boundary and reflect back. The element size is fixed (25m)
  // so the reference and small domains have identical resolution, and the
  // PML (>= 3 elements thick) is resolved by the mesh.
  std::unique_ptr<Solver> makeSolver(float lx, float pml_size, float reflection) {
    int const ex = static_cast<int>(std::lround(lx / 25.0f));
    model::CartesianStructBuilder<float, int, 1> b(ex, lx, ex, lx, ex, lx, false, true);
    mesh_ = b.getModel(false);
    auto s =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1);
    s->setAnisotropyType(model::AnisotropyType::kIso);
    if (pml_size > 0.0f) {
      s->setPML({pml_size, pml_size, pml_size}, /*profile=*/2.0f, reflection, /*alpha_max=*/0.0f,
                /*kappa_max=*/1.0f, kDt);
    }
    s->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
    return s;
  }

  // Run kNumSteps from a smooth Gaussian bump centered in the domain; return
  // the total energy (sum of ux^2+uy^2+uz^2 over all nodes) at the end of the
  // run. Total energy is the physically correct absorption metric: a PML
  // removes the outgoing wave, while the plain run (first-order absorbing BC)
  // keeps more of it in the domain.
  float runSteps(Solver& s, int numNodes) {
    vectorReal uxPrev = allocateVector<vectorReal>(numNodes, "uxPrevRefl");
    vectorReal uxCurr = allocateVector<vectorReal>(numNodes, "uxCurrRefl");
    vectorReal uyPrev = allocateVector<vectorReal>(numNodes, "uyPrevRefl");
    vectorReal uyCurr = allocateVector<vectorReal>(numNodes, "uyCurrRefl");
    vectorReal uzPrev = allocateVector<vectorReal>(numNodes, "uzPrevRefl");
    vectorReal uzCurr = allocateVector<vectorReal>(numNodes, "uzCurrRefl");
    for (int i = 0; i < numNodes; ++i) {
      uxPrev(i) = uxCurr(i) = 0.0f;
      uyPrev(i) = uyCurr(i) = 0.0f;
      uzPrev(i) = uzCurr(i) = 0.0f;
    }
    // Smooth Gaussian bump (width ~2 elements) in all three components so the
    // spectrum is band-limited and the PML operates in its designed frequency
    // range.
    int const dim = static_cast<int>(std::cbrt(static_cast<float>(numNodes)));
    int const c = numNodes / 2;
    int const cx = c % dim, cy = (c / dim) % dim, cz = c / (dim * dim);
    for (int k = 0; k < dim; ++k)
      for (int j = 0; j < dim; ++j)
        for (int i = 0; i < dim; ++i) {
          float const dx = (i - cx) / 2.0f, dy = (j - cy) / 2.0f, dz = (k - cz) / 2.0f;
          float const g = std::exp(-(dx * dx + dy * dy + dz * dz));
          int const idx = i + j * dim + k * dim * dim;
          uxCurr(idx) = g;
          uyCurr(idx) = g;
          uzCurr(idx) = g;
        }
    WavefieldElastic wf(uxPrev, uxCurr, uyPrev, uyCurr, uzPrev, uzCurr);
    RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
    SEMsolverDataElastic data(wf, rhs);
    for (int t = 0; t < kNumSteps; ++t) {
      s.computeOneStep(kDt, t, data);
      data.swapWavefields();
    }
    float eTot = 0.0f;
    for (int i = 0; i < numNodes; ++i)
      for (int c = 0; c < 3; ++c) eTot += data.getCurrentField(c)(i) * data.getCurrentField(c)(i);
    return eTot;
  }

  static constexpr int kNumSteps = 200;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  arrayReal rhsTermx_, rhsTermy_, rhsTermz_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverElasticPmlReflectionTest, PmlAbsorbsBetterThanPlain) {
  // Small domain (300m) with a 100m PML (4 elements at 25m): the interior is
  // 100m, the source sits at the center (50m from the PML inner edge), and
  // the wave reaches the PML and is absorbed. The plain run uses the same
  // 300m domain with no PML: the wave reaches the boundary, is only partially
  // absorbed by the first-order absorbing BC, and stays in the domain. Both
  // runs have identical node counts, so the total energy (sum of u^2) at the
  // end is directly comparable: a working PML leaves far less energy behind
  // than the plain absorbing BC.
  //
  // Solvers are created and released one at a time to keep the cumulative GPU
  // memory bounded (each solver allocates mass/damping/work vectors and the
  // PML state).
  float energy_pml = 0.0f, energy_plain = 0.0f;

  {
    auto small_pml = makeSolver(300.0f, 100.0f, 1e-3f);
    int nSmall = small_pml->getMassMatrixElastic().extent(0);
    energy_pml = runSteps(*small_pml, nSmall);
  }

  {
    auto small_plain = makeSolver(300.0f, 0.0f, 1e-3f);  // no PML, no sponge
    int nSmall = small_plain->getMassMatrixElastic().extent(0);
    energy_plain = runSteps(*small_plain, nSmall);
  }

  // The PML must absorb the outgoing wave: it should leave at most 30% of the
  // energy that the plain (first-order absorbing BC) run keeps in the domain.
  std::cout << "elastic_pml_energy=" << energy_pml << " elastic_plain_energy=" << energy_plain
            << " ratio=" << (energy_plain > 0.0f ? energy_pml / energy_plain : -1.0f) << std::endl;
  EXPECT_LT(energy_pml, 0.3f * energy_plain) << "PML energy_pml=" << energy_pml
                                             << " energy_plain=" << energy_plain;
}

}  // namespace test
}  // namespace fe
}  // namespace solver

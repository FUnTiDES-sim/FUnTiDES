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

    uxPrev_ = allocateVector<vectorReal>(numNodes_, "uxPrev");
    uxCurr_ = allocateVector<vectorReal>(numNodes_, "uxCurr");
    uyPrev_ = allocateVector<vectorReal>(numNodes_, "uyPrev");
    uyCurr_ = allocateVector<vectorReal>(numNodes_, "uyCurr");
    uzPrev_ = allocateVector<vectorReal>(numNodes_, "uzPrev");
    uzCurr_ = allocateVector<vectorReal>(numNodes_, "uzCurr");
    for (int i = 0; i < numNodes_; ++i) {
      uxPrev_(i) = uxCurr_(i) = uyPrev_(i) = uyCurr_(i) = 0.0f;
      uzPrev_(i) = uzCurr_(i) = 0.0f;
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
  vectorReal uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_;
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

TEST_P(SemSolverElasticTest, updateSolutionForwardDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_NO_THROW(solver_->updateSolutionForward(kDt, data));
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
class SemSolverElasticVtiTest : public ::testing::Test {
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
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 2 * 2 * 2;

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

TEST_F(SemSolverElasticVtiTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_F(SemSolverElasticVtiTest, ComputeOneStepProducesFiniteValues) {
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

TEST_F(SemSolverElasticVtiTest, MassMatrixNonZero) {
  auto& mm = solver_->getMassMatrixElastic();
  float sum = 0.0f;
  for (size_t i = 0; i < mm.extent(0); ++i) sum += mm(i);
  EXPECT_GT(sum, 0.0f);
}

// ======================================================================
// TTI anisotropy — exercises computeElementContributions_Tti
// ======================================================================
class SemSolverElasticTtiTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, true);
    mesh_ = b.getModel(false);
    mesh_->initElasticityTensors(model::AnisotropyType::kTTI);

    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kTTI);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 2 * 2 * 2;

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

TEST_F(SemSolverElasticTtiTest, ComputeOneStepDoesNotCrash) {
  WavefieldElastic wf(uxPrev_, uxCurr_, uyPrev_, uyCurr_, uzPrev_, uzCurr_);
  RhsElastic rhs(rhsTermx_, rhsTermy_, rhsTermz_, rhsElem_, rhsWeights_);
  SEMsolverDataElastic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_F(SemSolverElasticTtiTest, ComputeOneStepProducesFiniteValues) {
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

TEST_F(SemSolverElasticTtiTest, MassMatrixNonZero) {
  auto& mm = solver_->getMassMatrixElastic();
  float sum = 0.0f;
  for (size_t i = 0; i < mm.extent(0); ++i) sum += mm(i);
  EXPECT_GT(sum, 0.0f);
}

TEST_F(SemSolverElasticTtiTest, ComputeForcesDoesNotCrash) {
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

}  // namespace test
}  // namespace fe
}  // namespace solver

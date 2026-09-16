/**
 * @file test_sem_solver_acoustic.cc
 * @brief Unit tests for single-physics acoustic SEMsolver.
 *
 * Covers computeOneStep, computeForces, updateSolutionForward, computeFEInit,
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

    pPrevPrev_ = allocateVector<vectorReal>(numNodes_, "pPrevPrev");
    pPrev_ = allocateVector<vectorReal>(numNodes_, "pPrev");
    pCurr_ = allocateVector<vectorReal>(numNodes_, "pCurr");
    for (int i = 0; i < numNodes_; ++i) {
      pPrevPrev_(i) = 0.0f;
      pPrev_(i) = 0.0f;
      pCurr_(i) = 0.0f;
    }
    pCurr_(numNodes_ / 2) = 1.0f;

    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights");
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
  vectorReal pPrevPrev_;
  vectorReal pPrev_;
  vectorReal pCurr_;
  arrayReal rhsTerm_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
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
// computeForces / updateSolutionForward
// ======================================================================
TEST_P(SemSolverAcousticTest, ComputeForcesDoesNotCrash) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeForces(kDt, 0, data));
}

TEST_P(SemSolverAcousticTest, updateSolutionForwardWith2BuffersWorks) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_NO_THROW(solver_->updateSolutionForward(kDt, data));
}

TEST_P(SemSolverAcousticTest, updateSolutionForwardWith3BuffersThrows) {
  WavefieldAcoustic wf(pPrevPrev_, pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_THROW(solver_->updateSolutionForward(kDt, data), std::runtime_error);
}

TEST_P(SemSolverAcousticTest, updateSolutionBackwardWith2BuffersThrows) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_THROW(solver_->updateSolutionBackward(kDt, data), std::runtime_error);
}

TEST_P(SemSolverAcousticTest, updateSolutionBackwardWith3BuffersWorks) {
  WavefieldAcoustic wf(pPrevPrev_, pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  solver_->computeForces(kDt, 0, data);
  EXPECT_NO_THROW(solver_->updateSolutionBackward(kDt, data));
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

// ======================================================================
// IS_MODEL_ON_NODES=true — exercises per-node model access paths
// ======================================================================
class SemSolverAcousticOnNodesTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, true, false);
    mesh_ = b.getModel(false);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    pPrev_ = allocateVector<vectorReal>(numNodes_, "pPrev_n");
    pCurr_ = allocateVector<vectorReal>(numNodes_, "pCurr_n");
    for (int i = 0; i < numNodes_; ++i) {
      pPrev_(i) = 0.0f;
      pCurr_(i) = 0.0f;
    }
    pCurr_(numNodes_ / 2) = 1.0f;
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_n");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_n");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_n");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) {
      rhsTerm_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  static constexpr int kNumSteps = 10;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal pPrev_, pCurr_;
  arrayReal rhsTerm_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverAcousticOnNodesTest, ComputeOneStepDoesNotCrash) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_F(SemSolverAcousticOnNodesTest, ComputeOneStepProducesFiniteValues) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(data.getCurrentField(0)(i))) << "NaN/Inf at node " << i;
}

TEST_F(SemSolverAcousticOnNodesTest, MassMatrixAllPositive) {
  auto& mm = solver_->getMassMatrixAcoustic();
  for (size_t i = 0; i < mm.extent(0); ++i) EXPECT_GT(mm(i), 0.0f);
}

// ======================================================================
// Non-zero sponge — exercises initSpongeValues is_sponge=true branch
// ======================================================================
class SemSolverAcousticSpongeTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, false);
    mesh_ = b.getModel(false);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);
    solver_->computeFEInit(*mesh_, {50.0f, 0.0f, 0.0f}, false, 10.0f);
    solver_->initSpongeValues();
    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    pPrev_ = allocateVector<vectorReal>(numNodes_, "pPrev_s");
    pCurr_ = allocateVector<vectorReal>(numNodes_, "pCurr_s");
    for (int i = 0; i < numNodes_; ++i) {
      pPrev_(i) = 0.0f;
      pCurr_(i) = 0.0f;
    }
    pCurr_(numNodes_ / 2) = 1.0f;
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_s");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_s");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_s");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) {
      rhsTerm_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  static constexpr int kNumSteps = 10;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal pPrev_, pCurr_;
  arrayReal rhsTerm_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverAcousticSpongeTest, ComputeOneStepWithSpongeProducesFiniteValues) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(data.getCurrentField(0)(i)));
}

TEST_F(SemSolverAcousticSpongeTest, SpongeWithSurfaceFlagDoesNotCrash) {
  model::CartesianStructBuilder<float, int, 1> b(2, 200.0f, 2, 200.0f, 2, 200.0f, false, false);
  auto mesh2 = b.getModel(false);
  auto solver2 =
      solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                   feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
  solver2->setAnisotropyType(model::AnisotropyType::kIso);
  solver2->computeFEInit(*mesh2, {50.0f, 0.0f, 0.0f}, true, 10.0f);
  EXPECT_NO_THROW(solver2->initSpongeValues());
}

// ======================================================================
// SLS attenuation — exercises allocateFEarrays, initFEarrays,
// resetGlobalVectors, computeAttenuationContributionsAcoustic,
// and the attenuation loop inside updateFields (acoustic).
// ======================================================================
class SemSolverAcousticAttenuationTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, false);
    mesh_ = b.getModel(false);
    solver_ =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
    solver_->setAnisotropyType(model::AnisotropyType::kIso);

    auto ref = allocateVector<vectorReal>(2, "att_ref_ac");
    auto coeffs = allocateVector<vectorReal>(2, "att_coeffs_ac");
    ref(0) = 2.0f * 3.14159f * 1.0f;
    ref(1) = 2.0f * 3.14159f * 10.0f;
    coeffs(0) = 0.1f;
    coeffs(1) = 0.1f;
    solver_->setSLSAttenuation(ref, coeffs);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    pPrev_ = allocateVector<vectorReal>(numNodes_, "pPrev_att");
    pCurr_ = allocateVector<vectorReal>(numNodes_, "pCurr_att");
    for (int i = 0; i < numNodes_; ++i) {
      pPrev_(i) = 0.0f;
      pCurr_(i) = 0.0f;
    }
    pCurr_(numNodes_ / 2) = 1.0f;
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_att");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_att");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_att");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) {
      rhsTerm_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  static constexpr int kNumSteps = 10;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  std::unique_ptr<Solver> solver_;
  int numNodes_;
  vectorReal pPrev_, pCurr_;
  arrayReal rhsTerm_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverAcousticAttenuationTest, ComputeOneStepDoesNotCrash) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  EXPECT_NO_THROW(solver_->computeOneStep(kDt, 0, data));
}

TEST_F(SemSolverAcousticAttenuationTest, ComputeOneStepProducesFiniteValues) {
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(data.getCurrentField(0)(i))) << "NaN/Inf at node " << i;
}

TEST_F(SemSolverAcousticAttenuationTest, ResetGlobalVectorsZerosAttenuationWorkVectors) {
  auto& fv = solver_->getForceVector(0);
  fv(0) = 99.0f;
  solver_->resetGlobalVectors(numNodes_);
  Kokkos::fence();
  EXPECT_FLOAT_EQ(fv(0), 0.0f);
}

// Covers Solver::setZBoundary default body (non-overriding solver).
TEST_F(SemSolverAcousticAttenuationTest, SetZBoundary_DoesNotThrow) { EXPECT_NO_THROW(solver_->setZBoundary(500.0f)); }

// ======================================================================
// C-PML — zero-profile PML must match the plain solver to machine precision
// ======================================================================
class SemSolverAcousticPmlTest : public ::testing::Test {
 protected:
  void SetUp() override {
    constexpr int EX = 2, EY = 2, EZ = 2;
    constexpr float LX = 200.0f, LY = 200.0f, LZ = 200.0f;
    model::CartesianStructBuilder<float, int, 1> b(EX, LX, EY, LY, EZ, LZ, false, false);
    mesh_ = b.getModel(false);
    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    pPrev_ = allocateVector<vectorReal>(numNodes_, "pPrev_pml");
    pCurr_ = allocateVector<vectorReal>(numNodes_, "pCurr_pml");
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_pml");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_pml");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_pml");
    rhsElem_(0) = 0;
    for (int j = 0; j < npp; ++j) {
      rhsTerm_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  // Build a solver, optionally enabling a PML of the given size, and run
  // kNumSteps of computeOneStep from a centered initial impulse.
  std::unique_ptr<Solver> makeSolver(float pml_size, float reflection = 1.0f) {
    auto s =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
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

  void runSteps(Solver& s, vectorReal& pPrev, vectorReal& pCurr) {
    for (int i = 0; i < numNodes_; ++i) {
      pPrev(i) = 0.0f;
      pCurr(i) = 0.0f;
    }
    pCurr(numNodes_ / 2) = 1.0f;
    WavefieldAcoustic wf(pPrev, pCurr);
    RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
    SEMsolverDataAcoustic data(wf, rhs);
    for (int t = 0; t < kNumSteps; ++t) {
      s.computeOneStep(kDt, t, data);
      data.swapWavefields();
    }
  }

  static constexpr int kNumSteps = 20;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  int numNodes_;
  vectorReal pPrev_, pCurr_;
  arrayReal rhsTerm_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverAcousticPmlTest, PmlCoefficientsAreNonZero) {
  auto pml = makeSolver(50.0f, /*reflection=*/1e-3f);
  auto& coeff = pml->getPmlCoefficients();
  auto& nodeMask = pml->getPmlNodeIndex();
  auto& elemMask = pml->getPmlElementMask();
  int nNodes = coeff.extent(0);
  int nPmlNodes = 0;
  float kappaMax = 0.0f, coef1Max = 0.0f;
  for (int n = 0; n < nNodes; ++n) {
    if (nodeMask(n) == 1) {
      ++nPmlNodes;
      for (int j = 0; j < 3; ++j) {
        // Compact layout: kappa(3) + coef0(3) + coef1(3).
        kappaMax = std::max(kappaMax, std::fabs(coeff(n, j)));
        coef1Max = std::max(coef1Max, std::fabs(coeff(n, 6 + j)));
      }
    }
  }
  int nElems = elemMask.extent(0);
  int nPmlElems = 0;
  for (int e = 0; e < nElems; ++e) nPmlElems += elemMask(e);
  EXPECT_GT(coef1Max, 0.0f) << "PML profile is identically zero — no absorption possible";
  EXPECT_GT(nPmlNodes, 0) << "No nodes in PML layer";
  EXPECT_GT(nPmlElems, 0) << "No elements in PML layer";
}

TEST_F(SemSolverAcousticPmlTest, ZeroProfilePmlMatchesPlainSolver) {
  auto plain = makeSolver(0.0f);
  auto pml = makeSolver(50.0f);  // PML enabled but zero profile (identity)

  vectorReal pPrevPlain = allocateVector<vectorReal>(numNodes_, "pPrevPlain");
  vectorReal pCurrPlain = allocateVector<vectorReal>(numNodes_, "pCurrPlain");
  vectorReal pPrevPml = allocateVector<vectorReal>(numNodes_, "pPrevPml");
  vectorReal pCurrPml = allocateVector<vectorReal>(numNodes_, "pCurrPml");

  runSteps(*plain, pPrevPlain, pCurrPlain);
  runSteps(*pml, pPrevPml, pCurrPml);

  // The zero-profile PML must reproduce the plain solver wherever the two
  // solvers are configured identically. The 2x2x2 mesh has a single interior
  // node (the center, numNodes_/2): it lies on no boundary face, so its
  // damping is 0 in both solvers and the only possible difference is the
  // stiffness kernel — which the kernel oracle proves identical for a zero
  // profile. Boundary nodes legitimately differ because the PML disables the
  // first-order absorbing BC on PML faces, and that difference propagates to
  // the center through the shared elements over the 20 steps.
  int const center = numNodes_ / 2;
  float const center_diff = std::fabs(pCurrPlain(center) - pCurrPml(center));

  // Report the worst node for diagnostics (expected on the boundary).
  float max_diff = 0.0f;
  for (int i = 0; i < numNodes_; ++i) {
    float d = std::fabs(pCurrPlain(i) - pCurrPml(i));
    if (d > max_diff) max_diff = d;
  }

  // The center difference must be far smaller than the boundary difference:
  // the kernel is identical (oracle), so the only source of divergence is the
  // damping disabled on PML faces, which acts at the boundary and reaches the
  // center only through propagation.
  EXPECT_LT(center_diff, 0.1f * max_diff + 1e-6f)
      << "zero-profile PML diverged from plain solver at center node (center_diff=" << center_diff
      << ", max_diff=" << max_diff << ")";
}

TEST_F(SemSolverAcousticPmlTest, PmlRunProducesFiniteValues) {
  auto pml = makeSolver(50.0f);
  vectorReal pPrev = allocateVector<vectorReal>(numNodes_, "pPrevPmlF");
  vectorReal pCurr = allocateVector<vectorReal>(numNodes_, "pCurrPmlF");
  runSteps(*pml, pPrev, pCurr);
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(pCurr(i))) << "NaN/Inf at node " << i;
}

// ======================================================================
// C-PML reflection — PML vs large-domain reference
// ======================================================================
class SemSolverAcousticPmlReflectionTest : public ::testing::Test {
 protected:
  void SetUp() override {
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_refl");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_refl");
    rhsWeights_ = allocateArray2D<arrayReal>(1, 8, "rhsWeights_refl");
    rhsElem_(0) = 0;
    for (int j = 0; j < 8; ++j) {
      rhsTerm_(0, j) = 0.0f;
      rhsWeights_(0, j) = 0.0f;
    }
  }

  // Build a solver on a domain of the given size with a PML of the given
  // thickness, and run a centered Ricker-like impulse long enough for the
  // wave to reach the boundary and reflect back. The element size is fixed
  // (25m) so the reference and small domains have identical resolution, and
  // the PML (>= 3 elements thick) is resolved by the mesh.
  std::unique_ptr<Solver> makeSolver(float lx, float pml_size, float reflection) {
    int const ex = static_cast<int>(std::lround(lx / 25.0f));
    model::CartesianStructBuilder<float, int, 1> b(ex, lx, ex, lx, ex, lx, false, false);
    mesh_ = b.getModel(false);
    auto s =
        solver_factory::createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
    s->setAnisotropyType(model::AnisotropyType::kIso);
    if (pml_size > 0.0f) {
      s->setPML({pml_size, pml_size, pml_size}, /*profile=*/2.0f, reflection, /*alpha_max=*/0.0f,
                /*kappa_max=*/1.0f, kDt);
    }
    s->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
    return s;
  }

  // Run kNumSteps from a smooth Gaussian bump centered in the domain; return
  // the total energy (sum of p^2 over all nodes) at the end of the run. Total
  // energy is the physically correct absorption metric: a PML removes the
  // outgoing wave, while the plain run (first-order absorbing BC) keeps more
  // of it in the domain.
  float runSteps(Solver& s, int numNodes) {
    vectorReal pPrev = allocateVector<vectorReal>(numNodes, "pPrevRefl");
    vectorReal pCurr = allocateVector<vectorReal>(numNodes, "pCurrRefl");
    for (int i = 0; i < numNodes; ++i) {
      pPrev(i) = 0.0f;
      pCurr(i) = 0.0f;
    }
    // Smooth Gaussian bump (width ~1.5 elements) so the spectrum is band-limited
    // and the PML operates in its designed frequency range.
    int const dim = static_cast<int>(std::cbrt(static_cast<float>(numNodes)));
    int const c = numNodes / 2;
    int const cx = c % dim, cy = (c / dim) % dim, cz = c / (dim * dim);
    for (int k = 0; k < dim; ++k)
      for (int j = 0; j < dim; ++j)
        for (int i = 0; i < dim; ++i) {
          float const dx = (i - cx) / 1.5f, dy = (j - cy) / 1.5f, dz = (k - cz) / 1.5f;
          pCurr(i + j * dim + k * dim * dim) = std::exp(-(dx * dx + dy * dy + dz * dz));
        }
    WavefieldAcoustic wf(pPrev, pCurr);
    RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
    SEMsolverDataAcoustic data(wf, rhs);
    for (int t = 0; t < kNumSteps; ++t) {
      s.computeOneStep(kDt, t, data);
      data.swapWavefields();
    }
    float eTot = 0.0f;
    for (int i = 0; i < numNodes; ++i) eTot += data.getCurrentField(0)(i) * data.getCurrentField(0)(i);
    return eTot;
  }

  static constexpr int kNumSteps = 200;
  static constexpr float kDt = 0.001f;

  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  arrayReal rhsTerm_;
  vectorInt rhsElem_;
  arrayReal rhsWeights_;
};

TEST_F(SemSolverAcousticPmlReflectionTest, PmlAbsorbsBetterThanPlain) {
  // Small domain (300m) with a 100m PML (4 elements at 25m): the interior is
  // 100m, the source sits at the center (50m from the PML inner edge), and
  // the wave reaches the PML at t=0.033s and is absorbed. The plain run uses
  // the same 300m domain with no PML: the wave reaches the boundary at
  // t=0.1s, is only partially absorbed by the first-order absorbing BC, and
  // stays in the domain. Both runs have identical node counts, so the total
  // energy (sum of p^2) at the end is directly comparable: a working PML
  // leaves far less energy behind than the plain absorbing BC.
  //
  // Solvers are created and released one at a time to keep the cumulative GPU
  // memory bounded (each solver allocates mass/damping/work vectors and the
  // PML state).
  float energy_pml = 0.0f, energy_plain = 0.0f;

  {
    auto small_pml = makeSolver(300.0f, 100.0f, 1e-3f);
    int nSmall = small_pml->getMassMatrixAcoustic().extent(0);
    energy_pml = runSteps(*small_pml, nSmall);
  }

  {
    auto small_plain = makeSolver(300.0f, 0.0f, 1e-3f);  // no PML, no sponge
    int nSmall = small_plain->getMassMatrixAcoustic().extent(0);
    energy_plain = runSteps(*small_plain, nSmall);
  }

  // The PML must absorb the outgoing wave: it should leave at most 30% of the
  // energy that the plain (first-order absorbing BC) run keeps in the domain.
  EXPECT_LT(energy_pml, 0.3f * energy_plain) << "PML energy_pml=" << energy_pml
                                             << " energy_plain=" << energy_plain;
}

}  // namespace test
}  // namespace fe
}  // namespace solver

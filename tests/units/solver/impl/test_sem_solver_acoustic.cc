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
#include "device_view_access.h"
#include "rhs_acoustic.h"
#include "sem_solver.h"
#include "sem_solver_data.h"
#include "solver_factory.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromList;
using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

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
    Kokkos::deep_copy(pPrevPrev_, 0.0f);
    Kokkos::deep_copy(pPrev_, 0.0f);
    Kokkos::deep_copy(pCurr_, 0.0f);
    setAt(pCurr_, 1.0f, numNodes_ / 2);

    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights");
    Kokkos::deep_copy(rhsElem_, 0);
    Kokkos::deep_copy(rhsTerm_, 0.0f);
    Kokkos::deep_copy(rhsWeights_, 0.0f);
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
  auto mm = toHost(solver_->getMassMatrixAcoustic());
  ASSERT_GT(mm.extent(0), 0u);
  float sum = 0.0f;
  for (size_t i = 0; i < mm.extent(0); ++i) sum += mm(i);
  EXPECT_GT(sum, 0.0f);
}

TEST_P(SemSolverAcousticTest, MassMatrixAllPositive) {
  auto mm = toHost(solver_->getMassMatrixAcoustic());
  for (size_t i = 0; i < mm.extent(0); ++i) EXPECT_GT(mm(i), 0.0f) << "mass matrix zero at node " << i;
}

// ======================================================================
// Damping matrix
// ======================================================================
TEST_P(SemSolverAcousticTest, DampingMatrixNonNegative) {
  auto dm = toHost(solver_->getDampingMatrix(0));
  for (size_t i = 0; i < dm.extent(0); ++i) EXPECT_GE(dm(i), 0.0f) << "negative damping at node " << i;
}

TEST_P(SemSolverAcousticTest, DampingMatrixHasBoundaryContribution) {
  auto dm = toHost(solver_->getDampingMatrix(0));
  float sum = 0.0f;
  for (size_t i = 0; i < dm.extent(0); ++i) sum += dm(i);
  EXPECT_GT(sum, 0.0f);
}

// ======================================================================
// resetGlobalVectors
// ======================================================================
TEST_P(SemSolverAcousticTest, ResetGlobalVectorsZerosForceVector) {
  auto& fv = solver_->getForceVector(0);
  setAt(fv, 99.0f, 0);
  solver_->resetGlobalVectors(numNodes_);
  Kokkos::fence();
  EXPECT_FLOAT_EQ(getAt(fv, 0), 0.0f);
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
  auto current = toHost(data.getCurrentField(0));
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(current(i))) << "NaN/Inf at node " << i;
}

TEST_P(SemSolverAcousticTest, ComputeOneStepZeroSourceStaysBounded) {
  Kokkos::deep_copy(pPrev_, 0.0f);
  Kokkos::deep_copy(pCurr_, 0.0f);
  WavefieldAcoustic wf(pPrev_, pCurr_);
  RhsAcoustic rhs(rhsTerm_, rhsElem_, rhsWeights_);
  SEMsolverDataAcoustic data(wf, rhs);
  for (int t = 0; t < kNumSteps; ++t) {
    solver_->computeOneStep(kDt, t, data);
    data.swapWavefields();
  }
  float sum = 0.0f;
  auto current = toHost(data.getCurrentField(0));
  for (int i = 0; i < numNodes_; ++i) sum += std::fabs(current(i));
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
  auto current = toHost(data.getCurrentField(0));
  for (int i = 0; i < numNodes_; ++i) norm += current(i) * current(i);
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
    Kokkos::deep_copy(pPrev_, 0.0f);
    Kokkos::deep_copy(pCurr_, 0.0f);
    setAt(pCurr_, 1.0f, numNodes_ / 2);
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_n");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_n");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_n");
    Kokkos::deep_copy(rhsElem_, 0);
    Kokkos::deep_copy(rhsTerm_, 0.0f);
    Kokkos::deep_copy(rhsWeights_, 0.0f);
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
  auto current = toHost(data.getCurrentField(0));
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(current(i))) << "NaN/Inf at node " << i;
}

TEST_F(SemSolverAcousticOnNodesTest, MassMatrixAllPositive) {
  auto mm = toHost(solver_->getMassMatrixAcoustic());
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
    Kokkos::deep_copy(pPrev_, 0.0f);
    Kokkos::deep_copy(pCurr_, 0.0f);
    setAt(pCurr_, 1.0f, numNodes_ / 2);
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_s");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_s");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_s");
    Kokkos::deep_copy(rhsElem_, 0);
    Kokkos::deep_copy(rhsTerm_, 0.0f);
    Kokkos::deep_copy(rhsWeights_, 0.0f);
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
  auto current = toHost(data.getCurrentField(0));
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(current(i)));
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
    fillFromList(ref, {2.0f * 3.14159f * 1.0f, 2.0f * 3.14159f * 10.0f});
    fillFromList(coeffs, {0.1f, 0.1f});
    solver_->setSLSAttenuation(ref, coeffs);
    solver_->computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);

    numNodes_ = mesh_->getNumberOfNodes();
    constexpr int npp = 8;
    pPrev_ = allocateVector<vectorReal>(numNodes_, "pPrev_att");
    pCurr_ = allocateVector<vectorReal>(numNodes_, "pCurr_att");
    Kokkos::deep_copy(pPrev_, 0.0f);
    Kokkos::deep_copy(pCurr_, 0.0f);
    setAt(pCurr_, 1.0f, numNodes_ / 2);
    rhsTerm_ = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTerm_att");
    rhsElem_ = allocateVector<vectorInt>(1, "rhsElem_att");
    rhsWeights_ = allocateArray2D<arrayReal>(1, npp, "rhsWeights_att");
    Kokkos::deep_copy(rhsElem_, 0);
    Kokkos::deep_copy(rhsTerm_, 0.0f);
    Kokkos::deep_copy(rhsWeights_, 0.0f);
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
  auto current = toHost(data.getCurrentField(0));
  for (int i = 0; i < numNodes_; ++i) EXPECT_TRUE(std::isfinite(current(i))) << "NaN/Inf at node " << i;
}

TEST_F(SemSolverAcousticAttenuationTest, ResetGlobalVectorsZerosAttenuationWorkVectors) {
  auto& fv = solver_->getForceVector(0);
  setAt(fv, 99.0f, 0);
  solver_->resetGlobalVectors(numNodes_);
  Kokkos::fence();
  EXPECT_FLOAT_EQ(getAt(fv, 0), 0.0f);
}

// Covers Solver::setZBoundary default body (non-overriding solver).
TEST_F(SemSolverAcousticAttenuationTest, SetZBoundary_DoesNotThrow) { EXPECT_NO_THROW(solver_->setZBoundary(500.0f)); }

}  // namespace test
}  // namespace fe
}  // namespace solver

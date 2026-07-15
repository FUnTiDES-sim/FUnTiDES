/**
 * @file test_dg_padaptive_solver_acoustic.cc
 * @brief Integration and kernel tests for DGPAdaptiveSolver (acoustic, ORDER_MIN=1, ORDER_MAX=2).
 *
 * The solver is directly instantiated with ModelStruct<float,int,ORDER_MAX> so all
 * template specialisations in dg_padaptive_solver_impl.h are exercised without going
 * through the solver factory. A 2×2×2 element Cartesian mesh (domain 2000³)
 * is built via CartesianStructBuilder; each element is a 1000×1000×1000 cube.
 *
 * Key numerical invariants verified:
 *   - Zero field + zero source: field remains zero (linearity sanity check).
 *   - Uniform field p=1 + zero source: field remains 1 everywhere.
 *     Proof: for p^n=p^{n-1}=1, stiffness=0 (DG consistency on constants),
 *     interface penalty=0 (no jump). Verlet formula then gives
 *     (2M - (M - ½dt·D)) / (M + ½dt·D) = 1 identically.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>

#include "Integrals.h"
#include "cartesian_struct_builder.h"
#include "common_macros.h"
#include "data_type.h"
#include "dg_padaptive_solver_data.h"
#include "dg_padaptive_solver_impl.h"
#include "model_struct.h"

namespace solver {
namespace fe {
namespace test {

// ============================================================
// Type aliases for ORDER_Max=2 DG p-adaptive acoustic solver on structured mesh
// ============================================================

static constexpr int kOrder_max = 2;
static constexpr int kNDof_max = (kOrder_max + 1) * (kOrder_max + 1) * (kOrder_max + 1);  // 27
static constexpr int kOrder_min = 1;
static constexpr int kNDof_min = (kOrder_min + 1) * (kOrder_min + 1) * (kOrder_min + 1);  // 8

using MeshType = model::ModelStruct<float, int, kOrder_max>;
using DGPAdaptiveSolverT = DGPAdaptiveSolver<kOrder_min, kOrder_max, IntegralTypeSelector, IntegralType::MAKUTU, MeshType, false, physicType::kAcoustic>;

// ============================================================
// Fixture
// ============================================================

class DGPAdaptiveSolverAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    model::CartesianStructBuilder<float, int, kOrder_max> builder(2, 2000.0f, 2, 2000.0f, 2, 2000.0f, false, false, 0.0,
                                                              0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, false,
                                                              0.0f);      // default parameters from the builder
    mesh_ = builder.getModel(false);
    nElem_ = mesh_->getNumberOfElements();  // 8 here
    nNode_ = mesh_->getNumberOfNodes();     // 27 here
    solver_.computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
  }

  /**
   * @brief Build a DGSEMsolverData with uniform initial field values
   *        and an optional single point source.
   *
   * @param p_min_CurrVal  Uniform initial value for pnPMinDGCurr.
   * @param p_min_PrevVal  Uniform initial value for pnPMinDGPrev.
   * @param p_max_CurrVal  Uniform initial value for pnPMaxDGCurr.
   * @param p_max_PrevVal  Uniform initial value for pnPMaxDGPrev.
   * @param nSample   Number of time samples in the RHS term (0 = no source).
   * @param srcElem   Source element index (ignored when nSample == 0).
   * @param wavelet   Constant wavelet value for all time samples (ignored when nSample == 0).
   */
  DGPAdaptiveSolverData makeData(float p_min_CurrVal, float p_min_PrevVal, float p_max_CurrVal, float p_max_PrevVal,
                           int nSample = 0, int srcElem = 0, float wavelet = 0.0f) {
    auto p_min_Curr = allocateArray2D<arrayReal>(nElem_, kNDof_min, "pnPMinDGCurr");
    auto p_min_Prev = allocateArray2D<arrayReal>(nElem_, kNDof_min, "pnPMinDGPrev");
    auto p_max_Curr = allocateArray2D<arrayReal>(nElem_, kNDof_max, "pnPMaxDGCurr");
    auto p_max_Prev = allocateArray2D<arrayReal>(nElem_, kNDof_max, "pnPMaxDGPrev");

    for (int e = 0; e < nElem_; ++e) {
      for (int d = 0; d < kNDof_min; ++d) {
        p_min_Curr(e, d) = p_min_CurrVal;
        p_min_Prev(e, d) = p_min_PrevVal;
      }
      for (int d = 0; d < kNDof_max; ++d) {
        p_max_Curr(e, d) = p_max_CurrVal;
        p_max_Prev(e, d) = p_max_PrevVal;
      }
    }
    int const nSrc = (nSample > 0) ? 1 : 0;
    auto rhsTerm_pmin = allocateArray2D<arrayReal>(nSrc, (nSample > 0 ? nSample : 1), "rhsTermPMin");
    auto rhsTerm_pmax = allocateArray2D<arrayReal>(nSrc, (nSample > 0 ? nSample : 1), "rhsTermPMax");
    auto rhsElem = allocateVector<vectorInt>(nSrc, "rhsElem");
    auto rhsWeights_pmin = allocateArray2D<arrayReal>(nSrc, kNDof_min, "rhsWeightsPMin");
    auto rhsWeights_pmax = allocateArray2D<arrayReal>(nSrc, kNDof_max, "rhsWeightsPMax");

    if (nSrc > 0) {
      rhsElem(0) = srcElem;
      for (int t = 0; t < nSample; ++t) {
        rhsTerm_pmin(0, t) = wavelet;
        rhsTerm_pmax(0, t) = wavelet;
      }
      for (int d = 0; d < kNDof_min; ++d) rhsWeights_pmin(0, d) = 1.0f / kNDof_min;
      for (int d = 0; d < kNDof_max; ++d) rhsWeights_pmax(0, d) = 1.0f / kNDof_max;
    }

    DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
    DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem, rhsWeights_pmin, rhsWeights_pmax);
    return DGPAdaptiveSolverData(wavefield, rhs);
  }

  DGPAdaptiveSolverT solver_;
  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  int nElem_, nNode_;
};

// ============================================================
// computeFEInit
// ============================================================

TEST_F(DGPAdaptiveSolverAcousticTest, ComputeFEInit_Succeeds) {
  // SetUp already calls computeFEInit; if it threw, the fixture would fail.
  EXPECT_EQ(nElem_, 8);
}

TEST_F(DGPAdaptiveSolverAcousticTest, ComputeFEInit_IncompatibleMeshThrows) {
  // A ModelStruct<float,int,3> is a different C++ type from ModelStruct<float,int,2>.
  // The dynamic_cast in computeFEInit must fail and throw.
  model::CartesianStructBuilder<float, int, 3> builder2(2, 2000.0f, 2, 2000.0f, 2, 2000.0f, false, false, 0.0, 0.0, 0.0,
                                                        -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, false,
                                                        0.0f);      // default parameters from the builder
  auto mesh2 = builder2.getModel(false);

  DGPAdaptiveSolverT fresh_solver;
  EXPECT_THROW(fresh_solver.computeFEInit(*mesh2, {0.0f, 0.0f, 0.0f}, false, 0.0f), std::runtime_error);
}

// ============================================================
// computeOneStep control flow
// ============================================================

TEST_F(DGPAdaptiveSolverAcousticTest, ComputeOneStep_Serial_DoesNotThrow) {
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f);
  EXPECT_NO_THROW(solver_.computeOneStep(0.001f, 0, data));
}

TEST_F(DGPAdaptiveSolverAcousticTest, ComputeOneStep_DistributedThrows) {
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f);
  data.isDistributed = true;
  EXPECT_THROW(solver_.computeOneStep(0.001f, 0, data), std::runtime_error);
}

// ============================================================
// Numerical invariants
// ============================================================

TEST_F(DGPAdaptiveSolverAcousticTest, ZeroFieldZeroSource_StaysZero) {
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f);
  solver_.computeOneStep(0.001f, 0, data);

  // After one step from p=0 with no source, p must remain 0.
  for (int e = 0; e < nElem_; ++e) {
    for (int d = 0; d < kNDof_min; ++d)
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMinPreviousField(0)(e, d), 0.0f) << "elem=" << e << " dof=" << d;
    for (int d = 0; d < kNDof_max; ++d)
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxPreviousField(0)(e, d), 0.0f) << "elem=" << e << " dof=" << d;
  }
}

TEST_F(DGPAdaptiveSolverAcousticTest, UniformFieldZeroSource_StaysUniform) {
  // p^n = p^{n-1} = 1 everywhere, no source.
  // DG stiffness of a constant field = 0 (consistency).
  // Interface flux and penalty are both 0 (no jump across faces).
  // Verlet: (2M*1 - 0 - (M - ½dt·D)*1) / (M + ½dt·D) = 1 analytically.
  auto data = makeData(1.0f, 1.0f, 1.0f, 1.0f);
  solver_.computeOneStep(0.001f, 0, data);

  for (int e = 0; e < nElem_; ++e) {
    for (int d = 0; d < kNDof_min; ++d)
      EXPECT_NEAR(data.m_wavefield.getPMinPreviousField(0)(e, d), 1.0f, 1e-4f) << "elem=" << e << " dof=" << d;
    for (int d = 0; d < kNDof_max; ++d)
      EXPECT_NEAR(data.m_wavefield.getPMaxPreviousField(0)(e, d), 1.0f, 1e-4f) << "elem=" << e << " dof=" << d;
  }
}

TEST_F(DGPAdaptiveSolverAcousticTest, NonzeroSourcePMin_FieldChangesAtSourceElem) {
  // A non-zero RHS term at element 0 with uniform weights must move the
  // solution away from zero at every DOF of that element.
  constexpr float kWavelet = 1.0f;
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f, /*nSample=*/1, /*srcElem=*/0, kWavelet);
  solver_.computeOneStep(0.001f, 0, data);

  // rhs_elem(0,d) = -kWavelet * (1/kNDof) < 0
  // Verlet (curr=prev=0): new_prev = -dt^2 * rhs / (M + ...) > 0
  for (int d = 0; d < kNDof_min; ++d) EXPECT_GT(data.m_wavefield.getPMinPreviousField(0)(0, d), 0.0f) << "dof=" << d;
}

TEST_F(DGPAdaptiveSolverAcousticTest, NonzeroSourcePMax_FieldChangesAtSourceElem) {
  // A non-zero RHS term at element 0 with uniform weights must move the
  // solution away from zero at every DOF of that element.
  constexpr float kWavelet = 1.0f;
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f, /*nSample=*/1, /*srcElem=*/7, kWavelet);
  solver_.computeOneStep(0.001f, 0, data);

  // rhs_elem(0,d) = -kWavelet * (1/kNDof) < 0
  // Verlet (curr=prev=0): new_prev = -dt^2 * rhs / (M + ...) > 0
  for (int d = 0; d < kNDof_max; ++d) EXPECT_GT(data.m_wavefield.getPMaxPreviousField(0)(7, d), 0.0f) << "dof=" << d;
}

TEST_F(DGPAdaptiveSolverAcousticTest, NoSourceElements_UnaffectedBySource) {
  // Elements other than the source element must remain zero when starting
  // from a zero initial field (the stiffness of the zero field is zero).
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f, /*nSample=*/1, /*srcElem=*/0, 1.0f);  // source in DG
  solver_.computeOneStep(0.001f, 0, data);

  // Elements != 0 should have zero field (the source kernel only writes
  // to element 0; interface flux from zero field is also zero).
  for (int e = 1; e < nElem_; ++e) {
    for (int d = 0; d < kNDof_min; ++d)
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMinPreviousField(0)(e, d), 0.0f) << "elem=" << e << " dof=" << d;
    for (int d = 0; d < kNDof_max; ++d)
      EXPECT_FLOAT_EQ(data.m_wavefield.getPMaxPreviousField(0)(e, d), 0.0f) << "elem=" << e << " dof=" << d;
  }
}

// ============================================================
// Multi-step stability
// ============================================================

TEST_F(DGPAdaptiveSolverAcousticTest, MultipleSteps_NoNanOrInf) {
  constexpr int kNumSteps = 5;
  constexpr float kDt = 0.001f;

  auto p_min_Curr = allocateArray2D<arrayReal>(nElem_, kNDof_min, "pnPMinDGCurr");
  auto p_min_Prev = allocateArray2D<arrayReal>(nElem_, kNDof_min, "pnPMinDGPrev");
  auto p_max_Curr = allocateArray2D<arrayReal>(nNode_, kNDof_max, "pnPMaxDGCurr");
  auto p_max_Prev = allocateArray2D<arrayReal>(nNode_, kNDof_max, "pnPMaxDGPrev");
  auto rhsTerm_pmin = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermPMin");
  auto rhsTerm_pmax = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermPMax");
  auto rhsElem = allocateVector<vectorInt>(1, "rhsElem");
  auto rhsWeights_pmin = allocateArray2D<arrayReal>(1, kNDof_min, "rhsWeightsPMin");
  auto rhsWeights_pmax = allocateArray2D<arrayReal>(1, kNDof_max, "rhsWeightsPMax");

  rhsElem(0) = 0;
  for (int t = 0; t < kNumSteps; ++t) {
    rhsTerm_pmin(0, t) = std::sin(t * 0.1f);
    rhsTerm_pmax(0, t) = std::sin(t * 0.1f);
  }
  for (int d = 0; d < kNDof_min; ++d) rhsWeights_pmin(0, d) = 1.0f / kNDof_min;
  for (int d = 0; d < kNDof_max; ++d) rhsWeights_pmax(0, d) = 1.0f / kNDof_max;

  DGPAdaptiveWavefieldAcoustic wavefield(p_min_Prev, p_min_Curr, p_max_Prev, p_max_Curr);
  DGPAdaptiveRhsAcoustic rhs(rhsTerm_pmin, rhsTerm_pmax, rhsElem, rhsWeights_pmin, rhsWeights_pmax);

  DGPAdaptiveSolverData data(wavefield, rhs);

  for (int t = 0; t < kNumSteps; ++t) {
    ASSERT_NO_THROW(solver_.computeOneStep(kDt, t, data));
    data.swapWavefields();
  }

  for (int e = 0; e < nElem_; ++e) {
    for (int d = 0; d < kNDof_min; ++d)
      EXPECT_TRUE(std::isfinite(data.m_wavefield.getPMinCurrentField(0)(e, d)))
          << "NaN/Inf at elem=" << e << " dof=" << d;
    for (int d = 0; d < kNDof_max; ++d)
      EXPECT_TRUE(std::isfinite(data.m_wavefield.getPMaxCurrentField(0)(e, d)))
          << "NaN/Inf at elem=" << e << " dof=" << d;
  }
}

// ============================================================
// Stub getters: must throw (not implemented for DG p-adaptive)
// ============================================================

TEST_F(DGPAdaptiveSolverAcousticTest, GettersThrowAsExpected) {
  EXPECT_THROW(solver_.getMassMatrixAcoustic(), std::runtime_error);
  EXPECT_THROW(solver_.getMassMatrixElastic(), std::runtime_error);
  EXPECT_THROW(solver_.getDampingMatrix(0), std::runtime_error);
  EXPECT_THROW(solver_.getForceVector(0), std::runtime_error);
}

// ============================================================
// outputSolutionValues
// ============================================================

TEST_F(DGPAdaptiveSolverAcousticTest, OutputSolutionValues_ArrayView_DoesNotCrash) {
  auto field_pmin = allocateArray2D<arrayReal>(nElem_, kNDof_min, "field");
  int t = 0, e = 0;
  EXPECT_NO_THROW(solver_.outputSolutionValues(t, e, field_pmin, "pressurePMin"));
}

}  // namespace test
}  // namespace fe
}  // namespace solver

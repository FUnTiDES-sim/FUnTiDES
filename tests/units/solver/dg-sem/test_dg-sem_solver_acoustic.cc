/**
 * @file test_dg-sem_solver_acoustic.cc
 * @brief Integration and kernel tests for DGSEMsolver (acoustic, ORDER=1).
 *
 * The solver is directly instantiated with ModelStruct<float,int,1> so all
 * template specialisations in dg_solver_impl.h are exercised without going
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
#include "dg-sem_solver_data.h"
#include "dg-sem_solver_impl.h"
#include "model_struct.h"

namespace solver {
namespace fe {
namespace test {

// ============================================================
// Type aliases for ORDER=1 DG-SEM acoustic solver on structured mesh
// ============================================================

static constexpr int kOrder = 1;
static constexpr int kNDof = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);  // 8

using MeshType = model::ModelStruct<float, int, kOrder>;
using IntType = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
using DGSEMSolverT = DGSEMsolver<kOrder, IntType, MeshType, false, physicType::kAcoustic>;

// ============================================================
// Fixture
// ============================================================

class DGSEMsolverAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    model::CartesianStructBuilder<float, int, kOrder> builder(2, 2000.0f, 2, 2000.0f, 2, 2000.0f, false, false, 0.0,
                                                              0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, false,
                                                              0.0f,      // default parameters from the builder
                                                              1000.0f);  // DgSemBoundaryZ
    mesh_ = builder.getModel(false);
    nElem_ = mesh_->getNumberOfElements();  // 8 here
    nNode_ = mesh_->getNumberOfNodes();     // 27 here
    solver_.computeFEInit(*mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
  }

  /**
   * @brief Build a DGSEMsolverData with uniform initial field values
   *        and an optional single point source.
   *
   * @param pCurrVal  Uniform initial value for pnCurr.
   * @param pPrevVal  Uniform initial value for pnPrev.
   * @param nSample   Number of time samples in the RHS term (0 = no source).
   * @param srcElem   Source element index (ignored when nSample == 0).
   * @param wavelet   Constant wavelet value for all time samples (ignored when nSample == 0).
   */
  DGSEMsolverData makeData(float p_dg_CurrVal, float p_dg_PrevVal, float p_sem_CurrVal, float p_sem_PrevVal,
                           int nSample = 0, int srcElem = 0, float wavelet = 0.0f) {
    auto p_dg_Curr = allocateArray2D<arrayReal>(nElem_, kNDof, "pCurrDG");
    auto p_dg_Prev = allocateArray2D<arrayReal>(nElem_, kNDof, "pPrevDG");
    auto p_sem_Curr = allocateVector<vectorReal>(nNode_, "pCurrSEM");
    auto p_sem_Prev = allocateVector<vectorReal>(nNode_, "pPrevSEM");

    for (int e = 0; e < nElem_; ++e)
      for (int d = 0; d < kNDof; ++d) {
        p_dg_Curr(e, d) = p_dg_CurrVal;
        p_dg_Prev(e, d) = p_dg_PrevVal;
      }
    for (int n = 0; n < nNode_; ++n) {
      p_sem_Curr(n) = p_sem_CurrVal;
      p_sem_Prev(n) = p_sem_PrevVal;
    }
    int const nSrc = (nSample > 0) ? 1 : 0;
    auto rhsTerm_dg = allocateArray2D<arrayReal>(nSrc, (nSample > 0 ? nSample : 1), "rhsTermDG");
    auto rhsTerm_sem = allocateArray2D<arrayReal>(nSrc, (nSample > 0 ? nSample : 1), "rhsTermSEM");
    auto rhsElem = allocateVector<vectorInt>(nSrc, "rhsElem");
    auto rhsWeights = allocateArray2D<arrayReal>(nSrc, kNDof, "rhsWeights");

    if (nSrc > 0) {
      rhsElem(0) = srcElem;
      for (int t = 0; t < nSample; ++t) {
        rhsTerm_dg(0, t) = wavelet;
        rhsTerm_sem(0, t) = wavelet;
      }
      for (int d = 0; d < kNDof; ++d) rhsWeights(0, d) = 1.0f / kNDof;
    }

    DGSEMWavefieldAcoustic wavefield(p_dg_Prev, p_dg_Curr, p_sem_Prev, p_sem_Curr);
    DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem, rhsWeights);
    return DGSEMsolverData(wavefield, rhs);
  }

  DGSEMSolverT solver_;
  std::shared_ptr<model::ModelApi<float, int>> mesh_;
  int nElem_, nNode_;
};

// ============================================================
// computeFEInit
// ============================================================

TEST_F(DGSEMsolverAcousticTest, ComputeFEInit_Succeeds) {
  // SetUp already calls computeFEInit; if it threw, the fixture would fail.
  EXPECT_EQ(nElem_, 8);
}

TEST_F(DGSEMsolverAcousticTest, ComputeFEInit_IncompatibleMeshThrows) {
  // A ModelStruct<float,int,2> is a different C++ type from ModelStruct<float,int,1>.
  // The dynamic_cast in computeFEInit must fail and throw.
  model::CartesianStructBuilder<float, int, 2> builder2(2, 2000.0f, 2, 2000.0f, 2, 2000.0f, false, false, 0.0, 0.0, 0.0,
                                                        -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, false,
                                                        0.0f,      // default parameters from the builder
                                                        1000.0f);  // DgSemBoundaryZ
  auto mesh2 = builder2.getModel(false);

  DGSEMSolverT fresh_solver;
  EXPECT_THROW(fresh_solver.computeFEInit(*mesh2, {0.0f, 0.0f, 0.0f}, false, 0.0f), std::runtime_error);
}

// ============================================================
// computeOneStep control flow
// ============================================================

TEST_F(DGSEMsolverAcousticTest, ComputeOneStep_Serial_DoesNotThrow) {
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f);
  EXPECT_NO_THROW(solver_.computeOneStep(0.001f, 0, data));
}

TEST_F(DGSEMsolverAcousticTest, ComputeOneStep_DistributedThrows) {
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f);
  data.isDistributed = true;
  EXPECT_THROW(solver_.computeOneStep(0.001f, 0, data), std::runtime_error);
}

// ============================================================
// Numerical invariants
// ============================================================

TEST_F(DGSEMsolverAcousticTest, ZeroFieldZeroSource_StaysZero) {
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f);
  solver_.computeOneStep(0.001f, 0, data);

  // After one step from p=0 with no source, p must remain 0.
  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < kNDof; ++d)
      EXPECT_FLOAT_EQ(data.m_wavefield.getDGPreviousField(0)(e, d), 0.0f) << "elem=" << e << " dof=" << d;
  for (int n = 0; n < nNode_; ++n) EXPECT_FLOAT_EQ(data.m_wavefield.getSEMPreviousField(0)(n), 0.0f) << "node=" << n;
}

TEST_F(DGSEMsolverAcousticTest, UniformFieldZeroSource_StaysUniform) {
  // p^n = p^{n-1} = 1 everywhere, no source.
  // DG stiffness of a constant field = 0 (consistency).
  // Interface flux and penalty are both 0 (no jump across faces).
  // Verlet: (2M*1 - 0 - (M - ½dt·D)*1) / (M + ½dt·D) = 1 analytically.
  auto data = makeData(1.0f, 1.0f, 1.0f, 1.0f);
  solver_.computeOneStep(0.001f, 0, data);

  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < kNDof; ++d)
      EXPECT_NEAR(data.m_wavefield.getDGPreviousField(0)(e, d), 1.0f, 1e-4f) << "elem=" << e << " dof=" << d;
  for (int n = 0; n < nNode_; ++n) EXPECT_NEAR(data.m_wavefield.getSEMPreviousField(0)(n), 1.0f, 1e-4f) << "node=" << n;
}

TEST_F(DGSEMsolverAcousticTest, NonzeroSourceDG_FieldChangesAtSourceElem) {
  // A non-zero RHS term at element 0 with uniform weights must move the
  // solution away from zero at every DOF of that element.
  constexpr float kWavelet = 1.0f;
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f, /*nSample=*/1, /*srcElem=*/0, kWavelet);
  solver_.computeOneStep(0.001f, 0, data);

  // rhs_elem(0,d) = -kWavelet * (1/kNDof) < 0
  // Verlet (curr=prev=0): new_prev = -dt^2 * rhs / (M + ...) > 0
  for (int d = 0; d < kNDof; ++d) EXPECT_GT(data.m_wavefield.getDGPreviousField(0)(0, d), 0.0f) << "dof=" << d;
}

TEST_F(DGSEMsolverAcousticTest, NonzeroSourceSEM_FieldChangesAtSourceElem) {
  // A non-zero RHS term at element 0 with uniform weights must move the
  // solution away from zero at every DOF of that element.
  constexpr float kWavelet = 1.0f;
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f, /*nSample=*/1, /*srcElem=*/7, kWavelet);
  solver_.computeOneStep(0.001f, 0, data);

  // rhs_elem(0,d) = -kWavelet * (1/kNDof) < 0
  // Verlet (curr=prev=0): new_prev = -dt^2 * rhs / (M + ...) > 0
  EXPECT_GT(data.m_wavefield.getSEMPreviousField(0)(26), 0.0f);
}

TEST_F(DGSEMsolverAcousticTest, NoSourceElements_UnaffectedBySource) {
  // Elements other than the source element must remain zero when starting
  // from a zero initial field (the stiffness of the zero field is zero).
  auto data = makeData(0.0f, 0.0f, 0.0f, 0.0f, /*nSample=*/1, /*srcElem=*/0, 1.0f);  // source in DG
  solver_.computeOneStep(0.001f, 0, data);

  // Elements != 0 should have zero field (the source kernel only writes
  // to element 0; interface flux from zero field is also zero).
  for (int e = 1; e < nElem_; ++e)
    for (int d = 0; d < kNDof; ++d)
      EXPECT_FLOAT_EQ(data.m_wavefield.getDGPreviousField(0)(e, d), 0.0f) << "elem=" << e << " dof=" << d;
  // SEM nodes which are not on the source element interface should have zero field
  EXPECT_FLOAT_EQ(data.m_wavefield.getSEMPreviousField(0)(11), 0.0f) << "node=" << 11;
  for (int n = 14; n < nNode_; ++n)  // test only SEM Nodes
    EXPECT_FLOAT_EQ(data.m_wavefield.getSEMPreviousField(0)(n), 0.0f) << "node=" << n;
}

// ============================================================
// Multi-step stability
// ============================================================

TEST_F(DGSEMsolverAcousticTest, MultipleSteps_NoNanOrInf) {
  constexpr int kNumSteps = 5;
  constexpr float kDt = 0.001f;

  auto p_dg_Curr = allocateArray2D<arrayReal>(nElem_, kNDof, "pCurrDG");
  auto p_dg_Prev = allocateArray2D<arrayReal>(nElem_, kNDof, "pPrevDG");
  auto p_sem_Curr = allocateVector<vectorReal>(nNode_, "pCurrSEM");
  auto p_sem_Prev = allocateVector<vectorReal>(nNode_, "pPrevSEM");
  auto rhsTerm_dg = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermDG");
  auto rhsTerm_sem = allocateArray2D<arrayReal>(1, kNumSteps, "rhsTermSEM");
  auto rhsElem = allocateVector<vectorInt>(1, "rhsElem");
  auto rhsWeights = allocateArray2D<arrayReal>(1, kNDof, "rhsWeights");

  rhsElem(0) = 0;
  for (int t = 0; t < kNumSteps; ++t) {
    rhsTerm_dg(0, t) = std::sin(t * 0.1f);
    rhsTerm_sem(0, t) = std::sin(t * 0.1f);
  }
  for (int d = 0; d < kNDof; ++d) rhsWeights(0, d) = 1.0f / kNDof;

  DGSEMWavefieldAcoustic wavefield(p_dg_Prev, p_dg_Curr, p_sem_Prev, p_sem_Curr);
  DGSEMRhsAcoustic rhs(rhsTerm_dg, rhsTerm_sem, rhsElem, rhsWeights);

  DGSEMsolverData data(wavefield, rhs);

  for (int t = 0; t < kNumSteps; ++t) {
    ASSERT_NO_THROW(solver_.computeOneStep(kDt, t, data));
    data.swapWavefields();
  }

  for (int e = 0; e < nElem_; ++e)
    for (int d = 0; d < kNDof; ++d)
      EXPECT_TRUE(std::isfinite(data.m_wavefield.getDGCurrentField(0)(e, d)))
          << "NaN/Inf at elem=" << e << " dof=" << d;
  for (int n = 0; n < nNode_; ++n)
    EXPECT_TRUE(std::isfinite(data.m_wavefield.getSEMCurrentField(0)(n))) << "NaN/Inf at node=" << n;
}

// ============================================================
// Stub getters: must throw (not implemented for DG-SEM)
// ============================================================

TEST_F(DGSEMsolverAcousticTest, GettersThrowAsExpected) {
  EXPECT_THROW(solver_.getMassMatrixAcoustic(), std::runtime_error);
  EXPECT_THROW(solver_.getMassMatrixElastic(), std::runtime_error);
  EXPECT_THROW(solver_.getDampingMatrix(0), std::runtime_error);
  EXPECT_THROW(solver_.getForceVector(0), std::runtime_error);
}

// ============================================================
// outputSolutionValues
// ============================================================

TEST_F(DGSEMsolverAcousticTest, OutputSolutionValues_ArrayView_DoesNotCrash) {
  auto field_dg = allocateArray2D<arrayReal>(nElem_, kNDof, "field");
  int t = 0, e = 0;
  EXPECT_NO_THROW(solver_.outputSolutionValues(t, e, field_dg, "pressureDG"));
}

TEST_F(DGSEMsolverAcousticTest, OutputSolutionValues_VectorView_DoesNotCrash) {
  auto field_sem = allocateVector<vectorReal>(nNode_, "field_vec");
  int t = 0, e = 0;
  EXPECT_NO_THROW(solver_.outputSolutionValues(t, e, field_sem, "pressureSEM"));
}

}  // namespace test
}  // namespace fe
}  // namespace solver

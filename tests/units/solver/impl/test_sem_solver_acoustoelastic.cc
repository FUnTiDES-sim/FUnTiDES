#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "Integrals.h"
#include "common_macros.h"
#include "data_type.h"
#include "model_struct.h"
#include "sem_solver_acoustoelastic.h"

namespace solver {
namespace fe {
namespace test {

// =============================================================================
// Order wrappers
// =============================================================================

struct Order1 {
  static constexpr int kOrder = 1;
};
struct Order2 {
  static constexpr int kOrder = 2;
};

using AEOrderTypes = ::testing::Types<Order1, Order2>;

// =============================================================================
// Helper: 1×1×2 bilayer mesh, IS_MODEL_ON_NODES = false.
//
//   element 0  (z ∈ [0, dz])    elastic   — vs=1500, rho=2000  → mu ≫ kMuTolerance
//   element 1  (z ∈ [dz, 2dz])  acoustic  — vs=0,    rho=1000  → mu = 0
//
// The (ORDER+1)² nodes at z=dz are shared by both elements → interface.
// =============================================================================

template <int ORDER>
model::ModelStruct<float, int, ORDER> makeBilayerMesh() {
  constexpr int kNumElems = 2;

  model::ModelStructData<float, int> d;
  d.ex_ = 1;
  d.ey_ = 1;
  d.ez_ = 2;
  d.dx_ = 1.0f;
  d.dy_ = 1.0f;
  d.dz_ = 1.0f;
  d.ox_ = 0.0f;
  d.oy_ = 0.0f;
  d.oz_ = 0.0f;
  d.isModelOnNodes_ = false;
  d.isElastic_ = false;

  auto vp = allocateVector<vectorReal>(kNumElems, "vp");
  auto vs = allocateVector<vectorReal>(kNumElems, "vs");
  auto rho = allocateVector<vectorReal>(kNumElems, "rho");
  vp[0] = 3000.0f;
  vs[0] = 1500.0f;
  rho[0] = 2000.0f;  // elastic
  vp[1] = 1500.0f;
  vs[1] = 0.0f;
  rho[1] = 1000.0f;  // acoustic
  d.model_vp_element_ = vp;
  d.model_vs_element_ = vs;
  d.model_rho_element_ = rho;

  return model::ModelStruct<float, int, ORDER>(d);
}

// =============================================================================
// Fixture
// =============================================================================

template <typename OrderWrapper>
class AEsolverOnElemTest : public ::testing::Test {
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  /// DOF per element: (ORDER+1)^3
  static constexpr int kNdof = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr float kDt = 0.001f;
  static constexpr int kNumSamples = 10;

  using Mesh = model::ModelStruct<float, int, kOrder>;
  using Integral = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using Solver = SEMsolverAcoustoElastic<kOrder, Integral, Mesh, false>;

  void SetUp() override {
    mesh_ = makeBilayerMesh<kOrder>();
    solver_.computeFEInit(mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
    nNodes_ = mesh_.getNumberOfNodes();

    p_prev_ = allocateVector<vectorReal>(nNodes_, "pPrev");
    p_curr_ = allocateVector<vectorReal>(nNodes_, "pCurr");
    ux_prev_ = allocateVector<vectorReal>(nNodes_, "uxPrev");
    ux_curr_ = allocateVector<vectorReal>(nNodes_, "uxCurr");
    uy_prev_ = allocateVector<vectorReal>(nNodes_, "uyPrev");
    uy_curr_ = allocateVector<vectorReal>(nNodes_, "uyCurr");
    uz_prev_ = allocateVector<vectorReal>(nNodes_, "uzPrev");
    uz_curr_ = allocateVector<vectorReal>(nNodes_, "uzCurr");
    p_pp_ = allocateVector<vectorReal>(nNodes_, "pPrevPrev");
    ux_pp_ = allocateVector<vectorReal>(nNodes_, "uxPrevPrev");
    uy_pp_ = allocateVector<vectorReal>(nNodes_, "uyPrevPrev");
    uz_pp_ = allocateVector<vectorReal>(nNodes_, "uzPrevPrev");
    zeroWavefields();

    rhs_term_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTerm");
    rhs_termx_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermX");
    rhs_termy_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermY");
    rhs_termz_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermZ");
    rhs_elem_ = allocateVector<vectorInt>(1, "rhsElem");
    rhs_wts_ = allocateArray2D<arrayReal>(1, kNdof, "rhsWts");
    rhs_elem_(0) = 0;
    for (int t = 0; t < kNumSamples; ++t)
      rhs_term_(0, t) = rhs_termx_(0, t) = rhs_termy_(0, t) = rhs_termz_(0, t) = 0.0f;
    for (int k = 0; k < kNdof; ++k) rhs_wts_(0, k) = 0.0f;
  }

  void zeroWavefields() {
    for (int i = 0; i < nNodes_; ++i) {
      p_prev_(i) = p_curr_(i) = 0.0f;
      ux_prev_(i) = ux_curr_(i) = 0.0f;
      uy_prev_(i) = uy_curr_(i) = 0.0f;
      uz_prev_(i) = uz_curr_(i) = 0.0f;
      p_pp_(i) = ux_pp_(i) = uy_pp_(i) = uz_pp_(i) = 0.0f;
    }
  }

  /// Build a DataStruct whose views alias the fixture wavefields/rhs.
  SEMsolverDataAcoustoElastic makeData() const {
    WavefieldAcoustoElastic wf(p_prev_, p_curr_, ux_prev_, ux_curr_, uy_prev_, uy_curr_, uz_prev_, uz_curr_);
    RhsAcoustoElastic rhs(rhs_term_, rhs_elem_, rhs_wts_, rhs_termx_, rhs_termy_, rhs_termz_);
    return SEMsolverDataAcoustoElastic(wf, rhs);
  }

  /// Same, but with the three-buffer wavefield the adjoint mode requires.
  SEMsolverDataAcoustoElastic makeAdjointData() const {
    WavefieldAcoustoElastic wf(p_pp_, p_prev_, p_curr_, ux_pp_, ux_prev_, ux_curr_, uy_pp_, uy_prev_, uy_curr_, uz_pp_,
                               uz_prev_, uz_curr_);
    RhsAcoustoElastic rhs(rhs_term_, rhs_elem_, rhs_wts_, rhs_termx_, rhs_termy_, rhs_termz_);
    return SEMsolverDataAcoustoElastic(wf, rhs);
  }

  Solver solver_;
  Mesh mesh_;
  int nNodes_{0};

  vectorReal p_prev_, p_curr_;
  vectorReal ux_prev_, ux_curr_;
  vectorReal uy_prev_, uy_curr_;
  vectorReal uz_prev_, uz_curr_;
  vectorReal p_pp_, ux_pp_, uy_pp_, uz_pp_;

  arrayReal rhs_term_, rhs_termx_, rhs_termy_, rhs_termz_;
  vectorInt rhs_elem_;
  arrayReal rhs_wts_;
};

TYPED_TEST_SUITE(AEsolverOnElemTest, AEOrderTypes);

// =============================================================================
// TagElements / TagNodes
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, TagElementsClassifiesOneAcousticOneElastic) {
  EXPECT_EQ(this->solver_.getNumAcousticElements(), 1);
  EXPECT_EQ(this->solver_.getNumElasticElements(), 1);
}

TYPED_TEST(AEsolverOnElemTest, TagNodesDetectsCorrectInterfaceNodeCount) {
  // 1×1 shared face between element 0 and element 1: (ORDER+1)^2 interface nodes.
  constexpr int kExpected = (TestFixture::kOrder + 1) * (TestFixture::kOrder + 1);
  EXPECT_EQ(this->solver_.getNumInterfaceNodes(), kExpected);
}

// =============================================================================
// Mass / damping matrices
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, AcousticMassMatrixNonZeroInDomain) {
  auto& m = this->solver_.getMassMatrixAcoustic();
  float total = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) total += m[i];
  EXPECT_GT(total, 0.0f);
}

TYPED_TEST(AEsolverOnElemTest, ElasticMassMatrixNonZeroInDomain) {
  auto& m = this->solver_.getMassMatrixElastic();
  float total = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) total += m[i];
  EXPECT_GT(total, 0.0f);
}

TYPED_TEST(AEsolverOnElemTest, DampingMatrixNonNegativeAbsorbingBC) {
  // computeDampingMatrix assembles absorbing-BC coefficients on boundary faces,
  // not a sponge — entries are always >= 0 and at least one boundary node is > 0.
  for (int c = 0; c < 4; ++c) {
    auto& d = this->solver_.getDampingMatrix(c);
    float total = 0.0f;
    for (int i = 0; i < this->nNodes_; ++i) {
      EXPECT_GE(d[i], 0.0f) << "component " << c << " node " << i;
      total += d[i];
    }
    EXPECT_GT(total, 0.0f) << "no absorbing-BC damping for component " << c;
  }
}

// =============================================================================
// resetGlobalVectors
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, ResetGlobalVectorsZerosForceArrays) {
  for (int c = 0; c < 4; ++c) {
    auto& f = this->solver_.getForceVector(c);
    for (int i = 0; i < this->nNodes_; ++i) f[i] = 1.0f;
  }
  this->solver_.resetGlobalVectors(this->nNodes_);
  FENCE
  for (int c = 0; c < 4; ++c) {
    auto& f = this->solver_.getForceVector(c);
    for (int i = 0; i < this->nNodes_; ++i) EXPECT_FLOAT_EQ(f[i], 0.0f) << "comp=" << c << " node=" << i;
  }
}

// =============================================================================
// Coupling coefficients (tested indirectly via ApplyCouplingAcousticToElastic)
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, CouplingCoeffZNonZeroForHorizontalInterface) {
  // Horizontal z-interface: normal is ±z → coupling modifies uz, not ux/uy.
  // Apply A→E with uniform p=1 from zero uz; at least one node must change.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 1.0f;
  auto data = this->makeData();
  this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data, /*backward=*/false);
  FENCE
  float sum = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) sum += std::fabs(this->uz_prev_(i));
  EXPECT_GT(sum, 0.0f);
}

TYPED_TEST(AEsolverOnElemTest, CouplingCoeffXYZeroForHorizontalInterface) {
  // Horizontal z-interface: cx=cy=0 → ux and uy unchanged.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 1.0f;
  auto data = this->makeData();
  this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data, /*backward=*/false);
  FENCE
  float sum_ux = 0.0f, sum_uy = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) {
    sum_ux += std::fabs(this->ux_prev_(i));
    sum_uy += std::fabs(this->uy_prev_(i));
  }
  EXPECT_NEAR(sum_ux, 0.0f, 1e-6f);
  EXPECT_NEAR(sum_uy, 0.0f, 1e-6f);
}

// =============================================================================
// ApplyCouplingAcousticToElastic
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, CouplingAEZeroPressureNoDisplacement) {
  // p=0 everywhere → coupling formula gives 0 increment.
  for (int i = 0; i < this->nNodes_; ++i) {
    this->p_curr_(i) = 0.0f;
    this->uz_prev_(i) = 5.0f;
  }
  auto data = this->makeData();
  this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data, /*backward=*/false);
  FENCE
  for (int i = 0; i < this->nNodes_; ++i) EXPECT_FLOAT_EQ(this->uz_prev_(i), 5.0f);
}

TYPED_TEST(AEsolverOnElemTest, CouplingAEDtScalingIsBoundedByTheDampedLaw) {
  // u += -dt^2*c*p*taper / (M_e + dt/2*C_e), the same denominator the Verlet
  // update applies to the physical RHS.  Doubling dt therefore multiplies the
  // correction by 4*(M_e + dt/2*C_e)/(M_e + dt*C_e), which is 4 undamped and
  // tends to 2 as the damping dominates.  Every node of this fixture sits on an
  // absorbing face, so only the bracket can be asserted here.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 100.0f;

  {
    auto data = this->makeData();
    this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data, /*backward=*/false);
    FENCE
  }
  float delta1 = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) delta1 += this->uz_prev_(i);

  for (int i = 0; i < this->nNodes_; ++i) this->uz_prev_(i) = 0.0f;
  {
    auto data = this->makeData();
    this->solver_.ApplyCouplingAcousticToElastic(2.0f * this->kDt, data, /*backward=*/false);
    FENCE
  }
  float delta2 = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) delta2 += this->uz_prev_(i);

  if (std::fabs(delta1) <= 1e-10f) GTEST_SKIP() << "coupling coefficient is zero; cannot verify the dt scaling";

  float const ratio = delta2 / delta1;
  EXPECT_GT(ratio, 2.0f);
  EXPECT_LE(ratio, 4.0f + 1e-4f);
}

// =============================================================================
// Backward (adjoint) coupling
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, CouplingAEBackwardCorrectsPrevPrevAndLeavesPrevAlone) {
  // The backward Verlet writes the new level into the prevPrev buffer, so the
  // traction correction must land there.  Before the fix the adjoint solid was
  // never driven and stayed identically zero.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 1.0f;
  auto data = this->makeAdjointData();
  this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data, /*backward=*/true);
  FENCE
  float sum_pp = 0.0f, sum_prev = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) {
    sum_pp += std::fabs(this->uz_pp_(i));
    sum_prev += std::fabs(this->uz_prev_(i));
  }
  EXPECT_GT(sum_pp, 0.0f);
  EXPECT_FLOAT_EQ(sum_prev, 0.0f);
}

TYPED_TEST(AEsolverOnElemTest, CouplingAEBackwardMatchesForwardIncrement) {
  // Same physics, only a different destination buffer: the increment written
  // backward into prevPrev must equal the one written forward into previous.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 3.0f;

  {
    auto data = this->makeData();
    this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data, /*backward=*/false);
    FENCE
  }
  {
    auto data = this->makeAdjointData();
    this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data, /*backward=*/true);
    FENCE
  }
  for (int i = 0; i < this->nNodes_; ++i) EXPECT_FLOAT_EQ(this->uz_pp_(i), this->uz_prev_(i));
}

TYPED_TEST(AEsolverOnElemTest, CouplingEABackwardReadsPrevPrevDisplacement) {
  // Backward, the second difference is u^{n-1} - 2u^n + u^{n+1} with u^{n-1} in
  // prevPrev; a non-zero prevPrev displacement must therefore move the pressure
  // the backward Verlet has just written, also in prevPrev.
  for (int i = 0; i < this->nNodes_; ++i) {
    this->uz_pp_(i) = 1.0f;
    this->uz_curr_(i) = 0.0f;
  }
  auto data = this->makeAdjointData();
  this->solver_.ApplyCouplingElasticToAcoustic(this->kDt, data, /*backward=*/true);
  FENCE
  float sum_p_pp = 0.0f, sum_p_prev = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) {
    sum_p_pp += std::fabs(this->p_pp_(i));
    sum_p_prev += std::fabs(this->p_prev_(i));
  }
  EXPECT_GT(sum_p_pp, 0.0f);
  EXPECT_FLOAT_EQ(sum_p_prev, 0.0f);
}

TYPED_TEST(AEsolverOnElemTest, InterfaceCouplingBackwardExchangesEnergyBothWays) {
  // The composite step must move information in both directions, which is what
  // makes the coupled adjoint carry the fluid residual into the solid.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 2.0f;
  auto data = this->makeAdjointData();
  this->solver_.ApplyInterfaceCoupling(this->kDt, data, /*backward=*/true);
  FENCE
  float sum_u = 0.0f, sum_p = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) {
    sum_u += std::fabs(this->uz_pp_(i));
    sum_p += std::fabs(this->p_pp_(i));
  }
  EXPECT_GT(sum_u, 0.0f) << "fluid pressure did not drive the solid";
  EXPECT_GT(sum_p, 0.0f) << "solid acceleration did not react on the fluid";
}

// =============================================================================
// ApplyCouplingElasticToAcoustic
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, CouplingEAZeroDisplacementNoEffect) {
  // u^{n+1}=u^n=u^{n-1}=0 → FD(u)=0 → p unchanged.
  auto data = this->makeData();
  this->solver_.ApplyCouplingElasticToAcoustic(this->kDt, data, /*backward=*/false);
  FENCE
  for (int i = 0; i < this->nNodes_; ++i) EXPECT_FLOAT_EQ(this->p_prev_(i), 0.0f);
}

TYPED_TEST(AEsolverOnElemTest, CouplingEANonZeroAccelerationChangesP) {
  // uz_prev (= u^{n+1}) = 1, uz_curr (= u^n) = 0, nm1 = 0 (initialized by TagNodes)
  // FD_z = 1 - 0 + 0 = 1 ≠ 0 → pressure at interface nodes must change.
  for (int i = 0; i < this->nNodes_; ++i) {
    this->uz_prev_(i) = 1.0f;
    this->uz_curr_(i) = 0.0f;
  }
  auto data = this->makeData();
  this->solver_.ApplyCouplingElasticToAcoustic(this->kDt, data, /*backward=*/false);
  FENCE
  float sum_p = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) sum_p += std::fabs(this->p_prev_(i));
  EXPECT_GT(sum_p, 0.0f);
}

// =============================================================================
// computeOneStep
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, ComputeOneStepDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(this->solver_.computeOneStep(this->kDt, 0, data));
}

TYPED_TEST(AEsolverOnElemTest, ComputeOneStepZeroFieldsStayNearZero) {
  // No source, zero IC → fields remain near zero after one step.
  auto data = this->makeData();
  this->solver_.computeOneStep(this->kDt, 0, data);
  FENCE
  float max_p = 0.0f, max_uz = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) {
    max_p = std::fmax(max_p, std::fabs(this->p_prev_(i)));
    max_uz = std::fmax(max_uz, std::fabs(this->uz_prev_(i)));
  }
  EXPECT_NEAR(max_p, 0.0f, 1e-10f);
  EXPECT_NEAR(max_uz, 0.0f, 1e-10f);
}

TYPED_TEST(AEsolverOnElemTest, ComputeOneStepRunsMultipleSteps) {
  auto data = this->makeData();
  for (int t = 0; t < this->kNumSamples; ++t) EXPECT_NO_THROW(this->solver_.computeOneStep(this->kDt, t, data));
}

// =============================================================================
// computeForces / updateSolutionForward
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, ComputeForcesDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(this->solver_.computeForces(this->kDt, 0, data));
}

TYPED_TEST(AEsolverOnElemTest, updateSolutionForwardWith2BuffersWorks) {
  auto data = this->makeData();
  EXPECT_NO_THROW(this->solver_.updateSolutionForward(this->kDt, data));
}

TYPED_TEST(AEsolverOnElemTest, updateSolutionForwardWith3BuffersThrows) {
  // Allocate prevprev buffers for 3-buffer mode
  auto p_prevprev = allocateVector<vectorReal>(this->nNodes_, "pPrevPrev");
  auto ux_prevprev = allocateVector<vectorReal>(this->nNodes_, "uxPrevPrev");
  auto uy_prevprev = allocateVector<vectorReal>(this->nNodes_, "uyPrevPrev");
  auto uz_prevprev = allocateVector<vectorReal>(this->nNodes_, "uzPrevPrev");
  for (int i = 0; i < this->nNodes_; ++i) {
    p_prevprev(i) = ux_prevprev(i) = uy_prevprev(i) = uz_prevprev(i) = 0.0f;
  }
  WavefieldAcoustoElastic wf(p_prevprev, this->p_prev_, this->p_curr_, ux_prevprev, this->ux_prev_, this->ux_curr_,
                             uy_prevprev, this->uy_prev_, this->uy_curr_, uz_prevprev, this->uz_prev_, this->uz_curr_);
  RhsAcoustoElastic rhs(this->rhs_term_, this->rhs_elem_, this->rhs_wts_, this->rhs_termx_, this->rhs_termy_,
                        this->rhs_termz_);
  SEMsolverDataAcoustoElastic data(wf, rhs);
  EXPECT_THROW(this->solver_.updateSolutionForward(this->kDt, data), std::runtime_error);
}

TYPED_TEST(AEsolverOnElemTest, updateSolutionBackwardWith2BuffersThrows) {
  auto data = this->makeData();
  EXPECT_THROW(this->solver_.updateSolutionBackward(this->kDt, data), std::runtime_error);
}

TYPED_TEST(AEsolverOnElemTest, updateSolutionBackwardWith3BuffersWorks) {
  // Allocate prevprev buffers for 3-buffer mode
  auto p_prevprev = allocateVector<vectorReal>(this->nNodes_, "pPrevPrev");
  auto ux_prevprev = allocateVector<vectorReal>(this->nNodes_, "uxPrevPrev");
  auto uy_prevprev = allocateVector<vectorReal>(this->nNodes_, "uyPrevPrev");
  auto uz_prevprev = allocateVector<vectorReal>(this->nNodes_, "uzPrevPrev");
  for (int i = 0; i < this->nNodes_; ++i) {
    p_prevprev(i) = ux_prevprev(i) = uy_prevprev(i) = uz_prevprev(i) = 0.0f;
  }
  WavefieldAcoustoElastic wf(p_prevprev, this->p_prev_, this->p_curr_, ux_prevprev, this->ux_prev_, this->ux_curr_,
                             uy_prevprev, this->uy_prev_, this->uy_curr_, uz_prevprev, this->uz_prev_, this->uz_curr_);
  RhsAcoustoElastic rhs(this->rhs_term_, this->rhs_elem_, this->rhs_wts_, this->rhs_termx_, this->rhs_termy_,
                        this->rhs_termz_);
  SEMsolverDataAcoustoElastic data(wf, rhs);
  EXPECT_NO_THROW(this->solver_.updateSolutionBackward(this->kDt, data));
}

// =============================================================================
// Time reversibility of the coupled step
//
// The coupled Verlet is time-symmetric, so with no dissipation a backward step
// taken from (u^n, u^{n+1}) must return u^{n-1} identically:
//
//   forward :  u^{n+1} = 2u^n - u^{n-1} - dt^2 M^-1 K u^n + c(p^n)
//   backward:  u^{n-1} = 2u^n - u^{n+1} - dt^2 M^-1 K u^n + c(p^n)
//
// The same stiffness term and the same interface correction appear in both, so
// substituting one into the other gives back the original level.  The fluid
// obeys the same identity because the second difference the elastic-to-acoustic
// coupling feeds on, u^new - 2u^n + u^other, is symmetric in the two
// neighbouring levels: forward it is u^{n+1} - 2u^n + u^{n-1}, backward it is
// u^{n-1} - 2u^n + u^{n+1}.  Reading the wrong buffer in either direction, or
// dropping the coupling from the backward step, breaks the identity.
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, ForwardThenBackwardStepIsTimeReversible) {
  // The two sub-solvers step the nodes of their own node list, which is private,
  // and the elastic mass matrix is non-zero even on pure fluid nodes, so it
  // cannot serve as a mask.  Instead the prevprev buffers are stamped with a
  // value the scheme cannot produce; whatever still carries it afterwards was
  // never stepped and is not part of the identity being tested.
  constexpr float kUntouched = -98765.0f;

  // Every node of this fixture sits on an absorbing face, and dissipation is
  // not reversible, so the damping has to go before the identity can hold.
  for (int c = 0; c < 4; ++c) {
    auto damping = this->solver_.getDampingMatrix(c);
    for (int i = 0; i < this->nNodes_; ++i) damping(i) = 0.0f;
  }

  // computeForces is deliberately not called: with a zero stiffness term the
  // bulk update collapses to u^{n+1} = 2u^n - u^{n-1} and what remains on top of
  // it is exactly the interface coupling, which is what this test is about.
  // With the stiffness included the coupling sits orders of magnitude below the
  // bulk term and float noise would hide any error in it.  For the same reason
  // dt is taken large: without stiffness there is no CFL limit, and a large dt
  // lifts the O(dt^2) coupling correction clear of roundoff.
  constexpr float kDtRev = 1.0f;

  for (int i = 0; i < this->nNodes_; ++i) {
    float const s = static_cast<float>(i + 1);
    this->p_prev_(i) = 0.5f * std::sin(0.3f * s);
    this->p_curr_(i) = 0.5f * std::sin(0.3f * s + 0.7f);
    this->ux_prev_(i) = 0.2f * std::cos(0.4f * s);
    this->ux_curr_(i) = 0.2f * std::cos(0.4f * s + 0.5f);
    this->uy_prev_(i) = 0.3f * std::sin(0.2f * s);
    this->uy_curr_(i) = 0.3f * std::sin(0.2f * s + 0.6f);
    this->uz_prev_(i) = 0.4f * std::cos(0.5f * s);
    this->uz_curr_(i) = 0.4f * std::cos(0.5f * s + 0.4f);
    this->p_pp_(i) = this->ux_pp_(i) = this->uy_pp_(i) = this->uz_pp_(i) = kUntouched;
  }

  // u^{n-1}: what the backward step has to reconstruct.
  std::vector<float> p_ref(this->nNodes_), ux_ref(this->nNodes_), uy_ref(this->nNodes_), uz_ref(this->nNodes_);
  std::vector<float> p_n(this->nNodes_), uz_n(this->nNodes_);
  for (int i = 0; i < this->nNodes_; ++i) {
    p_ref[i] = this->p_prev_(i);
    ux_ref[i] = this->ux_prev_(i);
    uy_ref[i] = this->uy_prev_(i);
    uz_ref[i] = this->uz_prev_(i);
    p_n[i] = this->p_curr_(i);
    uz_n[i] = this->uz_curr_(i);
  }

  // One forward step; the previous buffers are overwritten in place with u^{n+1}.
  auto fwd = this->makeData();
  this->solver_.updateSolutionForward(kDtRev, fwd);

  // One backward step from (u^n, u^{n+1}); it writes into the prevprev buffers
  // and leaves the current and previous ones alone, so the forward result is
  // still available afterwards.
  auto bwd = this->makeAdjointData();
  this->solver_.updateSolutionBackward(kDtRev, bwd);

  // Measure the coupling correction the forward step actually applied, as its
  // departure from the uncoupled u^{n+1} = 2u^n - u^{n-1}.  It calibrates the
  // tolerance below and guarantees the assertions are not vacuous.
  float coupling_u = 0.0f, coupling_p = 0.0f;
  int checked_el = 0, checked_ac = 0;
  for (int i = 0; i < this->nNodes_; ++i) {
    if (this->uz_pp_(i) != kUntouched) {
      ++checked_el;
      float const d = std::fabs(this->uz_prev_(i) - (2.0f * uz_n[i] - uz_ref[i]));
      if (d > coupling_u) coupling_u = d;
    }
    if (this->p_pp_(i) != kUntouched) {
      ++checked_ac;
      float const d = std::fabs(this->p_prev_(i) - (2.0f * p_n[i] - p_ref[i]));
      if (d > coupling_p) coupling_p = d;
    }
  }
  EXPECT_GT(checked_el, 0) << "the backward step never stepped the solid";
  EXPECT_GT(checked_ac, 0) << "the backward step never stepped the fluid";
  EXPECT_GT(coupling_u, 0.0f) << "the fluid pressure never reached the solid, so this test proves nothing";
  EXPECT_GT(coupling_p, 0.0f) << "the solid acceleration never reached the fluid, so this test proves nothing";

  // Demand the reconstruction error stay well below the coupling correction
  // itself: that is what makes this sensitive to the backward coupling rather
  // than only to the bulk update.
  float const tol_u = 1e-3f * coupling_u;
  float const tol_p = 1e-3f * coupling_p;
  for (int i = 0; i < this->nNodes_; ++i) {
    if (this->uz_pp_(i) != kUntouched) {
      EXPECT_NEAR(this->ux_pp_(i), ux_ref[i], tol_u) << "ux not recovered at node " << i;
      EXPECT_NEAR(this->uy_pp_(i), uy_ref[i], tol_u) << "uy not recovered at node " << i;
      EXPECT_NEAR(this->uz_pp_(i), uz_ref[i], tol_u) << "uz not recovered at node " << i;
    }
    if (this->p_pp_(i) != kUntouched) {
      EXPECT_NEAR(this->p_pp_(i), p_ref[i], tol_p) << "pressure not recovered at node " << i;
    }
  }
}

// =============================================================================
// Misc methods
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, InitSpongeValuesDoesNotCrash) { EXPECT_NO_THROW(this->solver_.initSpongeValues()); }

TYPED_TEST(AEsolverOnElemTest, InitFEarraysDoesNotCrash) { EXPECT_NO_THROW(this->solver_.initFEarrays()); }

TYPED_TEST(AEsolverOnElemTest, OutputSolutionValuesDoesNotCrash) {
  int e = 1;  // acoustic element
  EXPECT_NO_THROW(this->solver_.outputSolutionValues(0, e, this->p_prev_, "pressure"));
}

TYPED_TEST(AEsolverOnElemTest, GetNumComponentsReturnsFour) { EXPECT_EQ(this->solver_.getNumComponents(), 4); }

TYPED_TEST(AEsolverOnElemTest, SetAnisotropyTypeDoesNotCrash) {
  EXPECT_NO_THROW(this->solver_.setAnisotropyType(model::AnisotropyType::kIso));
}

TYPED_TEST(AEsolverOnElemTest, SetSLSAttenuationDoesNotCrash) {
  auto ref = allocateVector<vectorReal>(2, "sls_ref_ae");
  ref(0) = 6.28f;
  ref(1) = 62.8f;
  EXPECT_NO_THROW(this->solver_.setSLSAttenuation(ref));
}

TYPED_TEST(AEsolverOnElemTest, DataStructPrintDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(data.print());
}

TYPED_TEST(AEsolverOnElemTest, DataStructSwapWavefieldsDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(data.swapWavefields());
}

// =============================================================================
// Bilayer mesh helper — IS_MODEL_ON_NODES = true
//
//   Same 1×1×2 topology as makeBilayerMesh, but properties stored on nodes.
//   Interface nodes (shared z-layer at iz=ORDER) carry fluid props (>= conv.)
// =============================================================================

// vs_on_interface is stored on the shared plane only: the element type is
// probed at the element centre, so putting shear anywhere else in the fluid
// would turn that element elastic instead of exercising the interface.
template <int ORDER>
model::ModelStruct<float, int, ORDER> makeBilayerMeshOnNodes(float vs_on_interface = 0.0f) {
  int const nx = ORDER + 1, ny = ORDER + 1, nz = 2 * ORDER + 1;
  int const nNodes = nx * ny * nz;

  model::ModelStructData<float, int> d;
  d.ex_ = 1;
  d.ey_ = 1;
  d.ez_ = 2;
  d.dx_ = 1.0f;
  d.dy_ = 1.0f;
  d.dz_ = 1.0f;
  d.ox_ = 0.0f;
  d.oy_ = 0.0f;
  d.oz_ = 0.0f;
  d.isModelOnNodes_ = true;
  d.isElastic_ = false;

  auto vp_node = allocateVector<vectorReal>(nNodes, "vp_node_n");
  auto vs_node = allocateVector<vectorReal>(nNodes, "vs_node_n");
  auto rho_node = allocateVector<vectorReal>(nNodes, "rho_node_n");

  for (int iz = 0; iz < nz; ++iz)
    for (int iy = 0; iy < ny; ++iy)
      for (int ix = 0; ix < nx; ++ix) {
        int const n = ix + nx * (iy + ny * iz);
        if (iz < ORDER) {
          vp_node(n) = 3000.0f;
          vs_node(n) = 1500.0f;
          rho_node(n) = 2000.0f;
        } else {
          vp_node(n) = 1500.0f;
          vs_node(n) = (iz == ORDER) ? vs_on_interface : 0.0f;
          rho_node(n) = 1000.0f;
        }
      }

  d.model_vp_node_ = vp_node;
  d.model_vs_node_ = vs_node;
  d.model_rho_node_ = rho_node;

  return model::ModelStruct<float, int, ORDER>(d);
}

// =============================================================================
// Fixture — IS_MODEL_ON_NODES = true
// =============================================================================

template <typename OrderWrapper>
class AEsolverOnNodesTest : public ::testing::Test {
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNdof = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr float kDt = 0.001f;
  static constexpr int kNumSamples = 10;

  using Mesh = model::ModelStruct<float, int, kOrder>;
  using Integral = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using Solver = SEMsolverAcoustoElastic<kOrder, Integral, Mesh, true>;

  void SetUp() override {
    mesh_ = makeBilayerMeshOnNodes<kOrder>();
    solver_.computeFEInit(mesh_, {0.0f, 0.0f, 0.0f}, false, 0.0f);
    nNodes_ = mesh_.getNumberOfNodes();

    p_prev_ = allocateVector<vectorReal>(nNodes_, "pPrev_n");
    p_curr_ = allocateVector<vectorReal>(nNodes_, "pCurr_n");
    ux_prev_ = allocateVector<vectorReal>(nNodes_, "uxPrev_n");
    ux_curr_ = allocateVector<vectorReal>(nNodes_, "uxCurr_n");
    uy_prev_ = allocateVector<vectorReal>(nNodes_, "uyPrev_n");
    uy_curr_ = allocateVector<vectorReal>(nNodes_, "uyCurr_n");
    uz_prev_ = allocateVector<vectorReal>(nNodes_, "uzPrev_n");
    uz_curr_ = allocateVector<vectorReal>(nNodes_, "uzCurr_n");
    for (int i = 0; i < nNodes_; ++i) {
      p_prev_(i) = p_curr_(i) = 0.0f;
      ux_prev_(i) = ux_curr_(i) = 0.0f;
      uy_prev_(i) = uy_curr_(i) = 0.0f;
      uz_prev_(i) = uz_curr_(i) = 0.0f;
    }

    rhs_term_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTerm_n");
    rhs_termx_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermX_n");
    rhs_termy_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermY_n");
    rhs_termz_ = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermZ_n");
    rhs_elem_ = allocateVector<vectorInt>(1, "rhsElem_n");
    rhs_wts_ = allocateArray2D<arrayReal>(1, kNdof, "rhsWts_n");
    rhs_elem_(0) = 0;
    for (int t = 0; t < kNumSamples; ++t)
      rhs_term_(0, t) = rhs_termx_(0, t) = rhs_termy_(0, t) = rhs_termz_(0, t) = 0.0f;
    for (int k = 0; k < kNdof; ++k) rhs_wts_(0, k) = 0.0f;
  }

  SEMsolverDataAcoustoElastic makeData() const {
    WavefieldAcoustoElastic wf(p_prev_, p_curr_, ux_prev_, ux_curr_, uy_prev_, uy_curr_, uz_prev_, uz_curr_);
    RhsAcoustoElastic rhs(rhs_term_, rhs_elem_, rhs_wts_, rhs_termx_, rhs_termy_, rhs_termz_);
    return SEMsolverDataAcoustoElastic(wf, rhs);
  }

  Solver solver_;
  Mesh mesh_;
  int nNodes_{0};

  vectorReal p_prev_, p_curr_;
  vectorReal ux_prev_, ux_curr_;
  vectorReal uy_prev_, uy_curr_;
  vectorReal uz_prev_, uz_curr_;

  arrayReal rhs_term_, rhs_termx_, rhs_termy_, rhs_termz_;
  vectorInt rhs_elem_;
  arrayReal rhs_wts_;
};

TYPED_TEST_SUITE(AEsolverOnNodesTest, AEOrderTypes);

TYPED_TEST(AEsolverOnNodesTest, TagElementsClassifiesOneAcousticOneElastic) {
  EXPECT_EQ(this->solver_.getNumAcousticElements(), 1);
  EXPECT_EQ(this->solver_.getNumElasticElements(), 1);
}

TYPED_TEST(AEsolverOnNodesTest, TagNodesDetectsCorrectInterfaceNodeCount) {
  constexpr int kExpected = (TestFixture::kOrder + 1) * (TestFixture::kOrder + 1);
  EXPECT_EQ(this->solver_.getNumInterfaceNodes(), kExpected);
}

TYPED_TEST(AEsolverOnNodesTest, ComputeOneStepDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(this->solver_.computeOneStep(this->kDt, 0, data));
}

TYPED_TEST(AEsolverOnNodesTest, ComputeOneStepZeroFieldsStayNearZero) {
  auto data = this->makeData();
  this->solver_.computeOneStep(this->kDt, 0, data);
  FENCE
  float max_p = 0.0f, max_uz = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) {
    max_p = std::fmax(max_p, std::fabs(this->p_prev_(i)));
    max_uz = std::fmax(max_uz, std::fabs(this->uz_prev_(i)));
  }
  EXPECT_NEAR(max_p, 0.0f, 1e-10f);
  EXPECT_NEAR(max_uz, 0.0f, 1e-10f);
}

TYPED_TEST(AEsolverOnNodesTest, ComputeForcesDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(this->solver_.computeForces(this->kDt, 0, data));
}

TYPED_TEST(AEsolverOnNodesTest, MassMatricesNonZero) {
  float total_a = 0.0f, total_e = 0.0f;
  auto& ma = this->solver_.getMassMatrixAcoustic();
  auto& me = this->solver_.getMassMatrixElastic();
  for (int i = 0; i < this->nNodes_; ++i) {
    total_a += ma[i];
    total_e += me[i];
  }
  EXPECT_GT(total_a, 0.0f);
  EXPECT_GT(total_e, 0.0f);
}

// =============================================================================
// Interface property convention and coupling scheme
// =============================================================================

namespace {

/// Run a few coupled steps from a seeded pressure pulse and return the
/// resulting pressure and vertical displacement.
template <int ORDER, typename SolverT, typename MeshT>
void runSeededSteps(float vs_on_interface, std::vector<float>& p_out, std::vector<float>& uz_out, float dt = 2.0e-4f) {
  constexpr int kNdof = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);
  constexpr int kNumSamples = 4;

  MeshT mesh = makeBilayerMeshOnNodes<ORDER>(vs_on_interface);
  SolverT solver;
  solver.setInterfacePropertyConvention(utils::enums::interfacePropertyConvention::kSharedOnInterfaceNodes);
  solver.computeFEInit(mesh, {0.0f, 0.0f, 0.0f}, false, 0.0f);

  int const n = mesh.getNumberOfNodes();
  auto p_prev = allocateVector<vectorReal>(n, "pPrev_c");
  auto p_curr = allocateVector<vectorReal>(n, "pCurr_c");
  auto ux_prev = allocateVector<vectorReal>(n, "uxPrev_c");
  auto ux_curr = allocateVector<vectorReal>(n, "uxCurr_c");
  auto uy_prev = allocateVector<vectorReal>(n, "uyPrev_c");
  auto uy_curr = allocateVector<vectorReal>(n, "uyCurr_c");
  auto uz_prev = allocateVector<vectorReal>(n, "uzPrev_c");
  auto uz_curr = allocateVector<vectorReal>(n, "uzCurr_c");
  for (int i = 0; i < n; ++i) {
    // A smooth, non-symmetric seed so every coupling term is exercised.
    float const s = static_cast<float>(i % 7) / 7.0f;
    p_prev(i) = p_curr(i) = s;
    ux_prev(i) = ux_curr(i) = 0.5f * s;
    uy_prev(i) = uy_curr(i) = -0.25f * s;
    uz_prev(i) = uz_curr(i) = 0.75f * s;
  }

  auto rhs_term = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTerm_c");
  auto rhs_termx = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermX_c");
  auto rhs_termy = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermY_c");
  auto rhs_termz = allocateArray2D<arrayReal>(1, kNumSamples, "rhsTermZ_c");
  auto rhs_elem = allocateVector<vectorInt>(1, "rhsElem_c");
  auto rhs_wts = allocateArray2D<arrayReal>(1, kNdof, "rhsWts_c");
  rhs_elem(0) = 0;
  for (int t = 0; t < kNumSamples; ++t) rhs_term(0, t) = rhs_termx(0, t) = rhs_termy(0, t) = rhs_termz(0, t) = 0.0f;
  for (int k = 0; k < kNdof; ++k) rhs_wts(0, k) = 0.0f;

  WavefieldAcoustoElastic wf(p_prev, p_curr, ux_prev, ux_curr, uy_prev, uy_curr, uz_prev, uz_curr);
  RhsAcoustoElastic rhs(rhs_term, rhs_elem, rhs_wts, rhs_termx, rhs_termy, rhs_termz);
  SEMsolverDataAcoustoElastic data(wf, rhs);

  solver.computeForces(dt, 0, data);
  FENCE
  solver.updateSolutionForward(dt, data);
  FENCE

  p_out.assign(n, 0.0f);
  uz_out.assign(n, 0.0f);
  for (int i = 0; i < n; ++i) {
    p_out[i] = p_prev(i);
    uz_out[i] = uz_prev(i);
  }
}

}  // namespace

// Under the shared convention the interface node holds one state used by both
// sides, so the shear speed stored there is what the elastic kernel sees.  It
// must reach the solution: zeroing it removes the tangential stiffness that
// holds the interface together and the coupled run diverges.
TYPED_TEST(AEsolverOnNodesTest, SharedConventionPassesInterfaceShearToTheElasticKernel) {
  using MeshT = typename TestFixture::Mesh;
  using SolverT = typename TestFixture::Solver;
  constexpr int kOrder = TestFixture::kOrder;

  if constexpr (kOrder < 2) {
    // At order 1 the element centre is a corner, i.e. an interface node, so
    // the shear cannot be varied without also changing the element type.
    GTEST_SKIP() << "needs an element centre that is not on the interface";
  } else {
    std::vector<float> p_no_shear, uz_no_shear, p_shear, uz_shear;
    runSeededSteps<kOrder, SolverT, MeshT>(0.0f, p_no_shear, uz_no_shear);
    runSeededSteps<kOrder, SolverT, MeshT>(800.0f, p_shear, uz_shear);

    ASSERT_EQ(p_no_shear.size(), p_shear.size());
    double diff = 0.0;
    for (size_t i = 0; i < uz_no_shear.size(); ++i) diff += std::fabs(uz_shear[i] - uz_no_shear[i]);
    EXPECT_GT(diff, 0.0) << "the interface shear speed never reached the elastic kernel";
  }
}

}  // namespace test
}  // namespace fe
}  // namespace solver

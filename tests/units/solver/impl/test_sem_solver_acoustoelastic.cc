#include <gtest/gtest.h>

#include <array>
#include <cmath>

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
    }
  }

  /// Build a DataStruct whose views alias the fixture wavefields/rhs.
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
  this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data);
  FENCE
  float sum = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) sum += std::fabs(this->uz_prev_(i));
  EXPECT_GT(sum, 0.0f);
}

TYPED_TEST(AEsolverOnElemTest, CouplingCoeffXYZeroForHorizontalInterface) {
  // Horizontal z-interface: cx=cy=0 → ux and uy unchanged.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 1.0f;
  auto data = this->makeData();
  this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data);
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
  this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data);
  FENCE
  for (int i = 0; i < this->nNodes_; ++i) EXPECT_FLOAT_EQ(this->uz_prev_(i), 5.0f);
}

TYPED_TEST(AEsolverOnElemTest, CouplingAEScalesWithDt2) {
  // Δu ∝ dt² (formula: u += dt²·c·(-p)/M_e). Doubling dt → 4× change.
  for (int i = 0; i < this->nNodes_; ++i) this->p_curr_(i) = 100.0f;

  {
    auto data = this->makeData();
    this->solver_.ApplyCouplingAcousticToElastic(this->kDt, data);
    FENCE
  }
  float delta1 = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) delta1 += this->uz_prev_(i);

  for (int i = 0; i < this->nNodes_; ++i) this->uz_prev_(i) = 0.0f;
  {
    auto data = this->makeData();
    this->solver_.ApplyCouplingAcousticToElastic(2.0f * this->kDt, data);
    FENCE
  }
  float delta2 = 0.0f;
  for (int i = 0; i < this->nNodes_; ++i) delta2 += this->uz_prev_(i);

  if (std::fabs(delta1) > 1e-10f)
    EXPECT_NEAR(delta2 / delta1, 4.0f, 1e-4f);
  else
    GTEST_SKIP() << "coupling coefficient is zero; cannot verify dt² scaling";
}

// =============================================================================
// ApplyCouplingElasticToAcoustic
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, CouplingEAZeroDisplacementNoEffect) {
  // u^{n+1}=u^n=u^{n-1}=0 → FD(u)=0 → p unchanged.
  auto data = this->makeData();
  this->solver_.ApplyCouplingElasticToAcoustic(data);
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
  this->solver_.ApplyCouplingElasticToAcoustic(data);
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
// computeForces / updateSolution
// =============================================================================

TYPED_TEST(AEsolverOnElemTest, ComputeForcesDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(this->solver_.computeForces(this->kDt, 0, data));
}

TYPED_TEST(AEsolverOnElemTest, UpdateSolutionDoesNotCrash) {
  auto data = this->makeData();
  EXPECT_NO_THROW(this->solver_.updateSolution(this->kDt, data));
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

template <int ORDER>
model::ModelStruct<float, int, ORDER> makeBilayerMeshOnNodes() {
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
          vs_node(n) = 0.0f;
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

}  // namespace test
}  // namespace fe
}  // namespace solver

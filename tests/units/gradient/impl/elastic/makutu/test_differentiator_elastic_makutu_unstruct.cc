#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "Integrals.h"
#include "data_type.h"
#include "device_view_access.h"
#include "differentiator_data_elastic.h"
#include "differentiator_elastic.h"
#include "model_unstruct.h"

namespace gradient {
namespace test {

using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

// =============================================================================
// Order wrappers
// =============================================================================

struct Order1U {
  static constexpr int kOrder = 1;
};
struct Order2U {
  static constexpr int kOrder = 2;
};
struct Order3U {
  static constexpr int kOrder = 3;
};

#ifndef MAX_DIFFERENTIATOR_ELASTIC_ORDER
#define MAX_DIFFERENTIATOR_ELASTIC_ORDER 3
#endif

// Création dynamique de la liste de tests selon l'ordre max
#if MAX_DIFFERENTIATOR_ELASTIC_ORDER == 1
using OrderTypesU = ::testing::Types<Order1U>;
#elif MAX_DIFFERENTIATOR_ELASTIC_ORDER == 2
using OrderTypesU = ::testing::Types<Order1U, Order2U>;
#elif MAX_DIFFERENTIATOR_ELASTIC_ORDER == 3
using OrderTypesU = ::testing::Types<Order1U, Order2U, Order3U>;
#endif

// =============================================================================
// Helper: 1-element unstructured unit cube mesh [0,1]³
// =============================================================================

template <int ORDER>
static model::ModelUnstruct<float, int> makeUnstructMesh1x1x1() {
  constexpr int npe = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  model::ModelUnstructData<float, int> data;
  data.order_ = ORDER;
  data.n_element_ = 1;
  data.n_node_ = npe;
  data.lx_ = data.ly_ = data.lz_ = 1.0f;
  data.isModelOnNodes_ = false;
  data.isElastic_ = true;

  data.global_node_index_ = allocateArray2D<dualArrayInt>(1, npe, "gni");
  data.nodes_coords_x_ = allocateVector<dualVectorReal>(npe, "cx");
  data.nodes_coords_y_ = allocateVector<dualVectorReal>(npe, "cy");
  data.nodes_coords_z_ = allocateVector<dualVectorReal>(npe, "cz");

  for (int lDof = 0; lDof < npe; ++lDof) {
    data.global_node_index_.view_host()(0, lDof) = lDof;

    int k = lDof / ((ORDER + 1) * (ORDER + 1));
    int j = (lDof / (ORDER + 1)) % (ORDER + 1);
    int i = lDof % (ORDER + 1);

    data.nodes_coords_x_.view_host()(lDof) = static_cast<float>(i) / ORDER;
    data.nodes_coords_y_.view_host()(lDof) = static_cast<float>(j) / ORDER;
    data.nodes_coords_z_.view_host()(lDof) = static_cast<float>(k) / ORDER;
  }

  data.global_node_index_.modify_host();
  data.nodes_coords_x_.modify_host();
  data.nodes_coords_y_.modify_host();
  data.nodes_coords_z_.modify_host();

  return model::ModelUnstruct<float, int>(data);
}

// =============================================================================
// DifferentiatorElasticElemUnstructTest  (IS_MODEL_ON_NODES = false)
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorElasticElemUnstructTest : public ::testing::Test {
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelUnstruct<float, int>;
  using Integral = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using Diff = DifferentiatorElastic<kOrder, Integral, Mesh, false>;

  void SetUp() override {
    ux_fwd = allocateVector<vectorReal>(kNumNodes, "ux_fwd");
    uy_fwd = allocateVector<vectorReal>(kNumNodes, "uy_fwd");
    uz_fwd = allocateVector<vectorReal>(kNumNodes, "uz_fwd");
    ux_adj = allocateVector<vectorReal>(kNumNodes, "ux_adj");
    uy_adj = allocateVector<vectorReal>(kNumNodes, "uy_adj");
    uz_adj = allocateVector<vectorReal>(kNumNodes, "uz_adj");
    ux_dt2 = allocateVector<vectorReal>(kNumNodes, "ux_dt2");
    uy_dt2 = allocateVector<vectorReal>(kNumNodes, "uy_dt2");
    uz_dt2 = allocateVector<vectorReal>(kNumNodes, "uz_dt2");
    gradRho = allocateVector<vectorReal>(kNumElements, "gradRho");
    gradLambda = allocateVector<vectorReal>(kNumElements, "gradLambda");
    gradMu = allocateVector<vectorReal>(kNumElements, "gradMu");

    Kokkos::deep_copy(ux_fwd, 0.0f);
    Kokkos::deep_copy(uy_fwd, 0.0f);
    Kokkos::deep_copy(uz_fwd, 0.0f);
    Kokkos::deep_copy(ux_adj, 0.0f);
    Kokkos::deep_copy(uy_adj, 0.0f);
    Kokkos::deep_copy(uz_adj, 0.0f);
    Kokkos::deep_copy(ux_dt2, 0.0f);
    Kokkos::deep_copy(uy_dt2, 0.0f);
    Kokkos::deep_copy(uz_dt2, 0.0f);
    setAt(gradRho, 0.0f, 0);
    setAt(gradLambda, 0.0f, 0);
    setAt(gradMu, 0.0f, 0);
  }

  vectorReal ux_fwd, uy_fwd, uz_fwd;
  vectorReal ux_adj, uy_adj, uz_adj;
  vectorReal ux_dt2, uy_dt2, uz_dt2;
  vectorReal gradRho, gradLambda, gradMu;
};

TYPED_TEST_SUITE(DifferentiatorElasticElemUnstructTest, OrderTypesU);

// --- Static constants ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest, OrderConstant) {
  EXPECT_EQ(TestFixture::Diff::kOrder, TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, IsModelOnNodesConstant) {
  EXPECT_FALSE(TestFixture::Diff::kIsModelOnNodes);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, PointsPerElementConstant) {
  EXPECT_EQ(TestFixture::Diff::kPointsPerElement, TestFixture::kNumNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest, GetOrderReturnsCorrectOrder) {
  typename TestFixture::Diff diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, IsModelOnNodesReturnsFalse) {
  typename TestFixture::Diff diff;
  EXPECT_FALSE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest, PrintDoesNotThrow) {
  typename TestFixture::Diff diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest, ZeroWavefieldsYieldZeroGradients) {
  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(getAt(this->gradRho, 0), 0.0f);
  EXPECT_FLOAT_EQ(getAt(this->gradLambda, 0), 0.0f);
  EXPECT_FLOAT_EQ(getAt(this->gradMu, 0), 0.0f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, UniformFieldGradRhoEqualsVolume) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradRho, 0), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, ConstantFieldGradLambdaIsZero) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_adj, 1.0f);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradLambda, 0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, ConstantFieldGradMuIsZero) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->uy_fwd, 1.0f);
  Kokkos::deep_copy(this->uz_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_adj, 1.0f);
  Kokkos::deep_copy(this->uy_adj, 1.0f);
  Kokkos::deep_copy(this->uz_adj, 1.0f);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradMu, 0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, GradLambdaNonZeroForLinearField) {
  // u_x = x => div(u) = 1 => gradLambda = ∫ 1 dΩ = 1.0
  constexpr int dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof) {
    int i = lDof % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    h_ux_fwd(lDof) = x;
    h_ux_adj(lDof) = x;
  }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradLambda, 0), 1.0f, 0.2f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, GradMuNonZeroForLinearField) {
  // u_x = x => ε_xx = 1 => 2·ε†:ε = 2 => gradMu = 2.0
  constexpr int dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof) {
    int i = lDof % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    h_ux_fwd(lDof) = x;
    h_ux_adj(lDof) = x;
  }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradMu, 0), 2.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, PureDilatationLambdaMuConsistency) {
  // u_x=x, u_y=y, u_z=z => div=3, ε_xx=ε_yy=ε_zz=1
  // gradLambda = 9.0, gradMu = 6.0
  constexpr int dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_uy_fwd = Kokkos::create_mirror_view(this->uy_fwd);
  auto h_uz_fwd = Kokkos::create_mirror_view(this->uz_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  auto h_uy_adj = Kokkos::create_mirror_view(this->uy_adj);
  auto h_uz_adj = Kokkos::create_mirror_view(this->uz_adj);
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof) {
    int k = lDof / (dim * dim);
    int j = (lDof / dim) % dim;
    int i = lDof % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    float y = static_cast<float>(j) / TestFixture::kOrder;
    float z = static_cast<float>(k) / TestFixture::kOrder;
    h_ux_fwd(lDof) = x;
    h_uy_fwd(lDof) = y;
    h_uz_fwd(lDof) = z;
    h_ux_adj(lDof) = x;
    h_uy_adj(lDof) = y;
    h_uz_adj(lDof) = z;
  }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->uy_fwd, h_uy_fwd);
  Kokkos::deep_copy(this->uz_fwd, h_uz_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);
  Kokkos::deep_copy(this->uy_adj, h_uy_adj);
  Kokkos::deep_copy(this->uz_adj, h_uz_adj);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradLambda, 0), 9.0f, 1.0f);
  EXPECT_NEAR(getAt(this->gradMu, 0), 6.0f, 1.0f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, LinearShearFieldProducesCorrectGradMu) {
  // Pure shear: u_x=y, u_y=x => div=0 => gradLambda=0
  // ε_xy = 1 => 2·ε:ε = 4·(ε_xy²) = 4 => gradMu = 4.0
  constexpr int dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_uy_fwd = Kokkos::create_mirror_view(this->uy_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  auto h_uy_adj = Kokkos::create_mirror_view(this->uy_adj);
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof) {
    int i = lDof % dim;
    int j = (lDof / dim) % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    float y = static_cast<float>(j) / TestFixture::kOrder;
    h_ux_fwd(lDof) = y;
    h_uy_fwd(lDof) = x;
    h_ux_adj(lDof) = y;
    h_uy_adj(lDof) = x;
  }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->uy_fwd, h_uy_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);
  Kokkos::deep_copy(this->uy_adj, h_uy_adj);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradLambda, 0), 0.0f, 1e-5f);
  EXPECT_NEAR(getAt(this->gradMu, 0), 4.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, ComputeAccumulatesIntoExistingGradient) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);
  setAt(this->gradRho, 5.0f, 0);

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradRho, 0), 6.0f, 1e-5f);
}

// Fills prev/prevprev levels and the matching dt2 = (prevprev - 2 prev + n) / dt^2 for every component.
static void fillTimeLevels(vectorReal const (&fwd)[3], vectorReal const (&n)[3], vectorReal const (&dt2)[3], int nNodes,
                           float dt, vectorReal (&prev)[3], vectorReal (&prevprev)[3]) {
  for (int c = 0; c < 3; ++c) {
    prev[c] = allocateVector<vectorReal>(nNodes, "prev");
    prevprev[c] = allocateVector<vectorReal>(nNodes, "prevprev");
    auto h_fwd = Kokkos::create_mirror_view(fwd[c]);
    auto h_n = Kokkos::create_mirror_view(n[c]);
    auto h_prev = Kokkos::create_mirror_view(prev[c]);
    auto h_prevprev = Kokkos::create_mirror_view(prevprev[c]);
    auto h_dt2 = Kokkos::create_mirror_view(dt2[c]);
    for (int i = 0; i < nNodes; ++i) {
      h_fwd(i) = 0.25f * (c + 1) + 0.125f * i;
      h_n(i) = 0.5f * (c + 1) - 0.0625f * i;
      h_prev(i) = 0.375f * (c + 2) + 0.03125f * i;
      h_prevprev(i) = 0.125f * (c + 3) - 0.015625f * i;
      h_dt2(i) = (h_prevprev(i) - 2.0f * h_prev(i) + h_n(i)) / (dt * dt);
    }
    Kokkos::deep_copy(fwd[c], h_fwd);
    Kokkos::deep_copy(n[c], h_n);
    Kokkos::deep_copy(prev[c], h_prev);
    Kokkos::deep_copy(prevprev[c], h_prevprev);
    Kokkos::deep_copy(dt2[c], h_dt2);
  }
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, TimeLevelsMatchPrecomputedDt2) {
  float const dt = 0.5f;
  vectorReal prev[3], prevprev[3];
  fillTimeLevels({this->ux_fwd, this->uy_fwd, this->uz_fwd}, {this->ux_adj, this->uy_adj, this->uz_adj},
                 {this->ux_dt2, this->uy_dt2, this->uz_dt2}, TestFixture::kNumNodes, dt, prev, prevprev);
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();
  typename TestFixture::Diff diff;
  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);

  WavefieldViewBackwardElastic bwdDt2(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2,
                                      this->uz_dt2);
  GradientDataElastic dataDt2(fwd, bwdDt2, GradientElastic(this->gradRho, this->gradLambda, this->gradMu));
  diff.compute(mesh, dataDt2, dt);
  float const rho = getAt(this->gradRho, 0), lambda = getAt(this->gradLambda, 0), mu = getAt(this->gradMu, 0);

  Kokkos::deep_copy(this->gradRho, 0.0f);
  Kokkos::deep_copy(this->gradLambda, 0.0f);
  Kokkos::deep_copy(this->gradMu, 0.0f);
  WavefieldViewBackwardElastic bwdLevels(this->ux_adj, this->uy_adj, this->uz_adj, prev[0], prev[1], prev[2],
                                         prevprev[0], prevprev[1], prevprev[2]);
  EXPECT_TRUE(bwdLevels.fromTimeLevels());
  EXPECT_FALSE(bwdDt2.fromTimeLevels());
  GradientDataElastic dataLevels(fwd, bwdLevels, GradientElastic(this->gradRho, this->gradLambda, this->gradMu));
  diff.compute(mesh, dataLevels, dt);

  EXPECT_NE(rho, 0.0f);
  EXPECT_NEAR(getAt(this->gradRho, 0), rho, 1e-5f * std::fabs(rho));
  EXPECT_FLOAT_EQ(getAt(this->gradLambda, 0), lambda);
  EXPECT_FLOAT_EQ(getAt(this->gradMu, 0), mu);
}

// --- Node-based unstructured ---

template <typename OrderWrapper>
class DifferentiatorElasticNodeUnstructTest : public ::testing::Test {
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);

  using Mesh = model::ModelUnstruct<float, int>;
  using Integral = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using DiffNode = DifferentiatorElastic<kOrder, Integral, Mesh, true>;
  using DiffElem = DifferentiatorElastic<kOrder, Integral, Mesh, false>;

  void SetUp() override {
    ux_fwd = allocateVector<vectorReal>(kNumNodes, "ux_fwd");
    uy_fwd = allocateVector<vectorReal>(kNumNodes, "uy_fwd");
    uz_fwd = allocateVector<vectorReal>(kNumNodes, "uz_fwd");
    ux_adj = allocateVector<vectorReal>(kNumNodes, "ux_adj");
    uy_adj = allocateVector<vectorReal>(kNumNodes, "uy_adj");
    uz_adj = allocateVector<vectorReal>(kNumNodes, "uz_adj");
    ux_dt2 = allocateVector<vectorReal>(kNumNodes, "ux_dt2");
    uy_dt2 = allocateVector<vectorReal>(kNumNodes, "uy_dt2");
    uz_dt2 = allocateVector<vectorReal>(kNumNodes, "uz_dt2");
    gradRho = allocateVector<vectorReal>(kNumNodes, "gradRho");
    gradLambda = allocateVector<vectorReal>(kNumNodes, "gradLambda");
    gradMu = allocateVector<vectorReal>(kNumNodes, "gradMu");

    Kokkos::deep_copy(ux_fwd, 0.0f);
    Kokkos::deep_copy(uy_fwd, 0.0f);
    Kokkos::deep_copy(uz_fwd, 0.0f);
    Kokkos::deep_copy(ux_adj, 0.0f);
    Kokkos::deep_copy(uy_adj, 0.0f);
    Kokkos::deep_copy(uz_adj, 0.0f);
    Kokkos::deep_copy(ux_dt2, 0.0f);
    Kokkos::deep_copy(uy_dt2, 0.0f);
    Kokkos::deep_copy(uz_dt2, 0.0f);
    Kokkos::deep_copy(gradRho, 0.0f);
    Kokkos::deep_copy(gradLambda, 0.0f);
    Kokkos::deep_copy(gradMu, 0.0f);
  }

  float sumGrad(vectorReal const& v) const {
    float s = 0.0f;
    auto h_v = toHost(v);
    for (int i = 0; i < kNumNodes; ++i) s += h_v(i);
    return s;
  }

  float sumGradRho() const { return sumGrad(gradRho); }
  float sumGradLambda() const { return sumGrad(gradLambda); }
  float sumGradMu() const { return sumGrad(gradMu); }

  vectorReal ux_fwd, uy_fwd, uz_fwd;
  vectorReal ux_adj, uy_adj, uz_adj;
  vectorReal ux_dt2, uy_dt2, uz_dt2;
  vectorReal gradRho, gradLambda, gradMu;
};

TYPED_TEST_SUITE(DifferentiatorElasticNodeUnstructTest, OrderTypesU);

TYPED_TEST(DifferentiatorElasticNodeUnstructTest, ZeroWavefieldsYieldZeroGradients) {
  typename TestFixture::DiffNode diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->sumGradRho(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradLambda(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradMu(), 0.0f);
}

TYPED_TEST(DifferentiatorElasticNodeUnstructTest, UniformFieldGradRhoSumsToVolume) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  typename TestFixture::DiffNode diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradRho(), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticNodeUnstructTest, NodeBasedSumEqualsElementBasedResult) {
  constexpr int dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_uy_fwd = Kokkos::create_mirror_view(this->uy_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  auto h_uy_adj = Kokkos::create_mirror_view(this->uy_adj);
  auto h_ux_dt2 = Kokkos::create_mirror_view(this->ux_dt2);
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof) {
    int i = lDof % dim;
    int j = (lDof / dim) % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    float y = static_cast<float>(j) / TestFixture::kOrder;
    h_ux_fwd(lDof) = x;
    h_uy_fwd(lDof) = y;
    h_ux_adj(lDof) = x;
    h_uy_adj(lDof) = y;
    h_ux_dt2(lDof) = 1.0f;
  }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->uy_fwd, h_uy_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);
  Kokkos::deep_copy(this->uy_adj, h_uy_adj);
  Kokkos::deep_copy(this->ux_dt2, h_ux_dt2);

  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  // Node-based
  typename TestFixture::DiffNode diffNode;
  {
    WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
    WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2,
                                     this->uz_dt2);
    GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
    GradientDataElastic data(fwd, bwd, grad);
    diffNode.compute(mesh, data, 0.001f);
  }
  float nodeRho = this->sumGradRho();
  float nodeLambda = this->sumGradLambda();
  float nodeMu = this->sumGradMu();

  // Element-based
  auto gradRhoElem = allocateVector<vectorReal>(1, "gradRhoElem");
  auto gradLambdaElem = allocateVector<vectorReal>(1, "gradLambdaElem");
  auto gradMuElem = allocateVector<vectorReal>(1, "gradMuElem");
  setAt(gradRhoElem, 0.0f, 0);
  setAt(gradLambdaElem, 0.0f, 0);
  setAt(gradMuElem, 0.0f, 0);

  typename TestFixture::DiffElem diffElem;
  {
    WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
    WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2,
                                     this->uz_dt2);
    GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
    GradientDataElastic data(fwd, bwd, grad);
    diffElem.compute(mesh, data, 0.001f);
  }

  EXPECT_NEAR(nodeRho, getAt(gradRhoElem, 0), 1e-4f);
  EXPECT_NEAR(nodeLambda, getAt(gradLambdaElem, 0), 1e-4f);
  EXPECT_NEAR(nodeMu, getAt(gradMuElem, 0), 1e-4f);
}

// --- Polymorphic interface ---

TYPED_TEST(DifferentiatorElasticNodeUnstructTest, PolymorphicInterface) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  std::unique_ptr<Differentiator> diff = std::make_unique<typename TestFixture::DiffNode>();
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  EXPECT_EQ(diff->getOrder(), TestFixture::kOrder);
  EXPECT_TRUE(diff->isModelOnNodes());

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  EXPECT_NO_THROW(diff->compute(mesh, data, 0.001f));
  EXPECT_GT(this->sumGradRho(), 0.0f);
}

TYPED_TEST(DifferentiatorElasticNodeUnstructTest, TimeLevelsMatchPrecomputedDt2) {
  float const dt = 0.5f;
  vectorReal prev[3], prevprev[3];
  fillTimeLevels({this->ux_fwd, this->uy_fwd, this->uz_fwd}, {this->ux_adj, this->uy_adj, this->uz_adj},
                 {this->ux_dt2, this->uy_dt2, this->uz_dt2}, TestFixture::kNumNodes, dt, prev, prevprev);
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();
  typename TestFixture::DiffNode diff;
  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);

  WavefieldViewBackwardElastic bwdDt2(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2,
                                      this->uz_dt2);
  GradientDataElastic dataDt2(fwd, bwdDt2, GradientElastic(this->gradRho, this->gradLambda, this->gradMu));
  diff.compute(mesh, dataDt2, dt);
  auto h_grad_rho = toHost(this->gradRho);
  std::vector<float> rho(TestFixture::kNumNodes);
  for (int i = 0; i < TestFixture::kNumNodes; ++i) rho[i] = h_grad_rho(i);
  Kokkos::deep_copy(this->gradRho, 0.0f);
  float const lambda = this->sumGradLambda(), mu = this->sumGradMu();
  Kokkos::deep_copy(this->gradLambda, 0.0f);
  Kokkos::deep_copy(this->gradMu, 0.0f);

  WavefieldViewBackwardElastic bwdLevels(this->ux_adj, this->uy_adj, this->uz_adj, prev[0], prev[1], prev[2],
                                         prevprev[0], prevprev[1], prevprev[2]);
  GradientDataElastic dataLevels(fwd, bwdLevels, GradientElastic(this->gradRho, this->gradLambda, this->gradMu));
  diff.compute(mesh, dataLevels, dt);

  h_grad_rho = toHost(this->gradRho);
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
    EXPECT_NEAR(h_grad_rho(i), rho[i], 1e-5f * std::fabs(rho[i]) + 1e-7f) << "node " << i;
  EXPECT_FLOAT_EQ(this->sumGradLambda(), lambda);
  EXPECT_FLOAT_EQ(this->sumGradMu(), mu);
}

}  // namespace test
}  // namespace gradient

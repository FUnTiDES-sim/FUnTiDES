#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "Integrals.h"
#include "data_type.h"
#include "device_view_access.h"
#include "differentiator_data_elastic.h"
#include "differentiator_elastic.h"
#include "model_struct.h"

namespace gradient {
namespace test {

using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

// =============================================================================
// Order wrappers
// =============================================================================

struct Order1 {
  static constexpr int kOrder = 1;
};
struct Order2 {
  static constexpr int kOrder = 2;
};
struct Order3 {
  static constexpr int kOrder = 3;
};

// Création dynamique de la liste de tests selon l'ordre max
#if MAX_DIFFERENTIATOR_ACOUSTIC_ORDER == 1
using OrderTypes = ::testing::Types<Order1>;
#elif MAX_DIFFERENTIATOR_ACOUSTIC_ORDER == 2
using OrderTypes = ::testing::Types<Order1, Order2>;
#elif MAX_DIFFERENTIATOR_ACOUSTIC_ORDER == 3
using OrderTypes = ::testing::Types<Order1, Order2, Order3>;
#endif

// =============================================================================
// Helper: 1×1×1 structured mesh (1 element, unit cube [0,1]³)
// =============================================================================

template <int ORDER>
static model::ModelStruct<float, int, ORDER> makeMesh1x1x1() {
  model::ModelStructData<float, int> data;
  data.ex_ = 1;
  data.ey_ = 1;
  data.ez_ = 1;
  data.dx_ = 1.0f;
  data.dy_ = 1.0f;
  data.dz_ = 1.0f;
  data.ox_ = 0.0f;
  data.oy_ = 0.0f;
  data.oz_ = 0.0f;
  data.isModelOnNodes_ = false;
  data.isElastic_ = true;
  return model::ModelStruct<float, int, ORDER>(data);
}

// =============================================================================
// DifferentiatorElasticElemTest  (IS_MODEL_ON_NODES = false)
//
// Mesh: 1 element, (ORDER+1)^3 nodes.
// Gradients are element-indexed (size = 1).
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorElasticElemTest : public ::testing::Test {
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelStruct<float, int, kOrder>;
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

  void makeDataAndCompute(Diff& diff, float dt = 0.001f) {
    auto mesh = makeMesh1x1x1<kOrder>();
    WavefieldViewForwardElastic fwd(ux_fwd, uy_fwd, uz_fwd);
    WavefieldViewBackwardElastic bwd(ux_adj, uy_adj, uz_adj, ux_dt2, uy_dt2, uz_dt2);
    GradientElastic grad(gradRho, gradLambda, gradMu);
    GradientDataElastic data(fwd, bwd, grad);
    diff.initGeometricMassMatrix(mesh);
    diff.compute(mesh, data, dt);
  }

  vectorReal ux_fwd, uy_fwd, uz_fwd;
  vectorReal ux_adj, uy_adj, uz_adj;
  vectorReal ux_dt2, uy_dt2, uz_dt2;
  vectorReal gradRho, gradLambda, gradMu;
};

TYPED_TEST_SUITE(DifferentiatorElasticElemTest, OrderTypes);

// --- Static constants ---

TYPED_TEST(DifferentiatorElasticElemTest, OrderConstant) { EXPECT_EQ(TestFixture::Diff::kOrder, TestFixture::kOrder); }

TYPED_TEST(DifferentiatorElasticElemTest, IsModelOnNodesConstant) { EXPECT_FALSE(TestFixture::Diff::kIsModelOnNodes); }

TYPED_TEST(DifferentiatorElasticElemTest, PointsPerElementConstant) {
  EXPECT_EQ(TestFixture::Diff::kPointsPerElement, TestFixture::kNumNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorElasticElemTest, GetOrderReturnsCorrectOrder) {
  typename TestFixture::Diff diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticElemTest, IsModelOnNodesReturnsFalse) {
  typename TestFixture::Diff diff;
  EXPECT_FALSE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorElasticElemTest, PrintDoesNotThrow) {
  typename TestFixture::Diff diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorElasticElemTest, ZeroWavefieldsYieldZeroGradients) {
  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_FLOAT_EQ(getAt(this->gradRho, 0), 0.0f);
  EXPECT_FLOAT_EQ(getAt(this->gradLambda, 0), 0.0f);
  EXPECT_FLOAT_EQ(getAt(this->gradMu, 0), 0.0f);
}

TYPED_TEST(DifferentiatorElasticElemTest, UniformFieldGradRhoEqualsVolume) {
  // ux_fwd = 1, ux_dt2 = 1 on all nodes =>
  //   gradRho = ∫_Ω (ü†_x · u_x) dΩ = ∫_Ω 1 dΩ = 1.0 (unit cube)
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(getAt(this->gradRho, 0), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradRhoAllThreeComponents) {
  // All forward and dt2 components active:
  //   gradRho = ∫(ü†_x·u_x + ü†_y·u_y + ü†_z·u_z) dΩ
  // u_x=1,u_y=2,u_z=3 with all dt2=1  =>  gradRho = (1+2+3) * volume = 6.0
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->uy_fwd, 2.0f);
  Kokkos::deep_copy(this->uz_fwd, 3.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);
  Kokkos::deep_copy(this->uy_dt2, 1.0f);
  Kokkos::deep_copy(this->uz_dt2, 1.0f);

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(getAt(this->gradRho, 0), 6.0f, 1e-4f);
}

TYPED_TEST(DifferentiatorElasticElemTest, ConstantFieldGradLambdaIsZero) {
  // Constant displacement => ∇u = 0 => div(u†) = div(u) = 0 => gradLambda = 0
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->uy_fwd, 1.0f);
  Kokkos::deep_copy(this->uz_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_adj, 1.0f);
  Kokkos::deep_copy(this->uy_adj, 1.0f);
  Kokkos::deep_copy(this->uz_adj, 1.0f);

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(getAt(this->gradLambda, 0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, ConstantFieldGradMuIsZero) {
  // Constant displacement => ε(u) = 0 => gradMu = 0
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->uy_fwd, 1.0f);
  Kokkos::deep_copy(this->uz_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_adj, 1.0f);
  Kokkos::deep_copy(this->uy_adj, 1.0f);
  Kokkos::deep_copy(this->uz_adj, 1.0f);

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(getAt(this->gradMu, 0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradRhoScalesWithAmplitude) {
  // Doubling ux_fwd doubles gradRho (linearity)
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);
  float single = getAt(this->gradRho, 0);

  // Reset and double
  setAt(this->gradRho, 0.0f, 0);
  Kokkos::deep_copy(this->ux_fwd, 2.0f);

  this->makeDataAndCompute(diff);

  EXPECT_NEAR(getAt(this->gradRho, 0), 2.0f * single, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, ComputeAccumulatesIntoExistingGradient) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);
  setAt(this->gradRho, 5.0f, 0);

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  // 5.0 (initial) + 1.0 (volume of unit cube) = 6.0
  EXPECT_NEAR(getAt(this->gradRho, 0), 6.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradLambdaNonZeroForLinearField) {
  // Linear displacement u_x = x => div(u) = ∂u_x/∂x = 1
  // If both forward and adjoint have u_x = x:
  //   gradLambda = ∫ div(u†)·div(u) dΩ = ∫ 1·1 dΩ = 1.0
  //
  // On structure mesh [0,1]³ with nodes at x = i/ORDER:
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i) {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        h_ux_fwd(gIdx) = x;
        h_ux_adj(gIdx) = x;
      }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  // ∫₀¹ ∫₀¹ ∫₀¹ 1·1 dx dy dz = 1.0
  // Tolerance is wider for Order 3 where uniform node spacing diverges from
  // GLL quadrature points, introducing discretization error.
  EXPECT_NEAR(getAt(this->gradLambda, 0), 1.0f, 0.2f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradMuNonZeroForLinearField) {
  // Linear displacement u_x = x => ε_xx = 1, all other ε_ij = 0
  // If both forward and adjoint have u_x = x:
  //   2·ε†:ε = 2·(1·1) = 2
  //   gradMu = ∫ 2 dΩ = 2.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i) {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        h_ux_fwd(gIdx) = x;
        h_ux_adj(gIdx) = x;
      }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  // 2·ε†:ε = 2·ε_xx²  = 2·1 = 2; ∫_Ω 2 dΩ = 2.0
  EXPECT_NEAR(getAt(this->gradMu, 0), 2.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemTest, LinearShearFieldProducesCorrectGradMu) {
  // Pure shear: u_x = y, u_y = x  => ε_xy = 1, all diagonal ε = 0,
  //   div(u) = 0
  // 2·ε†:ε = 4·(ε_xy†·ε_xy) = 4·1 = 4
  // gradMu = ∫ 4 dΩ = 4.0
  // gradLambda = ∫ 0·0 dΩ = 0.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_uy_fwd = Kokkos::create_mirror_view(this->uy_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  auto h_uy_adj = Kokkos::create_mirror_view(this->uy_adj);
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i) {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        float y = static_cast<float>(j) / TestFixture::kOrder;
        h_ux_fwd(gIdx) = y;
        h_uy_fwd(gIdx) = x;
        h_ux_adj(gIdx) = y;
        h_uy_adj(gIdx) = x;
      }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->uy_fwd, h_uy_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);
  Kokkos::deep_copy(this->uy_adj, h_uy_adj);

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradLambda, 0), 0.0f, 1e-5f);
  EXPECT_NEAR(getAt(this->gradMu, 0), 4.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemTest, PureDilatationLambdaMuConsistency) {
  // Pure dilatation: u_x = x, u_y = y, u_z = z
  //   div(u) = 3, ε_xx = ε_yy = ε_zz = 1, off-diag = 0
  //   gradLambda = ∫ 3·3 dΩ = 9.0
  //   2·ε†:ε = 2·(1+1+1) = 6
  //   gradMu = ∫ 6 dΩ = 6.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_uy_fwd = Kokkos::create_mirror_view(this->uy_fwd);
  auto h_uz_fwd = Kokkos::create_mirror_view(this->uz_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  auto h_uy_adj = Kokkos::create_mirror_view(this->uy_adj);
  auto h_uz_adj = Kokkos::create_mirror_view(this->uz_adj);
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i) {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        float y = static_cast<float>(j) / TestFixture::kOrder;
        float z = static_cast<float>(k) / TestFixture::kOrder;
        h_ux_fwd(gIdx) = x;
        h_uy_fwd(gIdx) = y;
        h_uz_fwd(gIdx) = z;
        h_ux_adj(gIdx) = x;
        h_uy_adj(gIdx) = y;
        h_uz_adj(gIdx) = z;
      }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->uy_fwd, h_uy_fwd);
  Kokkos::deep_copy(this->uz_fwd, h_uz_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);
  Kokkos::deep_copy(this->uy_adj, h_uy_adj);
  Kokkos::deep_copy(this->uz_adj, h_uz_adj);

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(getAt(this->gradLambda, 0), 9.0f, 1.0f);
  EXPECT_NEAR(getAt(this->gradMu, 0), 6.0f, 1.0f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradRhoIndependentOfStrainKernels) {
  // gradRho depends only on dt2 and fwd, not on adjoint displacement gradients
  // Set ux_fwd=1, ux_dt2=1 but also set non-trivial adjoint displacement
  // => gradRho should still be 1.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_ux_dt2 = Kokkos::create_mirror_view(this->ux_dt2);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  auto h_uy_adj = Kokkos::create_mirror_view(this->uy_adj);
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i) {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        h_ux_fwd(gIdx) = 1.0f;
        h_ux_dt2(gIdx) = 1.0f;
        // Non-trivial adjoint displacement for strain kernels
        h_ux_adj(gIdx) = x;
        h_uy_adj(gIdx) = x;
      }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->ux_dt2, h_ux_dt2);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);
  Kokkos::deep_copy(this->uy_adj, h_uy_adj);

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  // gradRho only depends on dt2·fwd dot product: 1 * 1 * volume = 1.0
  EXPECT_NEAR(getAt(this->gradRho, 0), 1.0f, 1e-5f);
}

// --- Polymorphic interface ---

TYPED_TEST(DifferentiatorElasticElemTest, PolymorphicInterface) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  std::unique_ptr<Differentiator> diff = std::make_unique<typename TestFixture::Diff>();
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  EXPECT_EQ(diff->getOrder(), TestFixture::kOrder);
  EXPECT_FALSE(diff->isModelOnNodes());

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  EXPECT_NO_THROW(diff->compute(mesh, data, 0.001f));
  EXPECT_GT(getAt(this->gradRho, 0), 0.0f);
}

// =============================================================================
// DifferentiatorElasticNodeTest  (IS_MODEL_ON_NODES = true)
//
// Mesh: 1 element, (ORDER+1)^3 nodes.
// Gradients are node-indexed; contributions are scattered with ATOMICADD.
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorElasticNodeTest : public ::testing::Test {
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelStruct<float, int, kOrder>;
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

TYPED_TEST_SUITE(DifferentiatorElasticNodeTest, OrderTypes);

// --- Static constants ---

TYPED_TEST(DifferentiatorElasticNodeTest, IsModelOnNodesConstant) {
  EXPECT_TRUE(TestFixture::DiffNode::kIsModelOnNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorElasticNodeTest, GetOrderReturnsCorrectOrder) {
  typename TestFixture::DiffNode diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticNodeTest, IsModelOnNodesReturnsTrue) {
  typename TestFixture::DiffNode diff;
  EXPECT_TRUE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorElasticNodeTest, PrintDoesNotThrow) {
  typename TestFixture::DiffNode diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorElasticNodeTest, ZeroWavefieldsYieldZeroGradients) {
  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->sumGradRho(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradLambda(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradMu(), 0.0f);
}

TYPED_TEST(DifferentiatorElasticNodeTest, UniformFieldGradRhoSumsToVolume) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradRho(), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticNodeTest, ConstantFieldGradLambdaSumsToZero) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_adj, 1.0f);

  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradLambda(), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticNodeTest, NodeBasedSumEqualsElementBasedResult) {
  // For single-element mesh, ∑ node_gradients == elem_gradient
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  auto h_ux_fwd = Kokkos::create_mirror_view(this->ux_fwd);
  auto h_uy_fwd = Kokkos::create_mirror_view(this->uy_fwd);
  auto h_ux_adj = Kokkos::create_mirror_view(this->ux_adj);
  auto h_uy_adj = Kokkos::create_mirror_view(this->uy_adj);
  auto h_ux_dt2 = Kokkos::create_mirror_view(this->ux_dt2);
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i) {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        float y = static_cast<float>(j) / TestFixture::kOrder;
        h_ux_fwd(gIdx) = x;
        h_uy_fwd(gIdx) = y;
        h_ux_adj(gIdx) = x;
        h_uy_adj(gIdx) = y;
        h_ux_dt2(gIdx) = 1.0f;
      }
  Kokkos::deep_copy(this->ux_fwd, h_ux_fwd);
  Kokkos::deep_copy(this->uy_fwd, h_uy_fwd);
  Kokkos::deep_copy(this->ux_adj, h_ux_adj);
  Kokkos::deep_copy(this->uy_adj, h_uy_adj);
  Kokkos::deep_copy(this->ux_dt2, h_ux_dt2);

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

TYPED_TEST(DifferentiatorElasticNodeTest, PolymorphicInterface) {
  Kokkos::deep_copy(this->ux_fwd, 1.0f);
  Kokkos::deep_copy(this->ux_dt2, 1.0f);

  std::unique_ptr<Differentiator> diff = std::make_unique<typename TestFixture::DiffNode>();
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  EXPECT_EQ(diff->getOrder(), TestFixture::kOrder);
  EXPECT_TRUE(diff->isModelOnNodes());

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  EXPECT_NO_THROW(diff->compute(mesh, data, 0.001f));
  EXPECT_GT(this->sumGradRho(), 0.0f);
}

// =============================================================================
// GeometricMassMatrix tests — DifferentiatorElasticElemTest
// (IS_MODEL_ON_NODES = false)
// =============================================================================

TYPED_TEST(DifferentiatorElasticElemTest, GeometricMassMatrixSizeAfterCompute) {
  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_EQ(diff.getGeometricMassMatrix().extent(0), static_cast<size_t>(TestFixture::kNumNodes));
}

TYPED_TEST(DifferentiatorElasticElemTest, GeometricMassMatrixAllPositive) {
  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  auto gmm = toHost(diff.getGeometricMassMatrix());
  for (int i = 0; i < TestFixture::kNumNodes; ++i) EXPECT_GT(gmm(i), 0.0f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GeometricMassMatrixSumsToVolume) {
  // For a unit cube [0,1]^3, sum of all nodal volumes = 1.0
  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  auto gmm = toHost(diff.getGeometricMassMatrix());
  float total = 0.0f;
  for (int i = 0; i < TestFixture::kNumNodes; ++i) total += gmm(i);
  EXPECT_NEAR(total, 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GeometricMassMatrixIndependentOfWavefields) {
  // GMM depends only on geometry, not wavefields
  typename TestFixture::Diff diff;

  // First compute: zero wavefields
  this->makeDataAndCompute(diff);
  auto& gmm0 = diff.getGeometricMassMatrix();
  vectorReal saved = allocateVector<vectorReal>(TestFixture::kNumNodes, "saved_gmm");
  Kokkos::deep_copy(saved, gmm0);

  // Second compute: non-zero wavefields
  Kokkos::deep_copy(this->ux_fwd, 3.0f);
  Kokkos::deep_copy(this->ux_dt2, 2.5f);
  setAt(this->gradRho, 0.0f, 0);
  setAt(this->gradLambda, 0.0f, 0);
  setAt(this->gradMu, 0.0f, 0);
  this->makeDataAndCompute(diff);

  auto gmm1 = toHost(diff.getGeometricMassMatrix());
  auto h_saved = toHost(saved);
  for (int i = 0; i < TestFixture::kNumNodes; ++i) EXPECT_FLOAT_EQ(gmm1(i), h_saved(i));
}

TYPED_TEST(DifferentiatorElasticElemTest, GeometricMassMatrixPolymorphicAccess) {
  std::unique_ptr<Differentiator> diff = std::make_unique<typename TestFixture::Diff>();
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff->initGeometricMassMatrix(mesh);
  diff->compute(mesh, data, 0.001f);

  auto gmm = toHost(diff->getGeometricMassMatrix());
  EXPECT_EQ(gmm.extent(0), static_cast<size_t>(TestFixture::kNumNodes));
  float total = 0.0f;
  for (int i = 0; i < TestFixture::kNumNodes; ++i) total += gmm(i);
  EXPECT_NEAR(total, 1.0f, 1e-5f);
}

// =============================================================================
// GeometricMassMatrix tests — DifferentiatorElasticNodeTest
// (IS_MODEL_ON_NODES = true)
// =============================================================================

TYPED_TEST(DifferentiatorElasticNodeTest, GeometricMassMatrixSizeAfterCompute) {
  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.initGeometricMassMatrix(mesh);
  diff.compute(mesh, data, 0.001f);

  EXPECT_EQ(diff.getGeometricMassMatrix().extent(0), static_cast<size_t>(TestFixture::kNumNodes));
}

TYPED_TEST(DifferentiatorElasticNodeTest, GeometricMassMatrixSumsToVolume) {
  // For a unit cube [0,1]^3, sum of all nodal volumes = 1.0
  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.initGeometricMassMatrix(mesh);
  diff.compute(mesh, data, 0.001f);

  auto gmm = toHost(diff.getGeometricMassMatrix());
  float total = 0.0f;
  for (int i = 0; i < TestFixture::kNumNodes; ++i) total += gmm(i);
  EXPECT_NEAR(total, 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticNodeTest, GeometricMassMatrixSameForNodeAndElemVariant) {
  // The geometric mass matrix must agree between element-based and node-based
  // differentiators on the same mesh (geometry-only quantity).
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  typename TestFixture::DiffNode diffNode;
  {
    WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
    WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj, this->ux_dt2, this->uy_dt2,
                                     this->uz_dt2);
    GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
    GradientDataElastic data(fwd, bwd, grad);
    diffNode.initGeometricMassMatrix(mesh);
    diffNode.compute(mesh, data, 0.001f);
  }

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
    diffElem.initGeometricMassMatrix(mesh);
    diffElem.compute(mesh, data, 0.001f);
  }

  auto gmmNode = toHost(diffNode.getGeometricMassMatrix());
  auto gmmElem = toHost(diffElem.getGeometricMassMatrix());
  for (int i = 0; i < TestFixture::kNumNodes; ++i) EXPECT_FLOAT_EQ(gmmNode(i), gmmElem(i));
}

}  // namespace test
}  // namespace gradient

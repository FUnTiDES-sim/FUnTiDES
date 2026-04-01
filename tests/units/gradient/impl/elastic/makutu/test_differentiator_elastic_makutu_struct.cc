#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "data_type.h"
#include "differentiator_elastic.h"
#include "differentiator_data_elastic.h"
#include "fe/Integrals.hpp"
#include "model_struct.h"

namespace gradient
{
namespace test
{

// =============================================================================
// Order wrappers
// =============================================================================

struct Order1
{
  static constexpr int kOrder = 1;
};
struct Order2
{
  static constexpr int kOrder = 2;
};
struct Order3
{
  static constexpr int kOrder = 3;
};

using OrderTypes = ::testing::Types<Order1, Order2, Order3>;

// =============================================================================
// Helper: 1×1×1 structured mesh (1 element, unit cube [0,1]³)
// =============================================================================

template <int ORDER>
static model::ModelStruct<float, int, ORDER> makeMesh1x1x1()
{
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
class DifferentiatorElasticElemTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelStruct<float, int, kOrder>;
  using Integral =
      typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using Diff = DifferentiatorElastic<kOrder, Integral, Mesh, false>;

  void SetUp() override
  {
    ux_fwd = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "ux_fwd");
    uy_fwd = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uy_fwd");
    uz_fwd = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uz_fwd");
    ux_adj = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "ux_adj");
    uy_adj = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uy_adj");
    uz_adj = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uz_adj");
    ux_dt2 = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "ux_dt2");
    uy_dt2 = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uy_dt2");
    uz_dt2 = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uz_dt2");
    gradRho = allocateVector<VECTOR_REAL_VIEW>(kNumElements, "gradRho");
    gradLambda = allocateVector<VECTOR_REAL_VIEW>(kNumElements, "gradLambda");
    gradMu = allocateVector<VECTOR_REAL_VIEW>(kNumElements, "gradMu");

    for (int i = 0; i < kNumNodes; ++i)
    {
      ux_fwd(i) = 0.0f;
      uy_fwd(i) = 0.0f;
      uz_fwd(i) = 0.0f;
      ux_adj(i) = 0.0f;
      uy_adj(i) = 0.0f;
      uz_adj(i) = 0.0f;
      ux_dt2(i) = 0.0f;
      uy_dt2(i) = 0.0f;
      uz_dt2(i) = 0.0f;
    }
    gradRho(0) = 0.0f;
    gradLambda(0) = 0.0f;
    gradMu(0) = 0.0f;
  }

  void makeDataAndCompute(Diff& diff, float dt = 0.001f)
  {
    auto mesh = makeMesh1x1x1<kOrder>();
    WavefieldViewForwardElastic fwd(ux_fwd, uy_fwd, uz_fwd);
    WavefieldViewBackwardElastic bwd(ux_adj, uy_adj, uz_adj, ux_dt2, uy_dt2,
                                     uz_dt2);
    GradientElastic grad(gradRho, gradLambda, gradMu);
    GradientDataElastic data(fwd, bwd, grad);
    diff.compute(mesh, data, dt);
  }

  VECTOR_REAL_VIEW ux_fwd, uy_fwd, uz_fwd;
  VECTOR_REAL_VIEW ux_adj, uy_adj, uz_adj;
  VECTOR_REAL_VIEW ux_dt2, uy_dt2, uz_dt2;
  VECTOR_REAL_VIEW gradRho, gradLambda, gradMu;
};

TYPED_TEST_SUITE(DifferentiatorElasticElemTest, OrderTypes);

// --- Static constants ---

TYPED_TEST(DifferentiatorElasticElemTest, OrderConstant)
{
  EXPECT_EQ(TestFixture::Diff::kOrder, TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticElemTest, IsModelOnNodesConstant)
{
  EXPECT_FALSE(TestFixture::Diff::kIsModelOnNodes);
}

TYPED_TEST(DifferentiatorElasticElemTest, PointsPerElementConstant)
{
  EXPECT_EQ(TestFixture::Diff::kPointsPerElement, TestFixture::kNumNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorElasticElemTest, GetOrderReturnsCorrectOrder)
{
  typename TestFixture::Diff diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticElemTest, IsModelOnNodesReturnsFalse)
{
  typename TestFixture::Diff diff;
  EXPECT_FALSE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorElasticElemTest, PrintDoesNotThrow)
{
  typename TestFixture::Diff diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorElasticElemTest, ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_FLOAT_EQ(this->gradRho(0), 0.0f);
  EXPECT_FLOAT_EQ(this->gradLambda(0), 0.0f);
  EXPECT_FLOAT_EQ(this->gradMu(0), 0.0f);
}

TYPED_TEST(DifferentiatorElasticElemTest, UniformFieldGradRhoEqualsVolume)
{
  // ux_fwd = 1, ux_dt2 = 1 on all nodes =>
  //   gradRho = ∫_Ω (ü†_x · u_x) dΩ = ∫_Ω 1 dΩ = 1.0 (unit cube)
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(this->gradRho(0), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradRhoAllThreeComponents)
{
  // All forward and dt2 components active:
  //   gradRho = ∫(ü†_x·u_x + ü†_y·u_y + ü†_z·u_z) dΩ
  // u_x=1,u_y=2,u_z=3 with all dt2=1  =>  gradRho = (1+2+3) * volume = 6.0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->uy_fwd(i) = 2.0f;
    this->uz_fwd(i) = 3.0f;
    this->ux_dt2(i) = 1.0f;
    this->uy_dt2(i) = 1.0f;
    this->uz_dt2(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(this->gradRho(0), 6.0f, 1e-4f);
}

TYPED_TEST(DifferentiatorElasticElemTest, ConstantFieldGradLambdaIsZero)
{
  // Constant displacement => ∇u = 0 => div(u†) = div(u) = 0 => gradLambda = 0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->uy_fwd(i) = 1.0f;
    this->uz_fwd(i) = 1.0f;
    this->ux_adj(i) = 1.0f;
    this->uy_adj(i) = 1.0f;
    this->uz_adj(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(this->gradLambda(0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, ConstantFieldGradMuIsZero)
{
  // Constant displacement => ε(u) = 0 => gradMu = 0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->uy_fwd(i) = 1.0f;
    this->uz_fwd(i) = 1.0f;
    this->ux_adj(i) = 1.0f;
    this->uy_adj(i) = 1.0f;
    this->uz_adj(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  EXPECT_NEAR(this->gradMu(0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradRhoScalesWithAmplitude)
{
  // Doubling ux_fwd doubles gradRho (linearity)
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);
  float single = this->gradRho(0);

  // Reset and double
  this->gradRho(0) = 0.0f;
  for (int i = 0; i < TestFixture::kNumNodes; ++i) this->ux_fwd(i) = 2.0f;

  this->makeDataAndCompute(diff);

  EXPECT_NEAR(this->gradRho(0), 2.0f * single, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest,
           ComputeAccumulatesIntoExistingGradient)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }
  this->gradRho(0) = 5.0f;

  typename TestFixture::Diff diff;
  this->makeDataAndCompute(diff);

  // 5.0 (initial) + 1.0 (volume of unit cube) = 6.0
  EXPECT_NEAR(this->gradRho(0), 6.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradLambdaNonZeroForLinearField)
{
  // Linear displacement u_x = x => div(u) = ∂u_x/∂x = 1
  // If both forward and adjoint have u_x = x:
  //   gradLambda = ∫ div(u†)·div(u) dΩ = ∫ 1·1 dΩ = 1.0
  //
  // On structure mesh [0,1]³ with nodes at x = i/ORDER:
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i)
      {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        this->ux_fwd(gIdx) = x;
        this->ux_adj(gIdx) = x;
      }

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  // ∫₀¹ ∫₀¹ ∫₀¹ 1·1 dx dy dz = 1.0
  // Tolerance is wider for Order 3 where uniform node spacing diverges from
  // GLL quadrature points, introducing discretization error.
  EXPECT_NEAR(this->gradLambda(0), 1.0f, 0.2f);
}

TYPED_TEST(DifferentiatorElasticElemTest, GradMuNonZeroForLinearField)
{
  // Linear displacement u_x = x => ε_xx = 1, all other ε_ij = 0
  // If both forward and adjoint have u_x = x:
  //   2·ε†:ε = 2·(1·1) = 2
  //   gradMu = ∫ 2 dΩ = 2.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i)
      {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        this->ux_fwd(gIdx) = x;
        this->ux_adj(gIdx) = x;
      }

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  // 2·ε†:ε = 2·ε_xx²  = 2·1 = 2; ∫_Ω 2 dΩ = 2.0
  EXPECT_NEAR(this->gradMu(0), 2.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemTest,
           LinearShearFieldProducesCorrectGradMu)
{
  // Pure shear: u_x = y, u_y = x  => ε_xy = 1, all diagonal ε = 0,
  //   div(u) = 0
  // 2·ε†:ε = 4·(ε_xy†·ε_xy) = 4·1 = 4
  // gradMu = ∫ 4 dΩ = 4.0
  // gradLambda = ∫ 0·0 dΩ = 0.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i)
      {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        float y = static_cast<float>(j) / TestFixture::kOrder;
        this->ux_fwd(gIdx) = y;
        this->uy_fwd(gIdx) = x;
        this->ux_adj(gIdx) = y;
        this->uy_adj(gIdx) = x;
      }

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradLambda(0), 0.0f, 1e-5f);
  EXPECT_NEAR(this->gradMu(0), 4.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemTest,
           PureDilatationLambdaMuConsistency)
{
  // Pure dilatation: u_x = x, u_y = y, u_z = z
  //   div(u) = 3, ε_xx = ε_yy = ε_zz = 1, off-diag = 0
  //   gradLambda = ∫ 3·3 dΩ = 9.0
  //   2·ε†:ε = 2·(1+1+1) = 6
  //   gradMu = ∫ 6 dΩ = 6.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i)
      {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        float y = static_cast<float>(j) / TestFixture::kOrder;
        float z = static_cast<float>(k) / TestFixture::kOrder;
        this->ux_fwd(gIdx) = x;
        this->uy_fwd(gIdx) = y;
        this->uz_fwd(gIdx) = z;
        this->ux_adj(gIdx) = x;
        this->uy_adj(gIdx) = y;
        this->uz_adj(gIdx) = z;
      }

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradLambda(0), 9.0f, 1.0f);
  EXPECT_NEAR(this->gradMu(0), 6.0f, 1.0f);
}

TYPED_TEST(DifferentiatorElasticElemTest,
           GradRhoIndependentOfStrainKernels)
{
  // gradRho depends only on dt2 and fwd, not on adjoint displacement gradients
  // Set ux_fwd=1, ux_dt2=1 but also set non-trivial adjoint displacement
  // => gradRho should still be 1.0
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i)
      {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        this->ux_fwd(gIdx) = 1.0f;
        this->ux_dt2(gIdx) = 1.0f;
        // Non-trivial adjoint displacement for strain kernels
        this->ux_adj(gIdx) = x;
        this->uy_adj(gIdx) = x;
      }

  typename TestFixture::Diff diff;

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  // gradRho only depends on dt2·fwd dot product: 1 * 1 * volume = 1.0
  EXPECT_NEAR(this->gradRho(0), 1.0f, 1e-5f);
}

// --- Polymorphic interface ---

TYPED_TEST(DifferentiatorElasticElemTest, PolymorphicInterface)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  std::unique_ptr<Differentiator> diff =
      std::make_unique<typename TestFixture::Diff>();
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  EXPECT_EQ(diff->getOrder(), TestFixture::kOrder);
  EXPECT_FALSE(diff->isModelOnNodes());

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  EXPECT_NO_THROW(diff->compute(mesh, data, 0.001f));
  EXPECT_GT(this->gradRho(0), 0.0f);
}

// =============================================================================
// DifferentiatorElasticNodeTest  (IS_MODEL_ON_NODES = true)
//
// Mesh: 1 element, (ORDER+1)^3 nodes.
// Gradients are node-indexed; contributions are scattered with ATOMICADD.
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorElasticNodeTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelStruct<float, int, kOrder>;
  using Integral =
      typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using DiffNode = DifferentiatorElastic<kOrder, Integral, Mesh, true>;
  using DiffElem = DifferentiatorElastic<kOrder, Integral, Mesh, false>;

  void SetUp() override
  {
    ux_fwd = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "ux_fwd");
    uy_fwd = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uy_fwd");
    uz_fwd = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uz_fwd");
    ux_adj = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "ux_adj");
    uy_adj = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uy_adj");
    uz_adj = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uz_adj");
    ux_dt2 = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "ux_dt2");
    uy_dt2 = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uy_dt2");
    uz_dt2 = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "uz_dt2");
    gradRho = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "gradRho");
    gradLambda = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "gradLambda");
    gradMu = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "gradMu");

    for (int i = 0; i < kNumNodes; ++i)
    {
      ux_fwd(i) = 0.0f;
      uy_fwd(i) = 0.0f;
      uz_fwd(i) = 0.0f;
      ux_adj(i) = 0.0f;
      uy_adj(i) = 0.0f;
      uz_adj(i) = 0.0f;
      ux_dt2(i) = 0.0f;
      uy_dt2(i) = 0.0f;
      uz_dt2(i) = 0.0f;
      gradRho(i) = 0.0f;
      gradLambda(i) = 0.0f;
      gradMu(i) = 0.0f;
    }
  }

  float sumGrad(VECTOR_REAL_VIEW const& v) const
  {
    float s = 0.0f;
    for (int i = 0; i < kNumNodes; ++i) s += v(i);
    return s;
  }

  float sumGradRho() const { return sumGrad(gradRho); }
  float sumGradLambda() const { return sumGrad(gradLambda); }
  float sumGradMu() const { return sumGrad(gradMu); }

  VECTOR_REAL_VIEW ux_fwd, uy_fwd, uz_fwd;
  VECTOR_REAL_VIEW ux_adj, uy_adj, uz_adj;
  VECTOR_REAL_VIEW ux_dt2, uy_dt2, uz_dt2;
  VECTOR_REAL_VIEW gradRho, gradLambda, gradMu;
};

TYPED_TEST_SUITE(DifferentiatorElasticNodeTest, OrderTypes);

// --- Static constants ---

TYPED_TEST(DifferentiatorElasticNodeTest, IsModelOnNodesConstant)
{
  EXPECT_TRUE(TestFixture::DiffNode::kIsModelOnNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorElasticNodeTest, GetOrderReturnsCorrectOrder)
{
  typename TestFixture::DiffNode diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticNodeTest, IsModelOnNodesReturnsTrue)
{
  typename TestFixture::DiffNode diff;
  EXPECT_TRUE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorElasticNodeTest, PrintDoesNotThrow)
{
  typename TestFixture::DiffNode diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorElasticNodeTest, ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->sumGradRho(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradLambda(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradMu(), 0.0f);
}

TYPED_TEST(DifferentiatorElasticNodeTest, UniformFieldGradRhoSumsToVolume)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradRho(), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticNodeTest,
           ConstantFieldGradLambdaSumsToZero)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_adj(i) = 1.0f;
  }

  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);
  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradLambda(), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticNodeTest,
           NodeBasedSumEqualsElementBasedResult)
{
  // For single-element mesh, ∑ node_gradients == elem_gradient
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();
  int const dim = TestFixture::kOrder + 1;
  for (int k = 0; k < dim; ++k)
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < dim; ++i)
      {
        int gIdx = mesh.globalNodeIndex(0, i, j, k);
        float x = static_cast<float>(i) / TestFixture::kOrder;
        float y = static_cast<float>(j) / TestFixture::kOrder;
        this->ux_fwd(gIdx) = x;
        this->uy_fwd(gIdx) = y;
        this->ux_adj(gIdx) = x;
        this->uy_adj(gIdx) = y;
        this->ux_dt2(gIdx) = 1.0f;
      }

  // Node-based
  typename TestFixture::DiffNode diffNode;
  {
    WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
    WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                     this->ux_dt2, this->uy_dt2, this->uz_dt2);
    GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
    GradientDataElastic data(fwd, bwd, grad);
    diffNode.compute(mesh, data, 0.001f);
  }
  float nodeRho = this->sumGradRho();
  float nodeLambda = this->sumGradLambda();
  float nodeMu = this->sumGradMu();

  // Element-based
  auto gradRhoElem = allocateVector<VECTOR_REAL_VIEW>(1, "gradRhoElem");
  auto gradLambdaElem =
      allocateVector<VECTOR_REAL_VIEW>(1, "gradLambdaElem");
  auto gradMuElem = allocateVector<VECTOR_REAL_VIEW>(1, "gradMuElem");
  gradRhoElem(0) = 0.0f;
  gradLambdaElem(0) = 0.0f;
  gradMuElem(0) = 0.0f;

  typename TestFixture::DiffElem diffElem;
  {
    WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
    WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                     this->ux_dt2, this->uy_dt2, this->uz_dt2);
    GradientElastic grad(gradRhoElem, gradLambdaElem, gradMuElem);
    GradientDataElastic data(fwd, bwd, grad);
    diffElem.compute(mesh, data, 0.001f);
  }

  EXPECT_NEAR(nodeRho, gradRhoElem(0), 1e-4f);
  EXPECT_NEAR(nodeLambda, gradLambdaElem(0), 1e-4f);
  EXPECT_NEAR(nodeMu, gradMuElem(0), 1e-4f);
}

// --- Polymorphic interface ---

TYPED_TEST(DifferentiatorElasticNodeTest, PolymorphicInterface)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  std::unique_ptr<Differentiator> diff =
      std::make_unique<typename TestFixture::DiffNode>();
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  EXPECT_EQ(diff->getOrder(), TestFixture::kOrder);
  EXPECT_TRUE(diff->isModelOnNodes());

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  EXPECT_NO_THROW(diff->compute(mesh, data, 0.001f));
  EXPECT_GT(this->sumGradRho(), 0.0f);
}

}  // namespace test
}  // namespace gradient

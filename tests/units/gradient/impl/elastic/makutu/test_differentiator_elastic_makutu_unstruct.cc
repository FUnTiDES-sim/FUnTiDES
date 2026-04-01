#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "data_type.h"
#include "differentiator_elastic.h"
#include "differentiator_data_elastic.h"
#include "fe/Integrals.hpp"
#include "model_unstruct.h"

namespace gradient
{
namespace test
{

// =============================================================================
// Order wrappers
// =============================================================================

struct Order1U
{
  static constexpr int kOrder = 1;
};
struct Order2U
{
  static constexpr int kOrder = 2;
};
struct Order3U
{
  static constexpr int kOrder = 3;
};

using OrderTypesU = ::testing::Types<Order1U, Order2U, Order3U>;

// =============================================================================
// Helper: 1-element unstructured unit cube mesh [0,1]³
// =============================================================================

template <int ORDER>
static model::ModelUnstruct<float, int> makeUnstructMesh1x1x1()
{
  constexpr int npe = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  model::ModelUnstructData<float, int> data;
  data.order_ = ORDER;
  data.n_element_ = 1;
  data.n_node_ = npe;
  data.lx_ = data.ly_ = data.lz_ = 1.0f;
  data.isModelOnNodes_ = false;
  data.isElastic_ = true;

  data.global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(1, npe, "gni");
  data.nodes_coords_x_ = allocateVector<VECTOR_REAL_VIEW>(npe, "cx");
  data.nodes_coords_y_ = allocateVector<VECTOR_REAL_VIEW>(npe, "cy");
  data.nodes_coords_z_ = allocateVector<VECTOR_REAL_VIEW>(npe, "cz");

  for (int lDof = 0; lDof < npe; ++lDof)
  {
    data.global_node_index_(0, lDof) = lDof;

    int k = lDof / ((ORDER + 1) * (ORDER + 1));
    int j = (lDof / (ORDER + 1)) % (ORDER + 1);
    int i = lDof % (ORDER + 1);

    data.nodes_coords_x_[lDof] = static_cast<float>(i) / ORDER;
    data.nodes_coords_y_[lDof] = static_cast<float>(j) / ORDER;
    data.nodes_coords_z_[lDof] = static_cast<float>(k) / ORDER;
  }

  return model::ModelUnstruct<float, int>(data);
}

// =============================================================================
// DifferentiatorElasticElemUnstructTest  (IS_MODEL_ON_NODES = false)
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorElasticElemUnstructTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelUnstruct<float, int>;
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

  VECTOR_REAL_VIEW ux_fwd, uy_fwd, uz_fwd;
  VECTOR_REAL_VIEW ux_adj, uy_adj, uz_adj;
  VECTOR_REAL_VIEW ux_dt2, uy_dt2, uz_dt2;
  VECTOR_REAL_VIEW gradRho, gradLambda, gradMu;
};

TYPED_TEST_SUITE(DifferentiatorElasticElemUnstructTest, OrderTypesU);

// --- Static constants ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest, OrderConstant)
{
  EXPECT_EQ(TestFixture::Diff::kOrder, TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, IsModelOnNodesConstant)
{
  EXPECT_FALSE(TestFixture::Diff::kIsModelOnNodes);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, PointsPerElementConstant)
{
  EXPECT_EQ(TestFixture::Diff::kPointsPerElement, TestFixture::kNumNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest, GetOrderReturnsCorrectOrder)
{
  typename TestFixture::Diff diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest, IsModelOnNodesReturnsFalse)
{
  typename TestFixture::Diff diff;
  EXPECT_FALSE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest, PrintDoesNotThrow)
{
  typename TestFixture::Diff diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->gradRho(0), 0.0f);
  EXPECT_FLOAT_EQ(this->gradLambda(0), 0.0f);
  EXPECT_FLOAT_EQ(this->gradMu(0), 0.0f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           UniformFieldGradRhoEqualsVolume)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradRho(0), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           ConstantFieldGradLambdaIsZero)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_adj(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradLambda(0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           ConstantFieldGradMuIsZero)
{
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
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradMu(0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           GradLambdaNonZeroForLinearField)
{
  // u_x = x => div(u) = 1 => gradLambda = ∫ 1 dΩ = 1.0
  constexpr int dim = TestFixture::kOrder + 1;
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof)
  {
    int i = lDof % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    this->ux_fwd(lDof) = x;
    this->ux_adj(lDof) = x;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradLambda(0), 1.0f, 0.2f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           GradMuNonZeroForLinearField)
{
  // u_x = x => ε_xx = 1 => 2·ε†:ε = 2 => gradMu = 2.0
  constexpr int dim = TestFixture::kOrder + 1;
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof)
  {
    int i = lDof % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    this->ux_fwd(lDof) = x;
    this->ux_adj(lDof) = x;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradMu(0), 2.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           PureDilatationLambdaMuConsistency)
{
  // u_x=x, u_y=y, u_z=z => div=3, ε_xx=ε_yy=ε_zz=1
  // gradLambda = 9.0, gradMu = 6.0
  constexpr int dim = TestFixture::kOrder + 1;
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof)
  {
    int k = lDof / (dim * dim);
    int j = (lDof / dim) % dim;
    int i = lDof % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    float y = static_cast<float>(j) / TestFixture::kOrder;
    float z = static_cast<float>(k) / TestFixture::kOrder;
    this->ux_fwd(lDof) = x;
    this->uy_fwd(lDof) = y;
    this->uz_fwd(lDof) = z;
    this->ux_adj(lDof) = x;
    this->uy_adj(lDof) = y;
    this->uz_adj(lDof) = z;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradLambda(0), 9.0f, 1.0f);
  EXPECT_NEAR(this->gradMu(0), 6.0f, 1.0f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           LinearShearFieldProducesCorrectGradMu)
{
  // Pure shear: u_x=y, u_y=x => div=0 => gradLambda=0
  // ε_xy = 1 => 2·ε:ε = 4·(ε_xy²) = 4 => gradMu = 4.0
  constexpr int dim = TestFixture::kOrder + 1;
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof)
  {
    int i = lDof % dim;
    int j = (lDof / dim) % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    float y = static_cast<float>(j) / TestFixture::kOrder;
    this->ux_fwd(lDof) = y;
    this->uy_fwd(lDof) = x;
    this->ux_adj(lDof) = y;
    this->uy_adj(lDof) = x;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradLambda(0), 0.0f, 1e-5f);
  EXPECT_NEAR(this->gradMu(0), 4.0f, 0.4f);
}

TYPED_TEST(DifferentiatorElasticElemUnstructTest,
           ComputeAccumulatesIntoExistingGradient)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }
  this->gradRho(0) = 5.0f;

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradRho(0), 6.0f, 1e-5f);
}

// --- Node-based unstructured ---

template <typename OrderWrapper>
class DifferentiatorElasticNodeUnstructTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);

  using Mesh = model::ModelUnstruct<float, int>;
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

TYPED_TEST_SUITE(DifferentiatorElasticNodeUnstructTest, OrderTypesU);

TYPED_TEST(DifferentiatorElasticNodeUnstructTest,
           ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::DiffNode diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

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

TYPED_TEST(DifferentiatorElasticNodeUnstructTest,
           UniformFieldGradRhoSumsToVolume)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  typename TestFixture::DiffNode diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardElastic fwd(this->ux_fwd, this->uy_fwd, this->uz_fwd);
  WavefieldViewBackwardElastic bwd(this->ux_adj, this->uy_adj, this->uz_adj,
                                   this->ux_dt2, this->uy_dt2, this->uz_dt2);
  GradientElastic grad(this->gradRho, this->gradLambda, this->gradMu);
  GradientDataElastic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradRho(), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorElasticNodeUnstructTest,
           NodeBasedSumEqualsElementBasedResult)
{
  constexpr int dim = TestFixture::kOrder + 1;
  for (int lDof = 0; lDof < TestFixture::kNumNodes; ++lDof)
  {
    int i = lDof % dim;
    int j = (lDof / dim) % dim;
    float x = static_cast<float>(i) / TestFixture::kOrder;
    float y = static_cast<float>(j) / TestFixture::kOrder;
    this->ux_fwd(lDof) = x;
    this->uy_fwd(lDof) = y;
    this->ux_adj(lDof) = x;
    this->uy_adj(lDof) = y;
    this->ux_dt2(lDof) = 1.0f;
  }

  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

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

TYPED_TEST(DifferentiatorElasticNodeUnstructTest, PolymorphicInterface)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->ux_fwd(i) = 1.0f;
    this->ux_dt2(i) = 1.0f;
  }

  std::unique_ptr<Differentiator> diff =
      std::make_unique<typename TestFixture::DiffNode>();
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

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

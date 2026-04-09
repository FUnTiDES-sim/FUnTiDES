#include <gtest/gtest.h>

#include <memory>

#include "data_type.h"
#include "differentiator_acoustic.h"
#include "differentiator_data_acoustic.h"
#include "Integrals.h"
#include "model_unstruct.h"

namespace gradient
{
namespace test
{

// =============================================================================
// Order wrappers — Google Test requires types, not non-type parameters.
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
//
// Nodes are laid out on a regular (ORDER+1)^3 grid.
// Corner nodes are always at {0,1}^3, so the element volume = 1.
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
  data.isElastic_ = false;

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
// DifferentiatorAcousticElemUnstructTest  (IS_MODEL_ON_NODES = false)
//
// Mesh: 1 unstructured element, (ORDER+1)^3 nodes.
// Gradients are element-indexed (size = 1).
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorAcousticElemUnstructTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelUnstruct<float, int>;
  using Integral =
      typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using Diff = DifferentiatorAcoustic<kOrder, Integral, Mesh, false>;

  void SetUp() override
  {
    pn = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "pn");
    qn = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "qn");
    qnPrev = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "qnPrev");
    qnPrevPrev = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "qnPrevPrev");
    gradKappa = allocateVector<VECTOR_REAL_VIEW>(kNumElements, "gradKappa");
    gradBuoyancy =
        allocateVector<VECTOR_REAL_VIEW>(kNumElements, "gradBuoyancy");

    for (int i = 0; i < kNumNodes; ++i)
    {
      pn(i) = 0.0f;
      qn(i) = 0.0f;
      qnPrev(i) = 0.0f;
      qnPrevPrev(i) = 0.0f;
    }
    gradKappa(0) = 0.0f;
    gradBuoyancy(0) = 0.0f;
  }

  VECTOR_REAL_VIEW pn, qn, qnPrev, qnPrevPrev;
  VECTOR_REAL_VIEW gradKappa, gradBuoyancy;
};

TYPED_TEST_SUITE(DifferentiatorAcousticElemUnstructTest, OrderTypesU);

// --- Static constants ---

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, OrderConstant)
{
  EXPECT_EQ(TestFixture::Diff::kOrder, TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, IsModelOnNodesConstant)
{
  EXPECT_FALSE(TestFixture::Diff::kIsModelOnNodes);
}

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, PointsPerElementConstant)
{
  EXPECT_EQ(TestFixture::Diff::kPointsPerElement, TestFixture::kNumNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, GetOrderReturnsCorrectOrder)
{
  typename TestFixture::Diff diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, IsModelOnNodesReturnsFalse)
{
  typename TestFixture::Diff diff;
  EXPECT_FALSE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, PrintDoesNotThrow)
{
  typename TestFixture::Diff diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorAcousticElemUnstructTest,
           ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->gradKappa(0), 0.0f);
  EXPECT_FLOAT_EQ(this->gradBuoyancy(0), 0.0f);
}

TYPED_TEST(DifferentiatorAcousticElemUnstructTest,
           UniformFieldGradKappaEqualsVolume)
{
  // pn = 1, qn = 1e-6 (qdt2 = 1 with dt=0.001)  =>  gradKappa = ∫_Ω 1 dΩ = 1.0 (unit cube)
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f; this->qnPrev(i) = 0.0f; this->qnPrevPrev(i) = 0.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradKappa(0), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticElemUnstructTest,
           ConstantFieldGradBuoyancyIsZero)
{
  // pn = qn = const  =>  ∇pn = ∇qn = 0  =>  stiffness contribution = 0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradBuoyancy(0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, GradKappaScalesWithAmplitude)
{
  // Doubling pn doubles gradKappa
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f; this->qnPrev(i) = 0.0f; this->qnPrevPrev(i) = 0.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  {
    WavefieldViewForwardAcoustic fwd(this->pn);
    WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
    GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
    GradientDataAcoustic data(fwd, bwd, grad);
    diff.compute(mesh, data, 0.001f);
  }
  float single = this->gradKappa(0);

  this->gradKappa(0) = 0.0f;
  for (int i = 0; i < TestFixture::kNumNodes; ++i) this->pn(i) = 2.0f;

  {
    WavefieldViewForwardAcoustic fwd(this->pn);
    WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
    GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
    GradientDataAcoustic data(fwd, bwd, grad);
    diff.compute(mesh, data, 0.001f);
  }

  EXPECT_NEAR(this->gradKappa(0), 2.0f * single, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticElemUnstructTest,
           ComputeAccumulatesIntoExistingGradient)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f; this->qnPrev(i) = 0.0f; this->qnPrevPrev(i) = 0.0f;
  }
  this->gradKappa(0) = 5.0f;

  typename TestFixture::Diff diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  // 5.0 (initial) + 1.0 (volume of unit cube) = 6.0
  EXPECT_NEAR(this->gradKappa(0), 6.0f, 1e-5f);
}

// --- Polymorphic interface ---

TYPED_TEST(DifferentiatorAcousticElemUnstructTest, PolymorphicInterface)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f; this->qnPrev(i) = 0.0f; this->qnPrevPrev(i) = 0.0f;
  }

  std::unique_ptr<Differentiator> diff =
      std::make_unique<typename TestFixture::Diff>();
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  EXPECT_EQ(diff->getOrder(), TestFixture::kOrder);
  EXPECT_FALSE(diff->isModelOnNodes());

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  EXPECT_NO_THROW(diff->compute(mesh, data, 0.001f));
  EXPECT_GT(this->gradKappa(0), 0.0f);
}

// =============================================================================
// DifferentiatorAcousticNodeUnstructTest  (IS_MODEL_ON_NODES = true)
//
// Mesh: 1 unstructured element, (ORDER+1)^3 nodes.
// Gradients are node-indexed; contributions are scattered atomically.
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorAcousticNodeUnstructTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelUnstruct<float, int>;
  using Integral =
      typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
  using DiffNode = DifferentiatorAcoustic<kOrder, Integral, Mesh, true>;
  using DiffElem = DifferentiatorAcoustic<kOrder, Integral, Mesh, false>;

  void SetUp() override
  {
    pn = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "pn");
    qn = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "qn");
    qnPrev = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "qnPrev");
    qnPrevPrev = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "qnPrevPrev");
    gradKappa = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "gradKappa");
    gradBuoyancy = allocateVector<VECTOR_REAL_VIEW>(kNumNodes, "gradBuoyancy");

    for (int i = 0; i < kNumNodes; ++i)
    {
      pn(i) = 0.0f;
      qn(i) = 0.0f;
      qnPrev(i) = 0.0f;
      qnPrevPrev(i) = 0.0f;
      gradKappa(i) = 0.0f;
      gradBuoyancy(i) = 0.0f;
    }
  }

  float sumGradKappa() const
  {
    float s = 0.0f;
    for (int i = 0; i < kNumNodes; ++i) s += gradKappa(i);
    return s;
  }

  float sumGradBuoyancy() const
  {
    float s = 0.0f;
    for (int i = 0; i < kNumNodes; ++i) s += gradBuoyancy(i);
    return s;
  }

  VECTOR_REAL_VIEW pn, qn, qnPrev, qnPrevPrev;
  VECTOR_REAL_VIEW gradKappa, gradBuoyancy;
};

TYPED_TEST_SUITE(DifferentiatorAcousticNodeUnstructTest, OrderTypesU);

// --- Static constants ---

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest, IsModelOnNodesConstant)
{
  EXPECT_TRUE(TestFixture::DiffNode::kIsModelOnNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest, GetOrderReturnsCorrectOrder)
{
  typename TestFixture::DiffNode diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest, IsModelOnNodesReturnsTrue)
{
  typename TestFixture::DiffNode diff;
  EXPECT_TRUE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest, PrintDoesNotThrow)
{
  typename TestFixture::DiffNode diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest,
           ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::DiffNode diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->sumGradKappa(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradBuoyancy(), 0.0f);
}

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest,
           UniformFieldGradKappaSumsToVolume)
{
  // pn = 1, qn = 1e-6 (qdt2 = 1 with dt=0.001)  =>  scattered contributions must sum to ∫_Ω 1 dΩ = 1.0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f; this->qnPrev(i) = 0.0f; this->qnPrevPrev(i) = 0.0f;
  }

  typename TestFixture::DiffNode diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradKappa(), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest,
           ConstantFieldGradBuoyancySumsToZero)
{
  // pn = qn = 1  =>  ∇p = ∇q = 0  =>  stiffness contribution sums to 0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1.0f;
  }

  typename TestFixture::DiffNode diff;
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradBuoyancy(), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest,
           NodeBasedSumEqualsElementBasedResult)
{
  // For a single-element mesh (no shared nodes), the sum of node-scattered
  // gradients must equal the single element-accumulated value.
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f; this->qnPrev(i) = 0.0f; this->qnPrevPrev(i) = 0.0f;
  }

  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);

  // Node-based
  typename TestFixture::DiffNode diffNode;
  GradientAcoustic gradNode(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic dataNode(fwd, bwd, gradNode);
  diffNode.compute(mesh, dataNode, 0.001f);
  float nodeSum = this->sumGradKappa();

  // Element-based
  auto gradKappaElem = allocateVector<VECTOR_REAL_VIEW>(1, "gradKappaElem");
  auto gradBuoyancyElem =
      allocateVector<VECTOR_REAL_VIEW>(1, "gradBuoyancyElem");
  gradKappaElem(0) = 0.0f;
  gradBuoyancyElem(0) = 0.0f;

  typename TestFixture::DiffElem diffElem;
  GradientAcoustic gradElem(gradKappaElem, gradBuoyancyElem);
  GradientDataAcoustic dataElem(fwd, bwd, gradElem);
  diffElem.compute(mesh, dataElem, 0.001f);

  EXPECT_NEAR(nodeSum, gradKappaElem(0), 1e-5f);
}

// --- Polymorphic interface ---

TYPED_TEST(DifferentiatorAcousticNodeUnstructTest, PolymorphicInterface)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f; this->qnPrev(i) = 0.0f; this->qnPrevPrev(i) = 0.0f;
  }

  std::unique_ptr<Differentiator> diff =
      std::make_unique<typename TestFixture::DiffNode>();
  auto mesh = makeUnstructMesh1x1x1<TestFixture::kOrder>();

  EXPECT_EQ(diff->getOrder(), TestFixture::kOrder);
  EXPECT_TRUE(diff->isModelOnNodes());

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  EXPECT_NO_THROW(diff->compute(mesh, data, 0.001f));
  EXPECT_GT(this->sumGradKappa(), 0.0f);
}

}  // namespace test
}  // namespace gradient

#include <gtest/gtest.h>

#include <memory>
#include <numeric>

#include "data_type.h"
#include "differentiator_acoustic.h"
#include "differentiator_data_acoustic.h"
#include "fe/Integrals.hpp"
#include "model_struct.h"

namespace gradient
{
namespace test
{

// =============================================================================
// Order wrappers — Google Test requires types, not non-type parameters.
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
  data.isElastic_ = false;
  return model::ModelStruct<float, int, ORDER>(data);
}

// =============================================================================
// DifferentiatorAcousticElemTest  (IS_MODEL_ON_NODES = false)
//
// Mesh: 1 element, (ORDER+1)^3 nodes.
// Gradients are accumulated, element-indexed (size = 1).
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorAcousticElemTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelStruct<float, int, kOrder>;
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

TYPED_TEST_SUITE(DifferentiatorAcousticElemTest, OrderTypes);

// --- Static constants ---

TYPED_TEST(DifferentiatorAcousticElemTest, OrderConstant)
{
  EXPECT_EQ(TestFixture::Diff::kOrder, TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorAcousticElemTest, IsModelOnNodesConstant)
{
  EXPECT_FALSE(TestFixture::Diff::kIsModelOnNodes);
}

TYPED_TEST(DifferentiatorAcousticElemTest, PointsPerElementConstant)
{
  EXPECT_EQ(TestFixture::Diff::kPointsPerElement, TestFixture::kNumNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorAcousticElemTest, GetOrderReturnsCorrectOrder)
{
  typename TestFixture::Diff diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorAcousticElemTest, IsModelOnNodesReturnsFalse)
{
  typename TestFixture::Diff diff;
  EXPECT_FALSE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorAcousticElemTest, PrintDoesNotThrow)
{
  typename TestFixture::Diff diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorAcousticElemTest, ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::Diff diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->gradKappa(0), 0.0f);
  EXPECT_FLOAT_EQ(this->gradBuoyancy(0), 0.0f);
}

TYPED_TEST(DifferentiatorAcousticElemTest, UniformFieldGradKappaEqualsVolume)
{
  // pn = 1, qn = 1e-6 (qdt2 = 1 with dt=0.001)  =>  gradKappa = ∑_q
  // mass_weight_q = ∫_Ω 1 dΩ = 1.0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f;
    this->qnPrev(i) = 0.0f;
    this->qnPrevPrev(i) = 0.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradKappa(0), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticElemTest, ConstantFieldGradBuoyancyIsZero)
{
  // pn = qn = const  =>  ∇pn = ∇qn = 0  =>  stiffness contribution = 0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->gradBuoyancy(0), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticElemTest, GradKappaScalesWithAmplitude)
{
  // Doubling pn doubles gradKappa
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f;
    this->qnPrev(i) = 0.0f;
    this->qnPrevPrev(i) = 0.0f;
  }

  typename TestFixture::Diff diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

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

TYPED_TEST(DifferentiatorAcousticElemTest,
           ComputeAccumulatesIntoExistingGradient)
{
  // gradKappa starts at some non-zero value; compute() must add to it
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f;
    this->qnPrev(i) = 0.0f;
    this->qnPrevPrev(i) = 0.0f;
  }
  this->gradKappa(0) = 5.0f;

  typename TestFixture::Diff diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  // 5.0 (initial) + 1.0 (volume of unit cube) = 6.0
  EXPECT_NEAR(this->gradKappa(0), 6.0f, 1e-5f);
}

// --- Polymorphic interface ---

TYPED_TEST(DifferentiatorAcousticElemTest, PolymorphicInterface)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f;
    this->qnPrev(i) = 0.0f;
    this->qnPrevPrev(i) = 0.0f;
  }

  std::unique_ptr<Differentiator> diff =
      std::make_unique<typename TestFixture::Diff>();
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

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
// DifferentiatorAcousticNodeTest  (IS_MODEL_ON_NODES = true)
//
// Mesh: 1 element, (ORDER+1)^3 nodes.
// Gradients are accumulated, node-indexed; contributions are scattered
// atomically from each element thread.
// =============================================================================

template <typename OrderWrapper>
class DifferentiatorAcousticNodeTest : public ::testing::Test
{
 protected:
  static constexpr int kOrder = OrderWrapper::kOrder;
  static constexpr int kNumNodes = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
  static constexpr int kNumElements = 1;

  using Mesh = model::ModelStruct<float, int, kOrder>;
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

TYPED_TEST_SUITE(DifferentiatorAcousticNodeTest, OrderTypes);

// --- Static constants ---

TYPED_TEST(DifferentiatorAcousticNodeTest, IsModelOnNodesConstant)
{
  EXPECT_TRUE(TestFixture::DiffNode::kIsModelOnNodes);
}

// --- Virtual getters ---

TYPED_TEST(DifferentiatorAcousticNodeTest, GetOrderReturnsCorrectOrder)
{
  typename TestFixture::DiffNode diff;
  EXPECT_EQ(diff.getOrder(), TestFixture::kOrder);
}

TYPED_TEST(DifferentiatorAcousticNodeTest, IsModelOnNodesReturnsTrue)
{
  typename TestFixture::DiffNode diff;
  EXPECT_TRUE(diff.isModelOnNodes());
}

// --- Print ---

TYPED_TEST(DifferentiatorAcousticNodeTest, PrintDoesNotThrow)
{
  typename TestFixture::DiffNode diff;
  EXPECT_NO_THROW(diff.print());
}

// --- compute() correctness ---

TYPED_TEST(DifferentiatorAcousticNodeTest, ZeroWavefieldsYieldZeroGradients)
{
  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_FLOAT_EQ(this->sumGradKappa(), 0.0f);
  EXPECT_FLOAT_EQ(this->sumGradBuoyancy(), 0.0f);
}

TYPED_TEST(DifferentiatorAcousticNodeTest, UniformFieldGradKappaSumsToVolume)
{
  // pn = 1, qn = 1e-6 (qdt2 = 1 with dt=0.001) on all nodes.
  // Node-scattered contributions must sum to ∫_Ω 1 dΩ = 1.0.
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f;
    this->qnPrev(i) = 0.0f;
    this->qnPrevPrev(i) = 0.0f;
  }

  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradKappa(), 1.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticNodeTest, ConstantFieldGradBuoyancySumsToZero)
{
  // pn = qn = 1  =>  ∇p = ∇q = 0  =>  stiffness contribution sums to 0
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1.0f;
  }

  typename TestFixture::DiffNode diff;
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

  WavefieldViewForwardAcoustic fwd(this->pn);
  WavefieldViewBackwardAcoustic bwd(this->qn, this->qnPrev, this->qnPrevPrev);
  GradientAcoustic grad(this->gradKappa, this->gradBuoyancy);
  GradientDataAcoustic data(fwd, bwd, grad);

  diff.compute(mesh, data, 0.001f);

  EXPECT_NEAR(this->sumGradBuoyancy(), 0.0f, 1e-5f);
}

TYPED_TEST(DifferentiatorAcousticNodeTest, NodeBasedSumEqualsElementBasedResult)
{
  // For a single-element mesh (no shared nodes), the sum of node-scattered
  // gradients must equal the single element-accumulated value.
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f;
    this->qnPrev(i) = 0.0f;
    this->qnPrevPrev(i) = 0.0f;
  }

  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

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

TYPED_TEST(DifferentiatorAcousticNodeTest, PolymorphicInterface)
{
  for (int i = 0; i < TestFixture::kNumNodes; ++i)
  {
    this->pn(i) = 1.0f;
    this->qn(i) = 1e-6f;
    this->qnPrev(i) = 0.0f;
    this->qnPrevPrev(i) = 0.0f;
  }

  std::unique_ptr<Differentiator> diff =
      std::make_unique<typename TestFixture::DiffNode>();
  auto mesh = makeMesh1x1x1<TestFixture::kOrder>();

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
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "data_type.h"
#include "gllpoints.h"
#include "sem_2d_local_time_step_plan.h"

namespace solver {
namespace test {
namespace {

constexpr float kVp = 1500.f;
constexpr float kCfl = 0.3f;

// A conforming strip of rectangular elements, with possibly nonuniform widths.
model::ModelUnstructData<float, int> strip(int order, std::vector<float> const& edges, bool nodal = false) {
  int const ex = static_cast<int>(edges.size()) - 1;
  int const nx = ex * order + 1;
  int const n1d = order + 1;
  int const nn = nx * n1d;
  model::ModelUnstructData<float, int> d;
  d.dim_ = 2;
  d.order_ = order;
  d.n_element_ = ex;
  d.n_node_ = nn;
  d.lx_ = edges.back() - edges.front();
  d.ly_ = 0;
  d.lz_ = 100;
  d.ox_ = d.oy_ = d.oz_ = 0;
  d.isModelOnNodes_ = nodal;
  d.isElastic_ = false;
  d.global_node_index_ = allocateArray2D<dualArrayInt>(ex, n1d * n1d, "strip nodes");
  d.nodes_coords_x_ = allocateVector<dualVectorReal>(nn, "strip x");
  d.nodes_coords_y_ = allocateVector<dualVectorReal>(nn, "strip y");
  d.nodes_coords_z_ = allocateVector<dualVectorReal>(nn, "strip z");
  for (int e = 0; e < ex; ++e)
    for (int k = 0; k < n1d; ++k)
      for (int i = 0; i < n1d; ++i) {
        int const node = e * order + i + k * nx;
        d.global_node_index_.view_host()(e, i + k * n1d) = node;
        d.nodes_coords_x_.view_host()(node) = edges[e] + (edges[e + 1] - edges[e]) * (GLLPoints::get(order, i) + 1) / 2;
        d.nodes_coords_z_.view_host()(node) = 50 * (GLLPoints::get(order, k) + 1);
      }
  auto& vp = nodal ? d.model_vp_node_ : d.model_vp_element_;
  auto& rho = nodal ? d.model_rho_node_ : d.model_rho_element_;
  vp = allocateVector<dualVectorReal>(nodal ? nn : ex, "strip vp");
  rho = allocateVector<dualVectorReal>(nodal ? nn : ex, "strip rho");
  Kokkos::deep_copy(vp.view_host(), kVp);
  Kokkos::deep_copy(rho.view_host(), 1.f);
  for (auto* view : {&d.nodes_coords_x_, &d.nodes_coords_y_, &d.nodes_coords_z_, &vp, &rho}) view->modify_host();
  d.global_node_index_.modify_host();
  return d;
}

// On a 100 x 100 element the smallest spacing is the first GLL gap times 50.
float uniformEstimate(int order) { return kCfl * 50 * (GLLPoints::get(order, 1) + 1) / kVp; }

template <typename T>
class SEM2DLocalTimeStepPlanTest : public ::testing::Test {};
using Orders = ::testing::Types<std::integral_constant<int, 1>, std::integral_constant<int, 2>,
                                std::integral_constant<int, 4>, std::integral_constant<int, 9>>;
TYPED_TEST_SUITE(SEM2DLocalTimeStepPlanTest, Orders);

TYPED_TEST(SEM2DLocalTimeStepPlanTest, UniformMeshNeedsNoFineRegion) {
  constexpr int order = TypeParam::value;
  model::ModelUnstruct<float, int> mesh(strip(order, {0, 100, 200}));
  auto const plan = buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.9f * uniformEstimate(order));
  ASSERT_EQ(plan.element_steps.size(), 2u);
  for (auto const& step : plan.element_steps) {
    EXPECT_NEAR(step.estimate, uniformEstimate(order), 1e-5f * uniformEstimate(order));
    EXPECT_EQ(step.max_vp, kVp);
  }
  EXPECT_EQ(plan.coarse_elements, (std::vector<int>{0, 1}));
  EXPECT_EQ(plan.coarse_nodes.size(), static_cast<size_t>(mesh.getNumberOfNodes()));
  EXPECT_TRUE(plan.fine_elements.empty());
  EXPECT_TRUE(plan.fine_nodes.empty());
  EXPECT_TRUE(plan.interface_nodes.empty());
  EXPECT_TRUE(plan.fine_support_elements.empty());
}

TYPED_TEST(SEM2DLocalTimeStepPlanTest, AllFineMeshHasNoCoarseInterface) {
  constexpr int order = TypeParam::value;
  model::ModelUnstruct<float, int> mesh(strip(order, {0, 100, 200}));
  auto const plan = buildSEM2DLocalTimeStepPlan(mesh, kCfl, 1.5f * uniformEstimate(order));
  EXPECT_EQ(plan.fine_dt, plan.coarse_dt / 2);
  EXPECT_TRUE(plan.coarse_elements.empty());
  EXPECT_TRUE(plan.coarse_nodes.empty());
  EXPECT_EQ(plan.fine_elements, (std::vector<int>{0, 1}));
  EXPECT_EQ(plan.fine_nodes.size(), static_cast<size_t>(mesh.getNumberOfNodes()));
  EXPECT_TRUE(plan.interface_nodes.empty());
  EXPECT_EQ(plan.fine_support_elements, plan.fine_elements);
}

TYPED_TEST(SEM2DLocalTimeStepPlanTest, RefinedCellIncludesCoarseNeighboursInSupport) {
  constexpr int order = TypeParam::value;
  model::ModelUnstruct<float, int> mesh(strip(order, {0, 100, 200, 250, 350, 450}));
  auto const plan = buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.9f * uniformEstimate(order));
  EXPECT_EQ(plan.coarse_elements, (std::vector<int>{0, 1, 3, 4}));
  EXPECT_EQ(plan.fine_elements, (std::vector<int>{2}));
  EXPECT_EQ(plan.fine_support_elements, (std::vector<int>{1, 2, 3}));
  EXPECT_EQ(plan.fine_nodes.size(), static_cast<size_t>((order + 1) * (order + 1)));
  EXPECT_EQ(plan.interface_nodes.size(), static_cast<size_t>(2 * (order + 1)));
  for (int n : plan.interface_nodes) EXPECT_TRUE(mesh.nodeCoord(n, 0) == 200 || mesh.nodeCoord(n, 0) == 250);
  EXPECT_EQ(plan.coarse_nodes.size() + plan.fine_nodes.size(), static_cast<size_t>(mesh.getNumberOfNodes()));
  EXPECT_NEAR(plan.element_steps[2].estimate, uniformEstimate(order) / 2, 1e-5f * uniformEstimate(order));
}

TYPED_TEST(SEM2DLocalTimeStepPlanTest, FasterMaterialRequiresFineSteps) {
  constexpr int order = TypeParam::value;
  auto d = strip(order, {0, 100, 200, 300});
  d.model_vp_element_.view_host()(1) = 2 * kVp;
  model::ModelUnstruct<float, int> mesh(d);
  auto const plan = buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.9f * uniformEstimate(order));
  EXPECT_EQ(plan.fine_elements, (std::vector<int>{1}));
  EXPECT_EQ(plan.fine_support_elements, (std::vector<int>{0, 1, 2}));
  EXPECT_EQ(plan.element_steps[1].max_vp, 2 * kVp);
  EXPECT_NEAR(plan.element_steps[1].estimate, plan.element_steps[0].estimate / 2, 1e-6f);
}

TYPED_TEST(SEM2DLocalTimeStepPlanTest, NodalModelUsesFastestNodeOfEachElement) {
  constexpr int order = TypeParam::value;
  auto d = strip(order, {0, 100, 200}, true);
  d.model_vp_node_.view_host()(2 * order) = 2 * kVp;  // rightmost bottom vertex
  model::ModelUnstruct<float, int> mesh(d);
  auto const plan = buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.9f * uniformEstimate(order));
  EXPECT_EQ(plan.fine_elements, (std::vector<int>{1}));
  EXPECT_EQ(plan.element_steps[0].max_vp, kVp);
  EXPECT_EQ(plan.element_steps[1].max_vp, 2 * kVp);
}

// nvcc rejects extended device lambdas in the (private) TestBody of a gtest: free function.
template <class View>
void setOnDevice(View const& view, int i, float value) {
  Kokkos::parallel_for(
      "update speed", 1, KOKKOS_LAMBDA(int) { view(i) = value; });
}

TEST(SEM2DLocalTimeStepPlan, ReadsVelocitiesUpdatedOnDevice) {
  auto d = strip(2, {0, 100, 200}, true);
  model::ModelUnstruct<float, int> mesh(d);
  float const dt = 0.9f * uniformEstimate(2);
  EXPECT_TRUE(buildSEM2DLocalTimeStepPlan(mesh, kCfl, dt).fine_elements.empty());
  setOnDevice(d.model_vp_node_.view_device(), 4, 2 * kVp);
  d.model_vp_node_.modify_device();
  Kokkos::fence();
  auto const plan = buildSEM2DLocalTimeStepPlan(mesh, kCfl, dt);
  EXPECT_EQ(plan.fine_elements, (std::vector<int>{1}));
  EXPECT_EQ(plan.element_steps[1].max_vp, 2 * kVp);
}

TEST(SEM2DLocalTimeStepPlan, RejectsStepsRequiringMoreThanTwoLevels) {
  model::ModelUnstruct<float, int> mesh(strip(2, {0, 100, 125, 225}));
  EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.9f * uniformEstimate(2)), std::invalid_argument);
}

TEST(SEM2DLocalTimeStepPlan, RejectsInvalidParameters) {
  model::ModelUnstruct<float, int> mesh(strip(2, {0, 100}));
  for (float invalid : {0.f, -1.f, std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()}) {
    EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh, invalid, 0.001f), std::invalid_argument);
    EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh, kCfl, invalid), std::invalid_argument);
  }
  EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh, kCfl, std::numeric_limits<float>::denorm_min()),
               std::invalid_argument);
}

TEST(SEM2DLocalTimeStepPlan, RejectsInvalidVelocitiesForBothLocations) {
  for (bool nodal : {false, true})
    for (float invalid : {0.f, -1.f, std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()}) {
      auto d = strip(2, {0, 100}, nodal);
      (nodal ? d.model_vp_node_ : d.model_vp_element_).view_host()(0) = invalid;
      model::ModelUnstruct<float, int> mesh(d);
      EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.001f), std::invalid_argument);
    }
}

TEST(SEM2DLocalTimeStepPlan, RejectsInvalidCoordinatesAndConnectivity) {
  for (int kind = 0; kind < 3; ++kind) {
    auto d = strip(2, {0, 100});
    if (kind == 0) d.nodes_coords_x_.view_host()(0) = std::numeric_limits<float>::quiet_NaN();
    if (kind == 1) d.global_node_index_.view_host()(0, 0) = -1;
    if (kind == 2) d.global_node_index_.view_host()(0, 0) = d.n_node_;
    model::ModelUnstruct<float, int> mesh(d);
    EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.001f), std::invalid_argument);
  }
}

TEST(SEM2DLocalTimeStepPlan, RejectsZeroSpacing) {
  auto d = strip(2, {0, 100});
  d.nodes_coords_x_.view_host()(1) = d.nodes_coords_x_.view_host()(0);
  model::ModelUnstruct<float, int> mesh(d);
  EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh, kCfl, 0.001f), std::invalid_argument);
}

TEST(SEM2DLocalTimeStepPlan, RejectsNon2DAndEmptyMeshes) {
  auto d = strip(2, {0, 100});
  d.dim_ = 3;
  model::ModelUnstruct<float, int> mesh3d(d);
  EXPECT_THROW(buildSEM2DLocalTimeStepPlan(mesh3d, kCfl, 0.001f), std::invalid_argument);
  d.dim_ = 2;
  d.n_element_ = 0;
  model::ModelUnstruct<float, int> empty(d);
  EXPECT_THROW(buildSEM2DLocalTimeStepPlan(empty, kCfl, 0.001f), std::invalid_argument);
}

}  // namespace
}  // namespace test
}  // namespace solver

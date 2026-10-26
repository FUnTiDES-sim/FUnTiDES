/**
 * @file test_sem_solver_zdeformed.cc
 * @brief The z-deformed stiffness kernels must assemble the same forces as the generic ones.
 *
 * A Cartesian mesh stored as unstructured keeps its xi and eta edges parallel to x and y, so
 * the z-deformed kernels accept it, and so does the same mesh with its inner nodes shifted along
 * z, where the z-tilt terms of the inverse Jacobian no longer vanish. On both meshes each kernel
 * is compared with the generic kernel of the same physics on the same field. The team kernels
 * need one thread per element node and fall back to simpler paths below that.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "Integrals.h"
#include "cartesian_params.h"
#include "cartesian_unstruct_builder.h"
#include "data_type.h"
#include "device_view_access.h"
#include "model_unstruct.h"
#include "rhs_acoustic.h"
#include "rhs_elastic.h"
#include "sem_solver.h"
#include "sem_solver_data.h"
#include "sem_solver_impl.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace solver {
namespace fe {
namespace test {

namespace {

constexpr int kOrder = 1;
constexpr int kPointsPerElement = (kOrder + 1) * (kOrder + 1) * (kOrder + 1);
constexpr int kElementsX = 2;
constexpr int kElementsY = 2;
constexpr int kElementsZ = 3;
constexpr int kNumElements = kElementsX * kElementsY * kElementsZ;
constexpr float kSizeX = 200.0f;
constexpr float kSizeY = 160.0f;
constexpr float kSizeZ = 180.0f;

using Integral = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
using Mesh = model::ModelUnstruct<float, int>;
template <bool kOnNodes>
using AcousticSolver = SEMsolver<kOrder, Integral, Mesh, kOnNodes, utils::enums::physicType::kAcoustic>;
using ElasticSolver = SEMsolver<kOrder, Integral, Mesh, true, utils::enums::physicType::kElastic>;

dualVectorReal MakeUniform(int size, float value, const char* name) {
  dualVectorReal values = allocateVector<dualVectorReal>(size, name);
  Kokkos::deep_copy(values.view_host(), value);
  values.modify_host();
  return values;
}

/**
 * @brief Mesh of the same box whose inner node layers are shifted along z, each node by its own amount.
 *
 * The outer faces stay flat, so the domain is still a box, and the middle layer of elements has
 * both its bottom and top faces tilted.
 */
std::shared_ptr<model::ModelApi<float, int>> MakeDeformedMesh(bool on_nodes, bool elastic) {
  static_assert(kOrder == 1, "the nodes of the mesh below are its vertices");
  constexpr int kNodesX = kElementsX + 1;
  constexpr int kNodesY = kElementsY + 1;
  constexpr int kNodesZ = kElementsZ + 1;
  constexpr int kNumNodes = kNodesX * kNodesY * kNodesZ;
  constexpr float kStepX = kSizeX / kElementsX;
  constexpr float kStepY = kSizeY / kElementsY;
  constexpr float kStepZ = kSizeZ / kElementsZ;

  model::ModelUnstructData<float, int> data;
  data.order_ = kOrder;
  data.n_element_ = kNumElements;
  data.n_node_ = kNumNodes;
  data.lx_ = kSizeX;
  data.ly_ = kSizeY;
  data.lz_ = kSizeZ;
  data.ox_ = data.oy_ = data.oz_ = 0.0f;
  data.isModelOnNodes_ = on_nodes;
  data.isElastic_ = elastic;

  data.nodes_coords_x_ = allocateVector<dualVectorReal>(kNumNodes, "nodes coords x");
  data.nodes_coords_y_ = allocateVector<dualVectorReal>(kNumNodes, "nodes coords y");
  data.nodes_coords_z_ = allocateVector<dualVectorReal>(kNumNodes, "nodes coords z");
  for (int k = 0; k < kNodesZ; ++k)
    for (int j = 0; j < kNodesY; ++j)
      for (int i = 0; i < kNodesX; ++i) {
        int const node = i + j * kNodesX + k * kNodesX * kNodesY;
        bool const inner = k > 0 && k < kNodesZ - 1;
        float const shift = inner ? 0.2f * kStepZ * std::sin(1.9f * static_cast<float>(node) + 0.4f) : 0.0f;
        data.nodes_coords_x_.view_host()(node) = static_cast<float>(i) * kStepX;
        data.nodes_coords_y_.view_host()(node) = static_cast<float>(j) * kStepY;
        data.nodes_coords_z_.view_host()(node) = static_cast<float>(k) * kStepZ + shift;
      }
  data.nodes_coords_x_.modify_host();
  data.nodes_coords_y_.modify_host();
  data.nodes_coords_z_.modify_host();

  data.global_node_index_ = allocateArray2D<dualArrayInt>(kNumElements, kPointsPerElement, "global node index");
  for (int n = 0; n < kElementsZ; ++n)
    for (int m = 0; m < kElementsY; ++m)
      for (int l = 0; l < kElementsX; ++l) {
        int const e = l + m * kElementsX + n * kElementsX * kElementsY;
        for (int k = 0; k < 2; ++k)
          for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i)
              data.global_node_index_.view_host()(e, i + 2 * j + 4 * k) =
                  (l + i) + (m + j) * kNodesX + (n + k) * kNodesX * kNodesY;
      }
  data.global_node_index_.modify_host();

  // Same uniform material as the Cartesian builder.
  if (on_nodes) {
    data.model_rho_node_ = MakeUniform(kNumNodes, 1.0f, "model rho node");
    data.model_vp_node_ = MakeUniform(kNumNodes, 1500.0f, "model vp node");
    if (elastic) {
      data.model_vs_node_ = MakeUniform(kNumNodes, 755.0f, "model vs node");
      data.model_delta_node_ = MakeUniform(kNumNodes, 0.0f, "model delta node");
      data.model_epsilon_node_ = MakeUniform(kNumNodes, 0.0f, "model epsilon node");
      data.model_gamma_node_ = MakeUniform(kNumNodes, 0.0f, "model gamma node");
      data.model_theta_node_ = MakeUniform(kNumNodes, 0.0f, "model theta node");
      data.model_phi_node_ = MakeUniform(kNumNodes, 0.0f, "model phi node");
    }
  } else {
    data.model_rho_element_ = MakeUniform(kNumElements, 1.0f, "model rho elem");
    data.model_vp_element_ = MakeUniform(kNumElements, 1500.0f, "model vp elem");
  }

  auto mesh = std::make_shared<Mesh>(data);
  mesh->buildFaceConnectivity();
  return mesh;
}

/// Box with a different element size along each axis, Cartesian or deformed along z.
std::shared_ptr<model::ModelApi<float, int>> MakeMesh(bool on_nodes, bool elastic, bool deformed) {
  if (deformed) return MakeDeformedMesh(on_nodes, elastic);
  model::CartesianParams<float, int> const params(kOrder, kElementsX, kElementsY, kElementsZ, kSizeX, kSizeY, kSizeZ,
                                                  on_nodes, elastic);
  return model::CartesianUnstructBuilder<float, int>(params).getModel(false);
}

vectorReal MakeField(int n_nodes, float phase, const char* name) {
  vectorReal field = allocateVector<vectorReal>(n_nodes, name);
  test_utils::fillFromHost(field, [&](const auto& host) {
    for (int i = 0; i < n_nodes; ++i) host(i) = std::sin(0.7f * static_cast<float>(i) + phase);
  });
  return field;
}

/// Forces of every component after @p kernel, starting from zeroed work vectors.
template <typename SolverType, typename Kernel>
std::vector<float> Assemble(SolverType& solver, int n_nodes, Kernel&& kernel) {
  solver.resetGlobalVectors(n_nodes);
  kernel();
  Kokkos::fence();
  std::vector<float> forces;
  for (int f = 0; f < SolverType::kNumFields; ++f) {
    auto const force = test_utils::toHost(solver.getForceVector(f));
    for (int i = 0; i < n_nodes; ++i) forces.push_back(force(i));
  }
  return forces;
}

void ExpectSameForces(std::vector<float> const& ref, std::vector<float> const& got, const char* what) {
  ASSERT_EQ(ref.size(), got.size());
  float scale = 0.0f;
  for (float const v : ref) scale = std::max(scale, std::fabs(v));
  ASSERT_GT(scale, 0.0f);
  for (size_t i = 0; i < ref.size(); ++i) EXPECT_NEAR(got[i], ref[i], 1e-4f * scale) << what << ", entry " << i;
}

template <bool kOnNodes>
void CheckAcoustic(bool deformed) {
  auto mesh = MakeMesh(kOnNodes, false, deformed);
  AcousticSolver<kOnNodes> solver;
  solver.computeFEInit(*mesh, {0.0f, 0.0f, 0.0f}, false, 0.0f);
  int const n_nodes = mesh->getNumberOfNodes();

  vectorReal p_prev = MakeField(n_nodes, 0.3f, "pPrev");
  vectorReal p_curr = MakeField(n_nodes, 1.1f, "pCurr");
  arrayReal rhs_term = allocateArray2D<arrayReal>(1, 1, "rhsTerm");
  vectorInt rhs_elem = allocateVector<vectorInt>(1, "rhsElem");
  arrayReal rhs_weights = allocateArray2D<arrayReal>(1, kPointsPerElement, "rhsWeights");
  Kokkos::deep_copy(rhs_term, 0.0f);
  Kokkos::deep_copy(rhs_elem, 0);
  Kokkos::deep_copy(rhs_weights, 0.0f);

  WavefieldAcoustic wf(p_prev, p_curr);
  RhsAcoustic rhs(rhs_term, rhs_elem, rhs_weights);
  typename AcousticSolver<kOnNodes>::DataType data(wf, rhs);

  solver.prepareZDeformedGeometry();
  auto const ref = Assemble(solver, n_nodes, [&] { solver.computeElementContributions_Acoustic_Flat(data); });

  ExpectSameForces(ref, Assemble(solver, n_nodes, [&] { solver.computeElementContributions_Acoustic_FlatZ(data); }),
                   "FlatZ");
  ExpectSameForces(ref, Assemble(solver, n_nodes, [&] { solver.computeElementContributions_Acoustic_TeamZ(data); }),
                   "TeamZ");

  // Every element once, in another order.
  vectorInt reversed = allocateVector<vectorInt>(kNumElements, "reversed");
  test_utils::fillFromHost(reversed, [&](const auto& host) {
    for (int i = 0; i < kNumElements; ++i) host(i) = kNumElements - 1 - i;
  });
  ExpectSameForces(
      ref, Assemble(solver, n_nodes, [&] { solver.computeElementContributionsFromList(data, reversed, kNumElements); }),
      "reversed list");

  // Two lists of the same size: the node table of the first must not be reused for the second.
  constexpr int kHalf = kNumElements / 2;
  vectorInt first = allocateVector<vectorInt>(kHalf, "first");
  vectorInt second = allocateVector<vectorInt>(kHalf, "second");
  test_utils::fillFromHost(first, [&](const auto& host) {
    for (int i = 0; i < kHalf; ++i) host(i) = i;
  });
  test_utils::fillFromHost(second, [&](const auto& host) {
    for (int i = 0; i < kHalf; ++i) host(i) = kHalf + i;
  });
  ExpectSameForces(ref,
                   Assemble(solver, n_nodes,
                            [&] {
                              solver.computeElementContributionsFromList(data, first, kHalf);
                              solver.computeElementContributionsFromList(data, second, kHalf);
                            }),
                   "two half lists");
}

void CheckTti(bool deformed) {
  auto mesh = MakeMesh(true, true, deformed);
  ElasticSolver solver;
  solver.setAnisotropyType(model::AnisotropyType::kTTI);
  solver.computeFEInit(*mesh, {0.0f, 0.0f, 0.0f}, false, 0.0f);
  int const n_nodes = mesh->getNumberOfNodes();

  vectorReal ux_prev = MakeField(n_nodes, 0.1f, "uxPrev");
  vectorReal ux_curr = MakeField(n_nodes, 0.5f, "uxCurr");
  vectorReal uy_prev = MakeField(n_nodes, 0.9f, "uyPrev");
  vectorReal uy_curr = MakeField(n_nodes, 1.3f, "uyCurr");
  vectorReal uz_prev = MakeField(n_nodes, 1.7f, "uzPrev");
  vectorReal uz_curr = MakeField(n_nodes, 2.1f, "uzCurr");
  arrayReal rhs_term_x = allocateArray2D<arrayReal>(1, 1, "rhsTermX");
  arrayReal rhs_term_y = allocateArray2D<arrayReal>(1, 1, "rhsTermY");
  arrayReal rhs_term_z = allocateArray2D<arrayReal>(1, 1, "rhsTermZ");
  vectorInt rhs_elem = allocateVector<vectorInt>(1, "rhsElem");
  arrayReal rhs_weights = allocateArray2D<arrayReal>(1, kPointsPerElement, "rhsWeights");
  Kokkos::deep_copy(rhs_term_x, 0.0f);
  Kokkos::deep_copy(rhs_term_y, 0.0f);
  Kokkos::deep_copy(rhs_term_z, 0.0f);
  Kokkos::deep_copy(rhs_elem, 0);
  Kokkos::deep_copy(rhs_weights, 0.0f);

  WavefieldElastic wf(ux_prev, ux_curr, uy_prev, uy_curr, uz_prev, uz_curr);
  RhsElastic rhs(rhs_term_x, rhs_term_y, rhs_term_z, rhs_elem, rhs_weights);
  ElasticSolver::DataType data(wf, rhs);

  solver.prepareZDeformedGeometry();
  auto const ref = Assemble(solver, n_nodes, [&] { solver.computeElementContributions_Tti_Flat(data); });
  ExpectSameForces(ref, Assemble(solver, n_nodes, [&] { solver.computeElementContributions_Tti_TeamZ(data); }),
                   "Tti TeamZ");
}

}  // namespace

TEST(SemSolverZDeformed, AcousticOnElementsMatchesGenericKernel) { CheckAcoustic<false>(false); }

TEST(SemSolverZDeformed, AcousticOnNodesMatchesGenericKernel) { CheckAcoustic<true>(false); }

TEST(SemSolverZDeformed, TtiOnNodesMatchesGenericKernel) { CheckTti(false); }

TEST(SemSolverZDeformed, AcousticOnElementsMatchesGenericKernelOnDeformedMesh) { CheckAcoustic<false>(true); }

TEST(SemSolverZDeformed, AcousticOnNodesMatchesGenericKernelOnDeformedMesh) { CheckAcoustic<true>(true); }

TEST(SemSolverZDeformed, TtiOnNodesMatchesGenericKernelOnDeformedMesh) { CheckTti(true); }

}  // namespace test
}  // namespace fe
}  // namespace solver

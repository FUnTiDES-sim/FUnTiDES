#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_

#include <data_type.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>

#include "Integrals.h"
#include "elastic_flux.h"
#include "mesh_type_traits.h"
#include "sem_solver.h"

namespace solver {
namespace fe {
namespace detail {
/// @brief True when T declares a nested type `TeamGemm`, i.e. the integral back-end provides the team GEMM kernels.
template <typename, typename = void>
struct has_team_gemm : std::false_type {};
template <typename T>
struct has_team_gemm<T, std::void_t<typename T::TeamGemm>> : std::true_type {};

/// @brief True when T declares a nested type `ZDeformedSumFact`, i.e. the integral back-end provides
/// the z-deformed sum-factorization kernels.
template <typename, typename = void>
struct has_z_deformed_sumfact : std::false_type {};
template <typename T>
struct has_z_deformed_sumfact<T, std::void_t<typename T::ZDeformedSumFact>> : std::true_type {};
}  // namespace detail


template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::prepareZDeformedGeometry() {
  if (zDeformedReady_) return;
  zDeformedReady_ = true;
  zDeformedEnabled_ = false;
  zDeformedGeom_ = vectorReal();
  zDeformedNodes_ = vectorInt();
  zDeformedNodesList_ = vectorInt();
  zDeformedInvRhoNodes_ = vectorReal();
  zDeformedBasisTab_ = vectorReal();

  bool const is_acoustic = (PHYSICS == utils::enums::physicType::kAcoustic);

  auto mesh_pc = m_mesh;
  int const nElems = mesh_pc.getNumberOfElements();
  vectorReal geom = allocateVector<vectorReal>(nElems * kZGeomStride, "zDeformedGeom");
  float* geomPtr = geom.data();

  int nonConforming = 0;
  Kokkos::parallel_reduce(
      "ZDeformed Geometry", Kokkos::RangePolicy<>(0, nElems),
      KOKKOS_LAMBDA(const int e, int& bad) {
        float X[8][3];
        auto const eIdx = mesh_pc.elementIndex(e);
        int I = 0;
        for (int kv = 0; kv < 2; ++kv)
          for (int jv = 0; jv < 2; ++jv)
            for (int iv = 0; iv < 2; ++iv) mesh_pc.vertexCoords(mesh_pc.globalVertexIndex(eIdx, iv, jv, kv), X[I++]);

        // Vertex k has indices (k % 2, (k % 4) / 2, k / 4): x may follow the first only, y the second only.
        float const hx = X[1][0] - X[0][0];
        float const hy = X[2][1] - X[0][1];
        float const tolX = 1e-4f * Kokkos::fabs(hx);
        float const tolY = 1e-4f * Kokkos::fabs(hy);
        for (int k = 0; k < 8; ++k) {
          if (Kokkos::fabs(X[k][0] - X[k % 2][0]) > tolX || Kokkos::fabs(X[k][1] - X[2 * ((k % 4) / 2)][1]) > tolY) {
            ++bad;
            break;
          }
        }

        geomPtr[e] = 0.5f * hx;
        geomPtr[nElems + e] = 0.5f * hy;
        for (int k = 0; k < 8; ++k) geomPtr[(2 + k) * nElems + e] = X[k][2];

        // With the model on nodes, 1/rho comes from zDeformedInvRhoNodeTable().
        if (is_acoustic && !IS_MODEL_ON_NODES)
          geomPtr[(kZGeomStride - 1) * nElems + e] = 1.0f / mesh_pc.getModelRhoOnElement(e);
        else
          geomPtr[(kZGeomStride - 1) * nElems + e] = 0.0f;
      },
      nonConforming);

  if (nonConforming > 0) return;

  zDeformedGeom_ = geom;
  zDeformedEnabled_ = true;
  if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    vectorReal basisTab = allocateVector<vectorReal>(INTEGRAL_TYPE::kBasisTableSize, "zDeformedBasisTab");
    auto h_basisTab = Kokkos::create_mirror_view(basisTab);
    for (int idx = 0; idx < INTEGRAL_TYPE::kBasisTableSize; ++idx)
      h_basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx);
    Kokkos::deep_copy(basisTab, h_basisTab);
    zDeformedBasisTab_ = basisTab;
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
int const* SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::zDeformedNodeTable(
    bool element_major) {
  bool const list_on = m_list_mode_;
  int const n_visited = list_on ? m_n_elem_list_ : m_mesh.getNumberOfElements();
  int const* list_key = list_on ? m_elem_list_.data() : nullptr;
  if (zDeformedNodes_.extent(0) == static_cast<size_t>(n_visited) * kPointsPerElement &&
      zDeformedNodesList_.data() == list_key && zDeformedNodesElementMajor_ == element_major)
    return zDeformedNodes_.data();

  auto mesh_pc = m_mesh;
  auto list_local = m_elem_list_;
  // Filled on the device, so that the kernels never read the mesh connectivity, whose pages the
  // host may pull back under unified memory.
  vectorInt nodes = allocateVector<vectorInt>(n_visited * kPointsPerElement, "zDeformedNodes");
  int* nodesPtr = nodes.data();
  Kokkos::parallel_for(
      "ZDeformed Node Table", Kokkos::RangePolicy<>(0, n_visited), KOKKOS_LAMBDA(const int i) {
        constexpr int dim = ORDER + 1;
        int const e = list_on ? list_local[i] : i;
        // The team kernels read the indices element by element, the flat kernel node by node.
        for (int q = 0; q < kPointsPerElement; ++q)
          nodesPtr[element_major ? i * kPointsPerElement + q : q * n_visited + i] =
              mesh_pc.globalNodeIndex(e, q % dim, (q / dim) % dim, q / (dim * dim));
      });

  zDeformedNodes_ = nodes;
  zDeformedNodesList_ = list_on ? m_elem_list_ : vectorInt();
  zDeformedNodesElementMajor_ = element_major;
  return zDeformedNodes_.data();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
float const* SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::zDeformedInvRhoNodeTable() {
  if constexpr (IS_MODEL_ON_NODES) {
    int const n_nodes = m_mesh.getNumberOfNodes();
    if (zDeformedInvRhoNodes_.extent(0) == static_cast<size_t>(n_nodes)) return zDeformedInvRhoNodes_.data();
    auto mesh_pc = m_mesh;
    vectorReal table = allocateVector<vectorReal>(n_nodes, "zDeformedInvRhoNodes");
    float* tablePtr = table.data();
    Kokkos::parallel_for(
        "ZDeformed InvRho Node Table", Kokkos::RangePolicy<>(0, n_nodes),
        KOKKOS_LAMBDA(const int g) { tablePtr[g] = 1.0f / mesh_pc.getModelRhoOnNodes(g); });
    zDeformedInvRhoNodes_ = table;
  }
  return zDeformedInvRhoNodes_.data();
}


template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributions(
    const DataType& data) {
  if (!attenuationEnabled_ || nSls_ <= 0) return;

  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic)
    computeAttenuationContributionsAcoustic(data);
  else
    computeAttenuationContributionsElastic(data);
}



}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_

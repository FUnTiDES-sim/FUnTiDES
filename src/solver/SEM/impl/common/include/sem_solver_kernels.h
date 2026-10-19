#ifndef FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_
#define FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_

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




// Builds the TTI tensor CTTI = M * CVTI * M^T, with CVTI the VTI tensor (symmetry axis z) in Voigt storage
// and M the rotation of the symmetry axis given by theta and phi.
// theta and phi are converted from degrees here.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
template <physicType P, typename>
PROXY_HOST_DEVICE void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeCMatrix(
    float const vp, float const vs, float const rho, float const delta, float const epsilon, float const gamma,
    float const phi, float const theta, float (&CTTI)[6][6]) {
  const float rho_vp2 = rho * vp * vp;
  const float rho_vs2 = rho * vs * vs;
  const float two_eps = 2.0f * epsilon;
  const float two_gam = 2.0f * gamma;

  float CVTI[6][6] = {0.0f};
  CVTI[0][0] = rho_vp2 * (1.0f + two_eps);
  CVTI[1][1] = CVTI[0][0];
  CVTI[2][2] = rho_vp2;
  CVTI[3][3] = rho_vs2;
  CVTI[4][4] = CVTI[3][3];
  CVTI[5][5] = rho_vs2 * (1.0f + two_gam);
  CVTI[0][1] = CVTI[0][0] - 2.0f * CVTI[5][5];
  CVTI[1][0] = CVTI[0][1];

  const float vp2_vs2 = vp * vp - vs * vs;
  const float sqrt_arg = vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2;
  CVTI[0][2] = rho * sqrtf(sqrt_arg) - rho_vs2;
  CVTI[1][2] = CVTI[0][2];
  CVTI[2][0] = CVTI[0][2];
  CVTI[2][1] = CVTI[0][2];

  constexpr float DEG_TO_RAD = 3.14159265358979323846f / 180.0f;
  const float theta_rad = theta * DEG_TO_RAD;
  const float phi_rad = phi * DEG_TO_RAD;

  const float ctheta = cosf(theta_rad);
  const float stheta = sinf(theta_rad);
  const float cphi = cosf(phi_rad);
  const float sphi = sinf(phi_rad);

  const float ct_cp = ctheta * cphi;
  const float ct_sp = ctheta * sphi;
  const float st_cp = stheta * cphi;
  const float st_sp = stheta * sphi;

  float R[3][3];
  // Rotates the VTI axis to (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)).
  R[0][0] = ct_cp;
  R[0][1] = -sphi;
  R[0][2] = st_cp;
  R[1][0] = ct_sp;
  R[1][1] = cphi;
  R[1][2] = st_sp;
  R[2][0] = -stheta;
  R[2][1] = 0.0f;
  R[2][2] = ctheta;

  const float R00_2 = R[0][0] * R[0][0];
  const float R01_2 = R[0][1] * R[0][1];
  const float R02_2 = R[0][2] * R[0][2];
  const float R10_2 = R[1][0] * R[1][0];
  const float R11_2 = R[1][1] * R[1][1];
  const float R12_2 = R[1][2] * R[1][2];
  const float R20_2 = R[2][0] * R[2][0];
  const float R21_2 = R[2][1] * R[2][1];
  const float R22_2 = R[2][2] * R[2][2];

  float M[6][6] = {0.0f};
  M[0][0] = R00_2;
  M[0][1] = R01_2;
  M[0][2] = R02_2;
  M[1][0] = R10_2;
  M[1][1] = R11_2;
  M[1][2] = R12_2;
  M[2][0] = R20_2;
  M[2][1] = R21_2;
  M[2][2] = R22_2;
  // M * CVTI * M^T needs the Bond matrix for stress, whose factors of two sit in
  // the upper right block; the lower left placement builds the one for strain.
  M[0][3] = 2.0f * R[0][1] * R[0][2];
  M[0][4] = 2.0f * R[0][0] * R[0][2];
  M[0][5] = 2.0f * R[0][0] * R[0][1];
  M[1][3] = 2.0f * R[1][1] * R[1][2];
  M[1][4] = 2.0f * R[1][0] * R[1][2];
  M[1][5] = 2.0f * R[1][0] * R[1][1];
  M[2][3] = 2.0f * R[2][1] * R[2][2];
  M[2][4] = 2.0f * R[2][0] * R[2][2];
  M[2][5] = 2.0f * R[2][0] * R[2][1];
  M[3][0] = R[1][0] * R[2][0];
  M[3][1] = R[1][1] * R[2][1];
  M[3][2] = R[1][2] * R[2][2];
  M[3][3] = R[1][1] * R[2][2] + R[1][2] * R[2][1];
  M[3][4] = R[1][0] * R[2][2] + R[1][2] * R[2][0];
  M[3][5] = R[1][0] * R[2][1] + R[1][1] * R[2][0];
  M[4][0] = R[0][0] * R[2][0];
  M[4][1] = R[0][1] * R[2][1];
  M[4][2] = R[0][2] * R[2][2];
  M[4][3] = R[0][1] * R[2][2] + R[0][2] * R[2][1];
  M[4][4] = R[0][0] * R[2][2] + R[0][2] * R[2][0];
  M[4][5] = R[0][0] * R[2][1] + R[0][1] * R[2][0];
  M[5][0] = R[0][0] * R[1][0];
  M[5][1] = R[0][1] * R[1][1];
  M[5][2] = R[0][2] * R[1][2];
  M[5][3] = R[0][1] * R[1][2] + R[0][2] * R[1][1];
  M[5][4] = R[0][0] * R[1][2] + R[0][2] * R[1][0];
  M[5][5] = R[0][0] * R[1][1] + R[0][1] * R[1][0];

  float temp[6][6];

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++) {
        sum += M[i][k] * CVTI[k][j];
      }
      temp[i][j] = sum;
    }
  }

  for (int i = 0; i < 6; i++) {
    for (int j = i; j < 6; j++) {
      float sum = 0.0f;
      for (int k = 0; k < 6; k++) {
        sum += temp[i][k] * M[j][k];
      }
      CTTI[i][j] = sum;
      if (i != j) CTTI[j][i] = sum;
    }
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
template <physicType P, typename>
PROXY_HOST_DEVICE void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeTtiCompact(
    float const vp, float const vs, float const rho, float const delta, float const epsilon, float const gamma,
    float const phi, float const theta, float (&p)[flux::kTtiCompactSize]) {
  // Same VTI coefficients as computeCMatrix().
  const float rho_vp2 = rho * vp * vp;
  const float rho_vs2 = rho * vs * vs;
  const float c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
  const float c66 = rho_vs2 * (1.0f + 2.0f * gamma);
  const float vp2_vs2 = vp * vp - vs * vs;
  const float c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;

  constexpr float DEG_TO_RAD = 3.14159265358979323846f / 180.0f;
  const float theta_rad = theta * DEG_TO_RAD;
  const float phi_rad = phi * DEG_TO_RAD;
  const float stheta = sinf(theta_rad);
  float const n[3] = {stheta * cosf(phi_rad), stheta * sinf(phi_rad), cosf(theta_rad)};
  flux::ttiCompactFromVti(c11, c13, rho_vp2, rho_vs2, c66, n, p);
}

}  // namespace solver

#endif  // FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_

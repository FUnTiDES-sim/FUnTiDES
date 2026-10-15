#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_ASSEMBLY_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_ASSEMBLY_H_

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
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrix() {
  auto mesh_local = m_mesh;
  auto local_massMatrixGlobal = massMatrixGlobal_;
  bool const mask_enabled = m_mask_enabled_;
  auto element_mask = m_element_mask_;
  int const mask_active_value = m_mask_active_value_;

  Kokkos::parallel_for(
      "Solver Compute GMatrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (mask_enabled && element_mask[elementNumber] != mask_active_value) return;
        float massMatrixLocal[kPointsPerElement] = {0};
        int const dim = mesh_local.getOrder() + 1;

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        INTEGRAL_TYPE::computeMassTerm(cornerCoords, [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

        // Acoustic: 1 / (vp^2 rho). Elastic: rho.
        real_t model_factor = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            model_factor =
                1.0f / (mesh_local.getModelVpOnElement(elementNumber) * mesh_local.getModelVpOnElement(elementNumber) *
                        mesh_local.getModelRhoOnElement(elementNumber));
          } else {
            model_factor = mesh_local.getModelRhoOnElement(elementNumber);
          }
        }

        for (int i = 0; i < mesh_local.getNumberOfPointsPerElement(); ++i) {
          int x = i % dim;
          int z = (i / dim) % dim;
          int y = i / (dim * dim);
          int const gIndex = mesh_local.globalNodeIndex(elementNumber, x, y, z);

          if constexpr (IS_MODEL_ON_NODES) {
            if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
              model_factor = 1.0f / (mesh_local.getModelVpOnNodes(gIndex) * mesh_local.getModelVpOnNodes(gIndex) *
                                     mesh_local.getModelRhoOnNodes(gIndex));
            } else {
              model_factor = mesh_local.getModelRhoOnNodes(gIndex);
            }
          }

          massMatrixLocal[i] *= model_factor;
          ATOMICADD(local_massMatrixGlobal[gIndex], massMatrixLocal[i]);
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrix() {
  auto mesh_local = m_mesh;
  bool const mask_enabled = m_mask_enabled_;
  auto element_mask = m_element_mask_;
  int const mask_active_value = m_mask_active_value_;

  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>, kNumFields> local_dampingMatrixGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_dampingMatrixGlobal[f] = dampingMatrixGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Compute Damping Matrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        (void)local_dampingMatrixGlobal;
        if (mask_enabled && element_mask[elementNumber] != mask_active_value) return;
        for (int i = 0; i < 6; ++i) {
          int f = mesh_local.getGlobalFace(elementNumber, static_cast<model::CubicFace>(i));

          // Only boundary faces carry an absorbing term.
          if (!mesh_local.isBoundaryFace(f)) continue;

          float coords[4][3];
          for (int j = 0; j < 4; ++j) {
            int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) {
              coords[j][d] = mesh_local.nodeCoord(globalNodeIndex, d);
            }
          }

          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            real_t model_rho = 0.0f;
            real_t model_vp = 0.0f;
            real_t alpha = 0.0f;

            if constexpr (!IS_MODEL_ON_NODES) {
              model_rho = mesh_local.getModelRhoOnElement(elementNumber);
              model_vp = mesh_local.getModelVpOnElement(elementNumber);
              alpha = 1.0 / (model_rho * model_vp);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Free-surface nodes are not damped.
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                model_rho = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                model_vp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                alpha = 1.0 / (model_rho * model_vp);
              }

              real_t localIncrement = alpha * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrement);
            }
          } else {
            float normal[3];
            mesh_local.faceNormal(elementNumber, static_cast<model::CubicFace>(i), normal);
            real_t nx = normal[0], ny = normal[1], nz = normal[2];

            real_t density, velocityVp, velocityVs;

            if constexpr (!IS_MODEL_ON_NODES) {
              density = mesh_local.getModelRhoOnElement(elementNumber);
              velocityVp = mesh_local.getModelVpOnElement(elementNumber);
              velocityVs = mesh_local.getModelVsOnElement(elementNumber);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Free-surface nodes are not damped.
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                density = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                velocityVp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                velocityVs = mesh_local.getModelVsOnNodes(globalNodeIndex);
              }

              real_t aux = density * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              real_t localIncrementx = aux * (velocityVp * std::fabs(nx) + velocityVs * std::sqrt(ny * ny + nz * nz));
              real_t localIncrementy = aux * (velocityVp * std::fabs(ny) + velocityVs * std::sqrt(nx * nx + nz * nz));
              real_t localIncrementz = aux * (velocityVp * std::fabs(nz) + velocityVs * std::sqrt(nx * nx + ny * ny));

              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrementx);
              ATOMICADD(local_dampingMatrixGlobal[1][globalNodeIndex], localIncrementy);
              ATOMICADD(local_dampingMatrixGlobal[2][globalNodeIndex], localIncrementz);
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::allocateFEarrays() {
  massMatrixGlobal_ = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), "massMatrixGlobal");

  static constexpr const char* dampingNames[3] = {"dampingX", "dampingY", "dampingZ"};
  for (int f = 0; f < kNumFields; ++f) {
    dampingMatrixGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), dampingNames[f]);
  }

  static constexpr const char* workVectorNames[3] = {"workVec0", "workVec1", "workVec2"};
  for (int f = 0; f < kNumFields; ++f) {
    workVectorsGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), workVectorNames[f]);
  }

  if (attenuationEnabled_ && nSls_ > 0) {
    static constexpr const char* attWorkNames[3] = {"attWorkVec0", "attWorkVec1", "attWorkVec2"};
    static constexpr const char* attMemNames[3] = {"attMemory0", "attMemory1", "attMemory2"};
    for (int f = 0; f < kNumFields; ++f) {
      attenuationWorkVectorsGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), attWorkNames[f]);
      attenuationMemoryVariables_[f] = allocateArray2D<arrayReal>(m_mesh.getNumberOfNodes(), nSls_, attMemNames[f]);
    }
  }

  spongeTaperCoeff_ = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), "spongeTaperCoeff");
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initFEarrays() {
  initSpongeValues();

  if (attenuationEnabled_ && nSls_ > 0) {
    for (int n = 0; n < m_mesh.getNumberOfNodes(); ++n) {
      for (int f = 0; f < kNumFields; ++f) {
        attenuationWorkVectorsGlobal_[f](n) = 0.0f;
        for (int l = 0; l < nSls_; ++l) {
          attenuationMemoryVariables_[f](n, l) = 0.0f;
        }
      }
    }
    FENCE
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initSpongeValues() {
  // @todo VERIFY: unit and origin of the peak damping factor sigma_max = 0.15 (dimensionless?).
  const double sigma_max = 0.15;

  for (int n = 0; n < m_mesh.getNumberOfNodes(); n++) {
    const double x = m_mesh.nodeCoord(n, 0);
    const double y = m_mesh.nodeCoord(n, 1);
    const double z = m_mesh.nodeCoord(n, 2);

    const double distToFrontierX = (surface_sponge_) ? m_mesh.domainSize(0) - x : std::min(m_mesh.domainSize(0) - x, x);
    const double distToFrontierY = std::min(m_mesh.domainSize(1) - y, y);
    const double distToFrontierZ = std::min(m_mesh.domainSize(2) - z, z);

    double minDistToFrontier = std::max(m_mesh.domainSize(0), std::max(m_mesh.domainSize(1), m_mesh.domainSize(2)));

    bool is_sponge = false;
    if (distToFrontierX < sponge_size_[0]) {
      is_sponge = true;
      minDistToFrontier = std::min(minDistToFrontier, distToFrontierX);
    }
    if (distToFrontierY < sponge_size_[1]) {
      is_sponge = true;
      minDistToFrontier = std::min(minDistToFrontier, distToFrontierY);
    }
    if (distToFrontierZ < sponge_size_[2]) {
      is_sponge = true;
      minDistToFrontier = std::min(minDistToFrontier, distToFrontierZ);
    }

    // Gaussian taper in the distance to the closest sponge frontier.
    if (is_sponge) {
      double d = minDistToFrontier;
      double delta = taper_delta_;
      double sigma = sigma_max * std::exp(-((d / delta) * (d / delta)));
      spongeTaperCoeff_(n) = 1.0 / (1.0 + sigma);
    } else {
      spongeTaperCoeff_(n) = 1.0;
    }
  }

  FENCE
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const vectorReal& fieldGlobal, const char* fieldName) {
  std::cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
            << " after computeOneStep = " << fieldGlobal(m_mesh.globalNodeIndex(e, 0, 0, 0)) << std::endl;
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
  R[0][0] = ct_cp;
  R[0][1] = ct_sp;
  R[0][2] = -stheta;
  R[1][0] = -sphi;
  R[1][1] = cphi;
  R[1][2] = 0.0f;
  R[2][0] = st_cp;
  R[2][1] = st_sp;
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
    [[maybe_unused]] float const phi, float const theta, float (&p)[flux::kTtiCompactSize]) {
  // Same VTI coefficients as computeCMatrix().
  const float rho_vp2 = rho * vp * vp;
  const float rho_vs2 = rho * vs * vs;
  const float c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
  const float c66 = rho_vs2 * (1.0f + 2.0f * gamma);
  const float vp2_vs2 = vp * vp - vs * vs;
  const float c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;

  constexpr float DEG_TO_RAD = 3.14159265358979323846f / 180.0f;
  const float theta_rad = theta * DEG_TO_RAD;
  float const n[3] = {-sinf(theta_rad), 0.0f, cosf(theta_rad)};
  flux::ttiCompactFromVti(c11, c13, rho_vp2, rho_vs2, c66, n, p);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrixMasked(
    const vectorInt& elem_mask, int active_value) {
  Kokkos::deep_copy(massMatrixGlobal_, 0.0f);
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeGlobalMassMatrix();
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrixMasked(
    const vectorInt& elem_mask, int active_value) {
  for (int f = 0; f < kNumFields; ++f) Kokkos::deep_copy(dampingMatrixGlobal_[f], 0.0f);
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeDampingMatrix();
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsMasked(
    const DataType& data, const vectorInt& elem_mask, int active_value) {
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeElementContributions(data);
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsFromList(
    const DataType& data, const vectorInt& elem_list, int n_elems) {
  m_elem_list_ = elem_list;
  m_n_elem_list_ = n_elems;
  m_list_mode_ = true;
  computeElementContributions(data);
  m_list_mode_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListForward(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_node_list_ = node_list;
  m_n_node_list_ = n_nodes;
  m_node_list_mode_ = true;
  updateFieldsForward(dt, data);
  m_node_list_mode_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListForwardAndReset(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_reset_forces_in_update_ = true;
  updateFieldsFromListForward(dt, data, node_list, n_nodes);
  m_reset_forces_in_update_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListBackward(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_node_list_ = node_list;
  m_n_node_list_ = n_nodes;
  m_node_list_mode_ = true;
  updateFieldsBackward(dt, data);
  m_node_list_mode_ = false;
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_ASSEMBLY_H_
